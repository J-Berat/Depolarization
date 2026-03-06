using CairoMakie
using Random
using Statistics
using Printf

include(joinpath(@__DIR__, "..", "lib", "DepolLib.jl"))
using .DepolLib

"""
    first_derivative(...)

    1D centered derivative.
"""
function first_derivative(f::AbstractVector, Δ::Real)
    n = length(f)
    n >= 3 || error("Need at least 3 points for derivative")
    out = similar(f, n - 2)
    inv2Δ = 1 / (2Δ)
    @inbounds for i in 2:(n-1)
        out[i-1] = (f[i+1] - f[i-1]) * inv2Δ
    end
    return out
end

"""
    count_reversals_in_window(...)

    Counts reversals in a window.
"""
count_reversals_in_window(ks::Vector{Int}, kL::Int, kR::Int) = count(k -> (k >= kL && k <= kR-1), ks)

"""
    tight_bounds_for_reversal(...)

    Estimates local bounds around a reversal.
"""
function tight_bounds_for_reversal(B::Vector{Float64}, Δ::Real, k0::Int, ks_all::Vector{Int};
                                   deriv_tol::Real=0.05,
                                   stop_at_next_reversal::Bool=true,
                                   max_half_width_pix::Int=18,
                                   isolate_fallback::Bool=true)
    n = length(B)
    d1 = first_derivative(B, Δ)
    i0 = clamp(k0 - 1, 1, length(d1))

    """
        includes_other_reversal(...)
    
        Returns `true` if another reversal index (different from `k0`) lies inside `[kL, kR-1]`.
    """
    function includes_other_reversal(kL::Int, kR::Int)
        for k in ks_all
            if k != k0 && (k >= kL && k <= kR - 1)
                return true
            end
        end
        return false
    end

    iL = i0
    while iL > 1
        if (i0 - (iL - 1)) > max_half_width_pix
            break
        end
        abs(d1[iL]) < deriv_tol && break
        if sign_eps(d1[iL-1]) != 0 && sign_eps(d1[iL]) != 0 && sign_eps(d1[iL-1]) != sign_eps(d1[iL])
            break
        end
        if stop_at_next_reversal
            k_left_candidate = iL
            k_right_current = i0 + 1
            includes_other_reversal(k_left_candidate, k_right_current) && break
        end
        iL -= 1
    end

    iR = i0
    while iR < length(d1)
        if ((iR + 1) - i0) > max_half_width_pix
            break
        end
        abs(d1[iR]) < deriv_tol && break
        if sign_eps(d1[iR]) != 0 && sign_eps(d1[iR+1]) != 0 && sign_eps(d1[iR]) != sign_eps(d1[iR+1])
            break
        end
        if stop_at_next_reversal
            k_left_current = i0 + 1
            k_right_candidate = iR + 2
            includes_other_reversal(k_left_current, k_right_candidate) && break
        end
        iR += 1
    end

    kL = clamp(iL + 1, 1, n)
    kR = clamp(iR + 1, 1, n)

    if isolate_fallback && count_reversals_in_window(ks_all, kL, kR) > 1
        best = nothing
        for w in 0:(n - 1)
            a = max(1, k0 - w)
            b = min(n, (k0 + 1) + w)
            if count_reversals_in_window(ks_all, a, b) == 1
                best = (a, b)
                break
            end
        end
        if best === nothing
            kL, kR = max(1, k0), min(n, k0 + 1)
        else
            kL, kR = best
        end
    end

    return kL, kR
end

"""
    merge_intervals_local(...)

    Merges adjacent or overlapping integer intervals.
"""
function merge_intervals_local(intervals::Vector{Tuple{Int,Int}})
    isempty(intervals) && return Tuple{Int,Int}[]
    ints = sort(intervals, by=first)
    merged = Tuple{Int,Int}[]
    a, b = ints[1]
    for (c, d) in ints[2:end]
        if c <= b + 1
            b = max(b, d)
        else
            push!(merged, (a, b))
            a, b = c, d
        end
    end
    push!(merged, (a, b))
    return merged
end

"""
    extract_fdf_spectrum(...)

    Extracts per-pixel FDF spectrum while auto-detecting the `phi` axis.
"""
function extract_fdf_spectrum(FDFcube, i::Int, j::Int, nphi::Int)
    sz = size(FDFcube)
    dims_phi = findall(d -> sz[d] == nphi, 1:ndims(FDFcube))
    isempty(dims_phi) && error("No FDF axis matches nphi=$nphi, size=$sz")
    dφ = dims_phi[1]

    if dφ == 3
        return vec(@view FDFcube[i, j, :])
    elseif dφ == 1
        return vec(@view FDFcube[:, i, j])
    else
        return vec(@view FDFcube[i, :, j])
    end
end

"""
    run_reversal_transition_job(...)

    Runs the merged reversal-map + transition analysis pipeline.
"""
function run_reversal_transition_job(cfg)::Dict{String,Any}
    enabled = task_enabled(cfg, ["tasks", "reversal_transition_job", "enabled"]; default=true)
    if !enabled
        return skipped_job_result("reversal_transition_job", "disabled by tasks.reversal_transition_job.enabled")
    end

    sim_root = resolve_simulations_root(cfg)
    sim = string(cfg_require(cfg, ["simulation", "name"]))
    los = require_los(string(cfg_require(cfg, ["simulation", "los"])))

    npix = Int(cfg_get(cfg, ["tasks", "reversals_map", "npix"]; default=256))
    lbox_pc = Float64(cfg_get(cfg, ["tasks", "reversals_map", "lbox_pc"]; default=50.0))
    b_scale = Float64(cfg_get(cfg, ["tasks", "reversals_map", "b_scale"]; default=1000.0))
    target_nrev = Int(cfg_get(cfg, ["tasks", "reversals_map", "target_nrev"]; default=24))

    smooth_win = Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "smooth_win"]; default=5))
    sign_eps_val = Float64(cfg_get(cfg, ["tasks", "reversal_transition_job", "sign_eps"]; default=0.0))
    deriv_tol = Float64(cfg_get(cfg, ["tasks", "reversal_transition_job", "deriv_tol"]; default=0.05))

    pmax_file = withfaraday_path(sim_root, sim, los, "Pmax.fits")

    bname, profile_fun = los_config(los)
    blos_file = simulation_field_path(sim_root, sim, bname)
    require_existing_files([pmax_file, blos_file]; context="reversal_transition_job")

    pmax = read_fits_f32(pmax_file)
    require_ndims(pmax, 2, "Pmax")
    size(pmax) == (npix, npix) || error("Pmax shape mismatch: $(size(pmax)) expected ($npix,$npix)")

    bcube = read_fits_f32(blos_file)
    require_ndims(bcube, 3, bname)
    size(bcube) == (npix, npix, npix) || error("$bname shape mismatch: $(size(bcube))")

    dist = collect(range(0.0, lbox_pc; length=npix))
    Δ = dist[2] - dist[1]

    # reversals map + target pixel
    nrev = Array{Int16}(undef, npix, npix)
    for j in 1:npix, i in 1:npix
        prof = b_scale .* profile_fun(i, j, bcube)
        prof_s = smooth_moving_average(prof, smooth_win)
        nrev[i, j] = Int16(length(reversal_indices(prof_s; eps=sign_eps_val)))
    end

    idx_exact = findall(nrev .== target_nrev)
    if isempty(idx_exact)
        best_i, best_j, best_d = 1, 1, typemax(Int)
        for j in 1:npix, i in 1:npix
            d = abs(Int(nrev[i, j]) - target_nrev)
            if d < best_d
                best_i, best_j, best_d = i, j, d
            end
        end
        i_t, j_t = best_i, best_j
    else
        i_t, j_t = Tuple(idx_exact[1])
    end

    # deciles from Pmax
    q10 = quantile(vec(pmax), 0.10)
    q90 = quantile(vec(pmax), 0.90)
    idxs_1 = findall(pmax .<= q10)
    idxs_10 = findall(pmax .>= q90)
    rng = MersenneTwister(Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "rng_seed"]; default=0)))
    i1, j1 = Tuple(rand(rng, idxs_1))
    i10, j10 = Tuple(rand(rng, idxs_10))

    prof_1 = b_scale .* profile_fun(i1, j1, bcube)
    prof_1_s = smooth_moving_average(prof_1, smooth_win)

    ks_1 = reversal_indices(Vector{Float64}(prof_1_s); eps=sign_eps_val)
    windows = Tuple{Int,Int,Int}[]
    for k0 in ks_1
        kL, kR = tight_bounds_for_reversal(Vector{Float64}(prof_1_s), Δ, k0, ks_1;
            deriv_tol=deriv_tol, max_half_width_pix=18)
        push!(windows, (k0, kL, kR))
    end
    win_merged = merge_intervals_local([(kL, kR) for (_, kL, kR) in windows])

    out_plot = standard_output_path(cfg, "reversal_transition_job", "overview", "pdf"; simu=sim, los=los)
    out_csv = standard_output_path(cfg, "reversal_transition_job", "windows", "csv"; simu=sim, los=los)
    out_log = standard_output_path(cfg, "reversal_transition_job", "summary", "log"; simu=sim, los=los)

    open(out_csv, "w") do io
        println(io, "k0,kmin,kmax,smin_pc,smax_pc,width_pc")
        for (k0, kL, kR) in windows
            println(io, @sprintf("%d,%d,%d,%.6f,%.6f,%.6f", k0, kL, kR, dist[kL], dist[kR], dist[kR]-dist[kL]))
        end
    end

    open(out_log, "w") do io
        println(io, "sim=$sim")
        println(io, "los=$los")
        println(io, "target_nrev=$target_nrev")
        println(io, "target_pixel=($i_t,$j_t)")
        println(io, "target_pixel_nrev=$(nrev[i_t,j_t])")
        println(io, "decile1_pixel=($i1,$j1)")
        println(io, "decile10_pixel=($i10,$j10)")
        println(io, "windows_merged=$(win_merged)")
    end

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1800, 780))

        ax1 = Axis(fig[1, 1],
            title="Nrev heatmap",
            xlabel="i", ylabel="j",
            xgridvisible=false, ygridvisible=false)
        hm1 = heatmap!(ax1, nrev; colormap=:viridis)
        scatter!(ax1, [i_t], [j_t]; color=:white, markersize=14, label="target nrev")
        scatter!(ax1, [i1], [j1]; color=:cyan, markersize=12, label="LOS profile")
        axislegend(ax1; position=:rt)
        Colorbar(fig[1, 2], hm1)

        ax2 = Axis(fig[1, 3],
            title=@sprintf("B_los LOS profile (i=%d, j=%d)", i1, j1),
            xlabel="Distance [pc]",
            ylabel="B_los [uG]",
            xgridvisible=false, ygridvisible=false)

        lines!(ax2, dist, prof_1_s; color=:black, linewidth=2)
        ytop = maximum(prof_1_s)
        yspan = max(maximum(prof_1_s) - minimum(prof_1_s), 1e-6)

        for (kL, kR) in win_merged
            xL = dist[kL]
            xR = dist[kR]
            width_pc = xR - xL
            vspan!(ax2, xL, xR; color=(:lightgray, 0.35))
            vlines!(ax2, [xL, xR]; color=:blue, linestyle=:dot, linewidth=2)
            text!(ax2, (xL + xR) / 2, ytop + 0.06 * yspan;
                  text=@sprintf("%.2f pc", width_pc),
                  align=(:center, :bottom),
                  fontsize=20,
                  color=:blue)
        end

        ylims!(ax2, minimum(prof_1_s) - 0.06 * yspan, ytop + 0.16 * yspan)

        save(out_plot, fig)
    end

    return Dict(
        "task" => "reversal_transition_job",
        "target_pixel" => Dict("i" => i_t, "j" => j_t, "nrev" => Int(nrev[i_t, j_t])),
        "decile_pixels" => Dict("d1" => [i1, j1], "d10" => [i10, j10]),
        "outputs" => Dict("plot" => out_plot, "csv" => out_csv, "log" => out_log),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_job_entrypoint("reversal_transition_job", run_reversal_transition_job)
end
