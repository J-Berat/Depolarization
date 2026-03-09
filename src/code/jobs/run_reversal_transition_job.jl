using CairoMakie
using LaTeXStrings
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
    resolve_phi_array(...)

    Resolves the `phi` axis from config (user-provided if present), otherwise default range.
"""
function resolve_phi_array(cfg, nphi::Int)
    phi_keys = [
        ["tasks", "reversal_transition_job", "phiArray"],
        ["tasks", "reversal_transition_job", "PhiArray"],
        ["tasks", "b_transition", "phiArray"],
        ["tasks", "b_transition", "PhiArray"],
        ["tasks", "instrumental_effect", "phiArray"],
        ["tasks", "instrumental_effect", "PhiArray"],
    ]

    for keys in phi_keys
        raw = cfg_get(cfg, keys; default=nothing)
        if raw !== nothing
            raw isa AbstractVector || error("$(join(keys, '.')) must be an array of numbers")
            phi = Float64.(collect(raw))
            isempty(phi) && error("$(join(keys, '.')) cannot be empty")
            return phi, true
        end
    end

    phi_min = cfg_get(cfg, ["tasks", "reversal_transition_job", "phi_min"]; default=nothing)
    phi_max = cfg_get(cfg, ["tasks", "reversal_transition_job", "phi_max"]; default=nothing)
    if phi_min !== nothing || phi_max !== nothing
        lo = Float64(phi_min === nothing ? -10.0 : phi_min)
        hi = Float64(phi_max === nothing ? 10.0 : phi_max)
        return collect(range(lo, hi; length=nphi)), true
    end

    return collect(range(-10.0, 10.0; length=nphi)), false
end

"""
    infer_nphi_from_fdf(...)

    Infers `nphi` from FDF cube shape when not explicitly provided.
"""
function infer_nphi_from_fdf(sz::NTuple{3,Int}, npix::Int, nphi_hint::Int)
    if nphi_hint > 0 && any(==(nphi_hint), sz)
        return nphi_hint
    end

    dims_not_npix = [d for d in sz if d != npix]
    if length(dims_not_npix) == 1
        return dims_not_npix[1]
    end

    for d in sort(collect(Set(sz)))
        if count(==(d), sz) == 1
            return d
        end
    end

    return minimum(sz)
end

"""
    compute_reversal_windows(...)

    Computes reversal indices, local windows, and merged windows for one LOS profile.
"""
function compute_reversal_windows(prof_s::Vector{Float64}, Δ::Real, sign_eps_val::Float64, deriv_tol::Float64)
    ks = reversal_indices(prof_s; eps=sign_eps_val)
    windows = Tuple{Int,Int,Int}[]
    for k0 in ks
        kL, kR = tight_bounds_for_reversal(prof_s, Δ, k0, ks;
            deriv_tol=deriv_tol, max_half_width_pix=18)
        push!(windows, (k0, kL, kR))
    end
    win_merged = merge_intervals_local([(kL, kR) for (_, kL, kR) in windows])
    return ks, windows, win_merged
end

"""
    plot_profile_with_windows!(...)

    Draws B profile with reversal windows and their distances.
"""
function plot_profile_with_windows!(ax, dist::Vector{Float64}, prof::Vector{Float64},
                                    win_merged::Vector{Tuple{Int,Int}}; color=:dodgerblue,
                                    annotation_size::Int=13)
    lines!(ax, dist, prof; color=color, linewidth=2.5)
    hlines!(ax, [0.0]; color=:blue, linestyle=:dash, linewidth=1.8)

    ylo = minimum(prof)
    yhi = maximum(prof)
    yspan = max(yhi - ylo, 1e-6)
    y_annot_base = yhi + 0.06 * yspan

    for (k, (kL, kR)) in enumerate(win_merged)
        xL = dist[kL]
        xR = dist[kR]
        width_pc = xR - xL
        y_annot = y_annot_base + 0.08 * yspan * ((k - 1) % 2)

        vspan!(ax, xL, xR; color=(:lightgray, 0.25))
        vlines!(ax, [xL, xR]; color=:gray30, linestyle=:dot, linewidth=1.4)
        lines!(ax, [xL, xR], [y_annot, y_annot]; color=:black, linewidth=1.3)
        text!(ax, (xL + xR) / 2, y_annot + 0.01 * yspan;
              text=latexstring(@sprintf("%.2f\\ \\mathrm{pc}", width_pc)),
              align=(:center, :bottom),
              fontsize=annotation_size,
              color=:black)
    end

    xlims!(ax, minimum(dist), maximum(dist))
    ylims!(ax, ylo - 0.08 * yspan, yhi + 0.26 * yspan)
    return nothing
end

"""
    run_reversal_transition_job(...)

    Runs the merged reversal-map + transition analysis pipeline.
"""
function run_reversal_transition_job(cfg)::Dict{String,Any}
    enabled = task_enabled(cfg, ["tasks", "reversal_transition_job", "enabled"]; default=true)
    if !enabled
        return skipped_job_result("reversal", "disabled by tasks.reversal_transition_job.enabled")
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
    nphi_hint = Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "nphi"];
                            default=cfg_get(cfg, ["tasks", "b_transition", "nphi"]; default=81)))

    pmax_file = withfaraday_path(sim_root, sim, los, "Pmax.fits")
    fdf_re_file = withfaraday_path(sim_root, sim, los, "realFDF.fits")
    fdf_im_file = withfaraday_path(sim_root, sim, los, "imagFDF.fits")

    bname, profile_fun = los_config(los)
    blos_file = simulation_field_path(sim_root, sim, bname)
    require_existing_files([pmax_file, blos_file, fdf_re_file, fdf_im_file]; context="reversal_transition_job")

    pmax = read_fits_f32(pmax_file)
    require_ndims(pmax, 2, "Pmax")
    size(pmax) == (npix, npix) || error("Pmax shape mismatch: $(size(pmax)) expected ($npix,$npix)")

    bcube = read_fits_f32(blos_file)
    require_ndims(bcube, 3, bname)
    size(bcube) == (npix, npix, npix) || error("$bname shape mismatch: $(size(bcube))")

    fdf_re = read_fits_f32(fdf_re_file)
    fdf_im = read_fits_f32(fdf_im_file)
    require_ndims(fdf_re, 3, "realFDF")
    require_ndims(fdf_im, 3, "imagFDF")
    size(fdf_re) == size(fdf_im) || error("realFDF/imagFDF size mismatch: $(size(fdf_re)) vs $(size(fdf_im))")

    nphi_auto = infer_nphi_from_fdf(size(fdf_re), npix, nphi_hint)
    phiArray, phi_from_cfg = resolve_phi_array(cfg, nphi_auto)
    if !any(==(length(phiArray)), size(fdf_re))
        if phi_from_cfg
            error("phiArray length=$(length(phiArray)) does not match any FDF axis size=$(size(fdf_re))")
        end
        nphi_auto = infer_nphi_from_fdf(size(fdf_re), npix, length(phiArray))
        phiArray = collect(range(-10.0, 10.0; length=nphi_auto))
    end
    nphi = length(phiArray)

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
    prof_10 = b_scale .* profile_fun(i10, j10, bcube)
    prof_1_s = smooth_moving_average(prof_1, smooth_win)
    prof_10_s = smooth_moving_average(prof_10, smooth_win)

    _, windows_1, win_merged_1 = compute_reversal_windows(Vector{Float64}(prof_1_s), Δ, sign_eps_val, deriv_tol)
    _, windows_10, win_merged_10 = compute_reversal_windows(Vector{Float64}(prof_10_s), Δ, sign_eps_val, deriv_tol)

    fdf_1 = hypot.(extract_fdf_spectrum(fdf_re, i1, j1, nphi), extract_fdf_spectrum(fdf_im, i1, j1, nphi))
    fdf_10 = hypot.(extract_fdf_spectrum(fdf_re, i10, j10, nphi), extract_fdf_spectrum(fdf_im, i10, j10, nphi))
    ipeak_1 = argmax(fdf_1)
    ipeak_10 = argmax(fdf_10)
    phi_peak_1 = phiArray[ipeak_1]
    phi_peak_10 = phiArray[ipeak_10]

    heatmap_label_size = Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "heatmap_label_size"];
                                     default=36))
    heatmap_tick_size = Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "heatmap_tick_size"];
                                    default=30))
    plot_label_size = Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "plot_label_size"];
                                  default=34))
    plot_tick_size = Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "plot_tick_size"];
                                 default=28))
    title_size = Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "title_size"]; default=38))
    annotation_size = Int(cfg_get(cfg, ["tasks", "reversal_transition_job", "annotation_size"]; default=18))

    out_plot = standard_output_path(cfg, "reversal", "overview", "pdf"; simu=sim, los=los)
    out_csv = standard_output_path(cfg, "reversal", "windows", "csv"; simu=sim, los=los)
    out_log = standard_output_path(cfg, "reversal", "summary", "log"; simu=sim, los=los)

    open(out_csv, "w") do io
        println(io, "decile,k0,kmin,kmax,smin_pc,smax_pc,width_pc")
        for (k0, kL, kR) in windows_1
            println(io, @sprintf("d1,%d,%d,%d,%.6f,%.6f,%.6f", k0, kL, kR, dist[kL], dist[kR], dist[kR]-dist[kL]))
        end
        for (k0, kL, kR) in windows_10
            println(io, @sprintf("d10,%d,%d,%d,%.6f,%.6f,%.6f", k0, kL, kR, dist[kL], dist[kR], dist[kR]-dist[kL]))
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
        println(io, "phi_length=$nphi")
        println(io, @sprintf("phi_peak_d1=%.6f", phi_peak_1))
        println(io, @sprintf("phi_peak_d10=%.6f", phi_peak_10))
        println(io, "windows_merged_d1=$(win_merged_1)")
        println(io, "windows_merged_d10=$(win_merged_10)")
    end
    @info "Reversal transition text outputs written" csv=out_csv log=out_log

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(2000, 1300))
        colgap!(fig.layout, 12)
        rowgap!(fig.layout, 18)

        pmin, pmaxv = extrema(finite_values(pmax))
        pmax_cr = (pmin, pmaxv)

        ax_h1 = Axis(fig[1, 1],
            title=latexstring(@sprintf("1^{\\mathrm{st}}\\ \\mathrm{decile}\\ (i=%d,\\ j=%d)", i1, j1)),
            xlabel=L"i", ylabel=L"j",
            titlesize=title_size,
            xlabelsize=heatmap_label_size, ylabelsize=heatmap_label_size,
            xticklabelsize=heatmap_tick_size, yticklabelsize=heatmap_tick_size,
            xgridvisible=false, ygridvisible=false)
        hm = heatmap!(ax_h1, pmax; colormap=:magma, colorrange=pmax_cr)
        vlines!(ax_h1, [i1]; color=:white, linestyle=:dash, linewidth=2.0)
        hlines!(ax_h1, [j1]; color=:white, linestyle=:dash, linewidth=2.0)
        xlims!(ax_h1, 1, npix)
        ylims!(ax_h1, 1, npix)

        ax_h10 = Axis(fig[2, 1],
            title=latexstring(@sprintf("10^{\\mathrm{th}}\\ \\mathrm{decile}\\ (i=%d,\\ j=%d)", i10, j10)),
            xlabel=L"i", ylabel=L"j",
            titlesize=title_size,
            xlabelsize=heatmap_label_size, ylabelsize=heatmap_label_size,
            xticklabelsize=heatmap_tick_size, yticklabelsize=heatmap_tick_size,
            xgridvisible=false, ygridvisible=false)
        heatmap!(ax_h10, pmax; colormap=:magma, colorrange=pmax_cr)
        vlines!(ax_h10, [i10]; color=:white, linestyle=:dash, linewidth=2.0)
        hlines!(ax_h10, [j10]; color=:white, linestyle=:dash, linewidth=2.0)
        xlims!(ax_h10, 1, npix)
        ylims!(ax_h10, 1, npix)

        Colorbar(fig[1:2, 2], hm; label=L"P_{\max}\ [\mathrm{K}]",
            labelsize=heatmap_label_size, ticklabelsize=heatmap_tick_size)

        g_right_top = GridLayout(fig[1, 3])
        rowgap!(g_right_top, 0)

        ax_b1 = Axis(g_right_top[1, 1],
            title="",
            xlabel=L"\mathrm{Distance}\ [\mathrm{pc}]",
            ylabel=L"B\ [\mu\mathrm{G}]",
            xaxisposition=:top,
            titlesize=title_size,
            xlabelsize=plot_label_size, ylabelsize=plot_label_size,
            xticklabelsize=plot_tick_size, yticklabelsize=plot_tick_size,
            xgridvisible=false, ygridvisible=false)
        plot_profile_with_windows!(ax_b1, dist, Vector{Float64}(prof_1_s), win_merged_1;
            color=:dodgerblue, annotation_size=annotation_size)

        ax_f1 = Axis(g_right_top[2, 1],
            title="",
            xlabel=L"\phi\ [\mathrm{rad}\ \mathrm{m}^{-2}]",
            ylabel=L"\left|F(\phi)\right|",
            xlabelsize=plot_label_size, ylabelsize=plot_label_size,
            xticklabelsize=plot_tick_size, yticklabelsize=plot_tick_size,
            xgridvisible=false, ygridvisible=false)
        lines!(ax_f1, phiArray, fdf_1; color=:dodgerblue, linewidth=2)
        xlims!(ax_f1, minimum(phiArray), maximum(phiArray))
        rowsize!(g_right_top, 1, Relative(0.6))
        rowsize!(g_right_top, 2, Relative(0.4))

        g_right_bottom = GridLayout(fig[2, 3])
        rowgap!(g_right_bottom, 0)

        ax_b10 = Axis(g_right_bottom[1, 1],
            title="",
            xlabel=L"\mathrm{Distance}\ [\mathrm{pc}]",
            ylabel=L"B\ [\mu\mathrm{G}]",
            xaxisposition=:top,
            titlesize=title_size,
            xlabelsize=plot_label_size, ylabelsize=plot_label_size,
            xticklabelsize=plot_tick_size, yticklabelsize=plot_tick_size,
            xgridvisible=false, ygridvisible=false)
        plot_profile_with_windows!(ax_b10, dist, Vector{Float64}(prof_10_s), win_merged_10;
            color=:red, annotation_size=annotation_size)

        ax_f10 = Axis(g_right_bottom[2, 1],
            title="",
            xlabel=L"\phi\ [\mathrm{rad}\ \mathrm{m}^{-2}]",
            ylabel=L"\left|F(\phi)\right|",
            xlabelsize=plot_label_size, ylabelsize=plot_label_size,
            xticklabelsize=plot_tick_size, yticklabelsize=plot_tick_size,
            xgridvisible=false, ygridvisible=false)
        lines!(ax_f10, phiArray, fdf_10; color=:dodgerblue, linewidth=2)
        xlims!(ax_f10, minimum(phiArray), maximum(phiArray))
        rowsize!(g_right_bottom, 1, Relative(0.6))
        rowsize!(g_right_bottom, 2, Relative(0.4))

        save(out_plot, fig)
    end
    @info "Reversal transition plot written" plot=out_plot

    return Dict(
        "task" => "reversal",
        "target_pixel" => Dict("i" => i_t, "j" => j_t, "nrev" => Int(nrev[i_t, j_t])),
        "decile_pixels" => Dict("d1" => [i1, j1], "d10" => [i10, j10]),
        "phi" => Dict("nphi" => nphi, "peak_d1" => phi_peak_1, "peak_d10" => phi_peak_10),
        "outputs" => Dict("plot" => out_plot, "csv" => out_csv, "log" => out_log),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_job_entrypoint("reversal", run_reversal_transition_job)
end
