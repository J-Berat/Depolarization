using CairoMakie
using LaTeXStrings
using Statistics
using StatsBase

include(joinpath(@__DIR__, "..", "lib", "DepolLib.jl"))
using .DepolLib

"""
    _sky_plane_labels(...)

Returns axis labels for the sky plane associated with LOS.
"""
function _sky_plane_labels(los::AbstractString)
    los = require_los(los)
    if los == "x"
        return "y", "z"
    elseif los == "y"
        return "x", "z"
    else
        return "x", "y"
    end
end

"""
    _bparallel_from_components(...)

Selects LOS magnetic component from `(Bx, By, Bz)`.
"""
function _bparallel_from_components(Bx::AbstractArray, By::AbstractArray, Bz::AbstractArray, los::AbstractString)
    los = require_los(los)
    if los == "x"
        return Bx
    elseif los == "y"
        return By
    else
        return Bz
    end
end

"""
    _cb_map(...)

Computes `Cb = |sum(Bpar)| / sum(|Bpar|)` along LOS.
"""
function _cb_map(Bpar::AbstractArray{<:Real,3}, los::AbstractString; eps::Real=1f-30)
    ax = los_axis(los)
    num = abs.(sum(Bpar; dims=ax))
    den = sum(abs.(Bpar); dims=ax) .+ eps
    return dropdims(num ./ den; dims=ax)
end

"""
    _cphi_map(...)

Computes `Cphi = |sum(ne*Bpar)| / sum(|ne*Bpar|)` along LOS.
"""
function _cphi_map(ne::AbstractArray{<:Real,3}, Bpar::AbstractArray{<:Real,3}, los::AbstractString; eps::Real=1f-30)
    size(ne) == size(Bpar) || error("ne and Bpar size mismatch: ne=$(size(ne)) Bpar=$(size(Bpar))")
    ax = los_axis(los)
    X = ne .* Bpar
    num = abs.(sum(X; dims=ax))
    den = sum(abs.(X); dims=ax) .+ eps
    return dropdims(num ./ den; dims=ax)
end

"""
    _ab_map(...)

Computes `<|Bpar|>` along LOS.
"""
function _ab_map(Bpar::AbstractArray{<:Real,3}, los::AbstractString)
    ax = los_axis(los)
    return dropdims(mean(abs.(Bpar); dims=ax); dims=ax)
end

"""
    _nrev_map(...)

Counts sign reversals along LOS for each sky-plane pixel.
"""
function _nrev_map(Bpar::AbstractArray{<:Real,3}, los::AbstractString; thresh::Real=0.0)
    los = require_los(los)
    nx, ny, nz = size(Bpar)
    if los == "x"
        out = zeros(Int16, ny, nz)
        @inbounds for j in 1:nz, i in 1:ny
            prev = 0
            n = 0
            for k in 1:nx
                v = Bpar[k, i, j]
                s = abs(v) <= thresh ? 0 : (v > 0 ? 1 : -1)
                if prev != 0 && s != 0 && s != prev
                    n += 1
                end
                if s != 0
                    prev = s
                end
            end
            out[i, j] = Int16(n)
        end
        return out
    elseif los == "y"
        out = zeros(Int16, nx, nz)
        @inbounds for j in 1:nz, i in 1:nx
            prev = 0
            n = 0
            for k in 1:ny
                v = Bpar[i, k, j]
                s = abs(v) <= thresh ? 0 : (v > 0 ? 1 : -1)
                if prev != 0 && s != 0 && s != prev
                    n += 1
                end
                if s != 0
                    prev = s
                end
            end
            out[i, j] = Int16(n)
        end
        return out
    else
        out = zeros(Int16, nx, ny)
        @inbounds for j in 1:ny, i in 1:nx
            prev = 0
            n = 0
            for k in 1:nz
                v = Bpar[i, j, k]
                s = abs(v) <= thresh ? 0 : (v > 0 ? 1 : -1)
                if prev != 0 && s != 0 && s != prev
                    n += 1
                end
                if s != 0
                    prev = s
                end
            end
            out[i, j] = Int16(n)
        end
        return out
    end
end

"""
    _canal_mask(...)

Returns low-Pmax mask and threshold from quantile.
"""
function _canal_mask(Pmax::AbstractArray{<:Real,2}; decile::Real=0.1)
    0.0 < decile < 1.0 || error("decile must be in (0,1), got $decile")
    vals = vec(Float64.(Pmax))
    vals = vals[isfinite.(vals)]
    isempty(vals) && error("Pmax has no finite values")
    pthr = quantile(vals, decile)
    return (Pmax .<= pthr), pthr
end

"""
    _pdf_hist!(...)

Plots normalized PDFs for canal and non-canal populations with shared edges.
"""
function _pdf_hist!(ax, data_canal::AbstractVector, data_out::AbstractVector; nbins::Int=60, xlabel::AbstractString="")
    d1 = data_canal[isfinite.(data_canal)]
    d2 = data_out[isfinite.(data_out)]
    (isempty(d1) || isempty(d2)) && return

    lo = min(minimum(d1), minimum(d2))
    hi = max(maximum(d1), maximum(d2))
    lo < hi || return

    edges = range(lo, hi; length=nbins + 1)
    h1 = fit(Histogram, d1, edges)
    h2 = fit(Histogram, d2, edges)
    s1 = sum(h1.weights)
    s2 = sum(h2.weights)
    (s1 > 0 && s2 > 0) || return

    xs = (edges[1:end-1] .+ edges[2:end]) ./ 2
    dx = edges[2] - edges[1]
    y1 = (h1.weights ./ s1) ./ dx
    y2 = (h2.weights ./ s2) ./ dx

    lines!(ax, xs, y1, label=latexstring("\\mathrm{canal}"))
    lines!(ax, xs, y2, label=latexstring("\\mathrm{non\\text{-}canal}"))
    ax.xlabel = xlabel
    ax.ylabel = "PDF"
    return nothing
end

"""
    _binned_stats(...)

Returns median and [16,84]% intervals of `y` in `x` bins.
"""
function _binned_stats(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; nbins::Int=30, mincount::Int=10)
    good = isfinite.(x) .& isfinite.(y)
    xv = Float64.(x[good])
    yv = Float64.(y[good])
    isempty(xv) && return Float64[], Float64[], Float64[], Float64[]

    xmin = minimum(xv)
    xmax = maximum(xv)
    xmin < xmax || return Float64[], Float64[], Float64[], Float64[]

    edges = range(xmin, xmax; length=nbins + 1)
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])

    med = fill(NaN, nbins)
    q16 = fill(NaN, nbins)
    q84 = fill(NaN, nbins)
    for b in 1:nbins
        lo, hi = edges[b], edges[b + 1]
        idx = (xv .>= lo) .& (xv .< hi)
        yy = yv[idx]
        if length(yy) >= mincount
            med[b] = quantile(yy, 0.50)
            q16[b] = quantile(yy, 0.16)
            q84[b] = quantile(yy, 0.84)
        end
    end

    return centers, med, q16, q84
end

"""
    _bin2d_median(...)

Returns median `z` in 2D `(x,y)` bins.
"""
function _bin2d_median(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, z::AbstractVector{<:Real};
                       nx::Int=35, ny::Int=35, mincount::Int=10)
    good = isfinite.(x) .& isfinite.(y) .& isfinite.(z)
    xv = Float64.(x[good])
    yv = Float64.(y[good])
    zv = Float64.(z[good])
    isempty(xv) && return Float64[], Float64[], fill(NaN, nx, ny)

    xmin = minimum(xv)
    xmax = maximum(xv)
    ymin = minimum(yv)
    ymax = maximum(yv)
    (xmin < xmax && ymin < ymax) || return Float64[], Float64[], fill(NaN, nx, ny)

    xedges = range(xmin, xmax; length=nx + 1)
    yedges = range(ymin, ymax; length=ny + 1)
    xc = 0.5 .* (xedges[1:end-1] .+ xedges[2:end])
    yc = 0.5 .* (yedges[1:end-1] .+ yedges[2:end])

    zmed = fill(NaN, nx, ny)
    buckets = [Float64[] for _ in 1:(nx * ny)]
    @inbounds for t in eachindex(xv)
        ix = searchsortedlast(xedges, xv[t]) - 1
        iy = searchsortedlast(yedges, yv[t]) - 1
        if 1 <= ix <= nx && 1 <= iy <= ny
            push!(buckets[(iy - 1) * nx + ix], zv[t])
        end
    end

    @inbounds for iy in 1:ny, ix in 1:nx
        b = buckets[(iy - 1) * nx + ix]
        if length(b) >= mincount
            zmed[ix, iy] = median(b)
        end
    end

    return xc, yc, zmed
end

"""
    run_canal_metrics_job(...)

Computes LOS coherence metrics (`Cb`, `Cphi`, `<|Bpar|>`, `Nrev`) and writes
map/PDF/trend figures with canal-vs-non-canal comparisons from `Pmax`.
"""
function run_canal_metrics_job(cfg)::Dict{String,Any}
    enabled = task_enabled(cfg, ["tasks", "canal_metrics", "enabled"]; default=true)
    if !enabled
        return skipped_job_result("canal_metrics", "disabled by tasks.canal_metrics.enabled")
    end

    sim_root = resolve_simulations_root(cfg)
    sim = string(cfg_require(cfg, ["simulation", "name"]))
    los = require_los(string(cfg_require(cfg, ["simulation", "los"])))

    decile = Float64(cfg_get(cfg, ["tasks", "canal_metrics", "decile"]; default=0.1))
    b_unit_factor = Float64(cfg_get(cfg, ["tasks", "canal_metrics", "b_unit_factor"];
                                    default=cfg_get(cfg, ["tasks", "reversals_map", "b_scale"]; default=1000.0)))
    lbox_pc = Float64(cfg_get(cfg, ["tasks", "canal_metrics", "lbox_pc"];
                              default=cfg_get(cfg, ["tasks", "reversals_map", "lbox_pc"]; default=50.0)))
    nrev_thresh = Float64(cfg_get(cfg, ["tasks", "canal_metrics", "nrev_thresh"]; default=0.0))
    nbins_pdf = Int(cfg_get(cfg, ["tasks", "canal_metrics", "nbins_pdf"]; default=60))
    nbins_nrev = Int(cfg_get(cfg, ["tasks", "canal_metrics", "nbins_nrev"]; default=30))
    nbins_stats = Int(cfg_get(cfg, ["tasks", "canal_metrics", "nbins_stats"]; default=30))
    nbins_matrix_x = Int(cfg_get(cfg, ["tasks", "canal_metrics", "nbins_matrix_x"]; default=35))
    nbins_matrix_y = Int(cfg_get(cfg, ["tasks", "canal_metrics", "nbins_matrix_y"]; default=35))
    mincount = Int(cfg_get(cfg, ["tasks", "canal_metrics", "mincount"]; default=10))

    pmax_file = withfaraday_path(sim_root, sim, los, "Pmax.fits")
    ne_file = synchrotron_path(sim_root, sim, los, "ne.fits")
    bx_file = simulation_field_path(sim_root, sim, "Bx")
    by_file = simulation_field_path(sim_root, sim, "By")
    bz_file = simulation_field_path(sim_root, sim, "Bz")
    require_existing_files([pmax_file, ne_file, bx_file, by_file, bz_file]; context="canal_metrics_job")

    Pmax = read_fits_f32(pmax_file)
    ne = read_fits_f32(ne_file)
    Bx = read_fits_f32(bx_file)
    By = read_fits_f32(by_file)
    Bz = read_fits_f32(bz_file)

    require_ndims(Pmax, 2, "Pmax")
    require_ndims(ne, 3, "ne")
    require_ndims(Bx, 3, "Bx")
    require_ndims(By, 3, "By")
    require_ndims(Bz, 3, "Bz")
    require_same_size([ne, Bx, By, Bz], ["ne", "Bx", "By", "Bz"])

    Bpar = _bparallel_from_components(Bx, By, Bz, los) .* Float32(b_unit_factor)
    Cb = _cb_map(Bpar, los)
    Cphi = _cphi_map(ne, Bpar, los)
    Ab = _ab_map(Bpar, los)
    Nrev = _nrev_map(Bpar, los; thresh=nrev_thresh)
    size(Pmax) == size(Cb) || error("Pmax/Cb shape mismatch: Pmax=$(size(Pmax)) Cb=$(size(Cb))")

    canal, pthr = _canal_mask(Pmax; decile=decile)

    vP = vec(Float64.(Pmax))
    vCb = vec(Float64.(Cb))
    vCphi = vec(Float64.(Cphi))
    vAb = vec(Float64.(Ab))
    vNr = vec(Float64.(Nrev))
    vcan = vec(canal)

    xlab, ylab = _sky_plane_labels(los)
    xticks = ticks_pc(size(Cb, 1); lbox_pc=lbox_pc, step_pc=10.0)
    yticks = ticks_pc(size(Cb, 2); lbox_pc=lbox_pc, step_pc=10.0)

    out_maps = standard_output_path(cfg, "canal_metrics", "maps_cb_cphi", "pdf"; simu=sim, los=los)
    out_pdfs = standard_output_path(cfg, "canal_metrics", "pdfs_2x2", "pdf"; simu=sim, los=los)
    out_matrix = standard_output_path(cfg, "canal_metrics", "cb_ab_matrix", "pdf"; simu=sim, los=los)

    with_theme(theme_latexfonts()) do
        fs_title = 26
        fs_label = 26
        fs_tick = 22

        figA = Figure(size=(1250, 560))

        ax1 = Axis(figA[1, 1],
            title=L"C_B=\frac{\left|\sum_k B_{\parallel}(k)\right|}{\sum_k \left|B_{\parallel}(k)\right|}",
            xlabel=string(xlab, " [pc]"), ylabel=string(ylab, " [pc]"),
            aspect=DataAspect(),
            xgridvisible=false, ygridvisible=false,
            titlesize=fs_title, xlabelsize=fs_label, ylabelsize=fs_label,
            xticklabelsize=fs_tick, yticklabelsize=fs_tick)
        ax1.xticks = xticks
        ax1.yticks = yticks
        hm1 = heatmap!(ax1, Cb; colorrange=(0, 1), colormap=:magma, interpolate=false)
        contour!(ax1, Float32.(canal); levels=[0.5], color=(:white, 0.8), linewidth=2, linestyle=:dash)
        Colorbar(figA[1, 2], hm1, label="")

        ax2 = Axis(figA[1, 3],
            title=L"C_{\phi}=\frac{\left|\sum_k n_e(k)\,B_{\parallel}(k)\right|}{\sum_k \left|n_e(k)\,B_{\parallel}(k)\right|}",
            xlabel=string(xlab, " [pc]"), ylabel=string(ylab, " [pc]"),
            aspect=DataAspect(),
            xgridvisible=false, ygridvisible=false,
            titlesize=fs_title, xlabelsize=fs_label, ylabelsize=fs_label,
            xticklabelsize=fs_tick, yticklabelsize=fs_tick)
        ax2.xticks = xticks
        ax2.yticks = yticks
        hm2 = heatmap!(ax2, Cphi; colorrange=(0, 1), colormap=:magma, interpolate=false)
        contour!(ax2, Float32.(canal); levels=[0.5], color=(:white, 0.8), linewidth=2, linestyle=:dash)
        Colorbar(figA[1, 4], hm2, label="")

        save(out_maps, figA)

        figB = Figure(size=(1150, 820))

        axB1 = Axis(figB[1, 1], title="Cb",
            xgridvisible=false, ygridvisible=false,
            titlesize=fs_title, xlabelsize=fs_label, ylabelsize=fs_label,
            xticklabelsize=fs_tick, yticklabelsize=fs_tick)
        _pdf_hist!(axB1, vCb[vcan], vCb[.!vcan]; nbins=nbins_pdf, xlabel="Cb")
        axislegend(axB1, position=:rt, framevisible=true)

        axB2 = Axis(figB[1, 2], title="Cphi",
            xgridvisible=false, ygridvisible=false,
            titlesize=fs_title, xlabelsize=fs_label, ylabelsize=fs_label,
            xticklabelsize=fs_tick, yticklabelsize=fs_tick)
        _pdf_hist!(axB2, vCphi[vcan], vCphi[.!vcan]; nbins=nbins_pdf, xlabel="Cphi")
        axislegend(axB2, position=:rt, framevisible=true)

        axB3 = Axis(figB[2, 1], title="<|Bpar|>",
            xgridvisible=false, ygridvisible=false,
            titlesize=fs_title, xlabelsize=fs_label, ylabelsize=fs_label,
            xticklabelsize=fs_tick, yticklabelsize=fs_tick)
        _pdf_hist!(axB3, vAb[vcan], vAb[.!vcan]; nbins=nbins_pdf, xlabel="<|Bpar|> [uG]")
        axislegend(axB3, position=:rt, framevisible=true)

        axB4 = Axis(figB[2, 2], title="Nrev",
            xgridvisible=false, ygridvisible=false,
            titlesize=fs_title, xlabelsize=fs_label, ylabelsize=fs_label,
            xticklabelsize=fs_tick, yticklabelsize=fs_tick)
        _pdf_hist!(axB4, vNr[vcan], vNr[.!vcan]; nbins=nbins_nrev, xlabel="Nrev")
        axislegend(axB4, position=:rt, framevisible=true)

        save(out_pdfs, figB)

        figC = Figure(size=(1000, 800))

        axC1 = Axis(figC[1, 1],
            xlabel="Cb", ylabel="<|Bpar|> [uG]",
            title="median(Pmax) in (Cb, <|Bpar|>) bins",
            xgridvisible=false, ygridvisible=false,
            titlesize=fs_title, xlabelsize=fs_label, ylabelsize=fs_label,
            xticklabelsize=fs_tick, yticklabelsize=fs_tick)
        xc, yc, Zmed = _bin2d_median(vCb, vAb, vP; nx=nbins_matrix_x, ny=nbins_matrix_y, mincount=mincount)
        if !isempty(xc) && !isempty(yc)
            hmC = heatmap!(axC1, xc, yc, Zmed; nan_color=:transparent, colormap=:magma)
            Colorbar(figC[1, 2], hmC, label="median(Pmax)")
        end

        axC2 = Axis(figC[2, 1],
            xlabel="Cb", ylabel="Pmax",
            title="Pmax vs Cb (median and [16,84]%)",
            xgridvisible=false, ygridvisible=false,
            titlesize=fs_title, xlabelsize=fs_label, ylabelsize=fs_label,
            xticklabelsize=fs_tick, yticklabelsize=fs_tick)
        centers, med, q16, q84 = _binned_stats(vCb, vP; nbins=nbins_stats, mincount=mincount)
        if !isempty(centers)
            good = isfinite.(med) .& isfinite.(q16) .& isfinite.(q84)
            lines!(axC2, centers[good], med[good])
            band!(axC2, centers[good], q16[good], q84[good])
        end

        save(out_matrix, figC)
    end

    return Dict(
        "task" => "canal_metrics",
        "simulation" => sim,
        "los" => los,
        "thresholds" => Dict(
            "pmax_decile" => decile,
            "pmax_threshold" => pthr,
            "nrev_thresh" => nrev_thresh,
        ),
        "outputs" => Dict(
            "maps_cb_cphi" => out_maps,
            "pdfs_2x2" => out_pdfs,
            "cb_ab_matrix" => out_matrix,
        ),
        "summary" => Dict(
            "npix_total" => length(vcan),
            "npix_canal" => count(vcan),
            "canal_fraction" => mean(Float64.(vcan)),
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_job_entrypoint("canal_metrics", run_canal_metrics_job)
end
