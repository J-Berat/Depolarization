using TOML
using CairoMakie
using LaTeXStrings
using Statistics
using Printf
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "code", "instrumental_effect", "InstrumentalEffect.jl"))

using .InstrumentalEffect

const CFG_PATH = joinpath(@__DIR__, "..", "config", "default.toml")
const MIN_AREA_PIX = 12

function _load_paths(cfg_path::AbstractString; los_override::Union{Nothing,String}=nothing)
    cfg = TOML.parsefile(cfg_path)
    sim_root = String(cfg["paths"]["simulations_root"])
    sim = String(cfg["simulation"]["name"])
    los = isnothing(los_override) ? String(cfg["simulation"]["los"]) : String(los_override)
    llarge_list = Float64.(cfg["tasks"]["instrumental_effect"]["llarge_list"])
    lbox_pc = Float64(get(cfg["tasks"]["canal_metrics"], "lbox_pc", 50.0))

    base_faraday = joinpath(sim_root, sim, los, "Synchrotron", "WithFaraday")
    base_out = joinpath(@__DIR__, "..", "outputs", "instrumental", sim, "LOS$(los)")
    return (
        sim_root=sim_root,
        sim=sim,
        los=los,
        llarge_list=llarge_list,
        lbox_pc=lbox_pc,
        q_in=joinpath(base_faraday, "Qnu.fits"),
        u_in=joinpath(base_faraday, "Unu.fits"),
        pmax0_in=joinpath(base_faraday, "Pmax.fits"),
        q_phi_in=joinpath(base_faraday, "realFDF.fits"),
        u_phi_in=joinpath(base_faraday, "imagFDF.fits"),
        bx_in=joinpath(sim_root, sim, "Bx.fits"),
        by_in=joinpath(sim_root, sim, "By.fits"),
        bz_in=joinpath(sim_root, sim, "Bz.fits"),
        dens_in=joinpath(sim_root, sim, "density.fits"),
        base_out=normpath(base_out),
    )
end

function _build_cfg(paths)
    return InstrumentalConfig(
        Q_in = paths.q_in,
        U_in = paths.u_in,
        Pmax_nofilter_path = paths.pmax0_in,
        Q_in_phi = paths.q_phi_in,
        U_in_phi = paths.u_phi_in,
        base_out = paths.base_out,
        Bx_in = paths.bx_in,
        By_in = paths.by_in,
        Bz_in = paths.bz_in,
        dens_in = paths.dens_in,
        los = paths.los,
    )
end

function _map_label(llarge_eff_pc::Real)
    return string(round(Int, llarge_eff_pc))
end

function _pc_ticks(lbox_pc::Real)
    ticks = collect(0.0:10.0:Float64(lbox_pc))
    last(ticks) == Float64(lbox_pc) || push!(ticks, Float64(lbox_pc))
    return (ticks, string.(Int.(round.(ticks))))
end

function _line_mean_orientation(θvals::AbstractVector{<:Real})
    isempty(θvals) && return NaN
    c2 = mean(cos.(2.0 .* θvals))
    s2 = mean(sin.(2.0 .* θvals))
    return mod(0.5 * atan(s2, c2), π)
end

function _connected_components(mask::AbstractMatrix{Bool})
    nx, ny = size(mask)
    labels = zeros(Int, nx, ny)
    comps = Vector{Vector{CartesianIndex{2}}}()
    neighbors = (
        CartesianIndex(-1, -1), CartesianIndex(-1, 0), CartesianIndex(-1, 1),
        CartesianIndex(0, -1),                         CartesianIndex(0, 1),
        CartesianIndex(1, -1),  CartesianIndex(1, 0), CartesianIndex(1, 1),
    )

    next_label = 0
    for j in 1:ny, i in 1:nx
        if !mask[i, j] || labels[i, j] != 0
            continue
        end
        next_label += 1
        queue = CartesianIndex{2}[CartesianIndex(i, j)]
        labels[i, j] = next_label
        comp = CartesianIndex{2}[]
        qhead = 1
        while qhead <= length(queue)
            idx = queue[qhead]
            qhead += 1
            push!(comp, idx)
            for δ in neighbors
                nb = idx + δ
                ii, jj = Tuple(nb)
                if 1 <= ii <= nx && 1 <= jj <= ny && mask[ii, jj] && labels[ii, jj] == 0
                    labels[ii, jj] = next_label
                    push!(queue, nb)
                end
            end
        end
        push!(comps, comp)
    end
    return comps
end

function _component_orientation(comp::Vector{CartesianIndex{2}}, x_pc::AbstractVector, y_pc::AbstractVector)
    xs = Float64[x_pc[I[1]] for I in comp]
    ys = Float64[y_pc[I[2]] for I in comp]
    xc = mean(xs)
    yc = mean(ys)
    dx = xs .- xc
    dy = ys .- yc
    Cxx = mean(dx .^ 2)
    Cxy = mean(dx .* dy)
    Cyy = mean(dy .^ 2)
    vals, vecs = eigen(Symmetric([Cxx Cxy; Cxy Cyy]))
    λmin, λmax = vals[1], vals[2]
    v = vecs[:, 2]
    θ = mod(atan(v[2], v[1]), π)
    aspect = sqrt(max(λmax, 1e-12) / max(λmin, 1e-12))
    return (; xc, yc, θ, λmax, λmin, aspect)
end

function _zone_summary(Pmap::AbstractMatrix, θB::AbstractMatrix, x_pc::AbstractVector, y_pc::AbstractVector;
        ridge_smooth_radius_pix::Int=InstrumentalEffect._CANAL_RIDGE_SMOOTH_RADIUS_PIX)
    Δx = x_pc[2] - x_pc[1]
    Δy = y_pc[2] - y_pc[1]
    score = InstrumentalEffect._canal_score_map(Float64.(Pmap))
    Ds = InstrumentalEffect._box_mean(score, ridge_smooth_radius_pix)
    thr = InstrumentalEffect._finite_quantile(Ds, InstrumentalEffect._CANAL_SCORE_QUANTILE)
    zone_mask = isfinite.(Ds) .& (Ds .>= thr)
    comps = _connected_components(zone_mask)
    raw_zone_count = length(comps)

    rows = NamedTuple[]
    overlay_segments = Tuple{Point2f, Point2f}[]
    filtered_mask = falses(size(zone_mask))

    for comp in comps
        length(comp) >= MIN_AREA_PIX || continue
        geom = _component_orientation(comp, x_pc, y_pc)

        θBvals = Float64[]
        for I in comp
            θ = θB[I]
            isfinite(θ) && push!(θBvals, θ)
        end
        isempty(θBvals) && continue
        θBzone = _line_mean_orientation(θBvals)
        delta_deg = abs(InstrumentalEffect._line_orientation_delta_deg(geom.θ, θBzone))
        area_pix = length(comp)
        area_pc2 = area_pix * Δx * Δy
        len_pc = clamp(3.0 * sqrt(max(geom.λmax, 1e-8)), 1.0, 5.0)
        ux = cos(geom.θ)
        uy = sin(geom.θ)
        p1 = Point2f(geom.xc - 0.5 * len_pc * ux, geom.yc - 0.5 * len_pc * uy)
        p2 = Point2f(geom.xc + 0.5 * len_pc * ux, geom.yc + 0.5 * len_pc * uy)
        push!(overlay_segments, (p1, p2))

        for I in comp
            filtered_mask[I] = true
        end

        push!(rows, (
            area_pix=area_pix,
            area_pc2=area_pc2,
            aspect=geom.aspect,
            lambda_min=geom.λmin,
            lambda_max=geom.λmax,
            mu_min=sqrt(max(geom.λmin, 0.0)),
            mu_max=sqrt(max(geom.λmax, 0.0)),
            mu_ratio=sqrt(max(geom.λmin, 1e-12) / max(geom.λmax, 1e-12)),
            theta_zone_deg=rad2deg(geom.θ),
            theta_B_deg=rad2deg(θBzone),
            abs_delta_deg=delta_deg,
            xc=geom.xc,
            yc=geom.yc,
        ))
    end

    deltas = [r.abs_delta_deg for r in rows]
    aspects = [r.aspect for r in rows]
    areas = [r.area_pix for r in rows]
    return (
        score=score,
        Ds=Ds,
        zone_mask=filtered_mask,
        rows=rows,
        raw_zone_count=raw_zone_count,
        zone_count=length(rows),
        area_fraction=mean(filtered_mask),
        median_area_pix=isempty(areas) ? NaN : median(areas),
        median_aspect=isempty(aspects) ? NaN : median(aspects),
        mean_abs_delta_deg=isempty(deltas) ? NaN : mean(deltas),
        median_abs_delta_deg=isempty(deltas) ? NaN : median(deltas),
        overlay_segments=overlay_segments,
    )
end

function _write_csv(path::AbstractString, labels, summaries)
    open(path, "w") do io
        println(io, "map_tag,raw_zone_count,zone_count,area_fraction,median_area_pix,median_aspect,mean_abs_delta_deg,median_abs_delta_deg")
        for tag in labels
            s = summaries[tag]
            @printf(io, "%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                tag, s.raw_zone_count, s.zone_count, s.area_fraction, s.median_area_pix,
                s.median_aspect, s.mean_abs_delta_deg, s.median_abs_delta_deg)
        end
    end
end

function _pmax_colorrange(maps)
    vals = Float64[]
    for P in values(maps)
        for v in P
            isfinite(v) && push!(vals, Float64(v))
        end
    end
    return (quantile(vals, 0.01), quantile(vals, 0.99))
end

function _plot_zone_test(pdf_path::AbstractString, png_path::AbstractString, labels, maps, summaries, lbox_pc::Real)
    firstmap = first(values(maps))
    x_edges = collect(range(0.0, Float64(lbox_pc), length=size(firstmap, 1) + 1))
    y_edges = collect(range(0.0, Float64(lbox_pc), length=size(firstmap, 2) + 1))
    pmax_colorrange = _pmax_colorrange(maps)
    ticks = _pc_ticks(lbox_pc)
    zone_cmap = cgrad([RGBAf(1,1,1,0.0), RGBAf(1,1,1,0.23)])

    fig = Figure(size=(3000, 1280), figure_padding=(16, 16, 12, 12))
    for (col, tag) in enumerate(labels)
        for row in 1:2
            ax = Axis(
                fig[row, col],
                xlabel = row == 2 ? L"z\ [\mathrm{pc}]" : "",
                ylabel = col == 1 ? L"x\ [\mathrm{pc}]" : "",
                xlabelsize = 38,
                ylabelsize = 38,
                xticklabelsize = 30,
                yticklabelsize = 30,
                xticks = ticks,
                yticks = ticks,
                aspect = DataAspect(),
                limits = (0, Float64(lbox_pc), 0, Float64(lbox_pc)),
                xautolimitmargin=(0.0f0, 0.0f0),
                yautolimitmargin=(0.0f0, 0.0f0),
            )
            hm = heatmap!(ax, x_edges, y_edges, maps[tag]; colormap=:magma, colorrange=pmax_colorrange)
            if row == 2
                zmask = fill(NaN, size(summaries[tag].zone_mask))
                zmask[summaries[tag].zone_mask] .= 1.0
                heatmap!(ax, x_edges, y_edges, zmask; colormap=zone_cmap, colorrange=(0, 1))
                segpts = Point2f[]
                for (p1, p2) in summaries[tag].overlay_segments
                    push!(segpts, p1)
                    push!(segpts, p2)
                end
                if !isempty(segpts)
                    linesegments!(ax, segpts; color=(:black, 0.85), linewidth=4.5)
                    linesegments!(ax, segpts; color=(:white, 0.98), linewidth=2.3)
                end
                text!(ax, 2.0, 47.5;
                    text="Nzones = $(summaries[tag].zone_count)",
                    align=(:left, :top), color=:white, fontsize=24)
            end
            if row == 1
                text!(ax, 47.0, 2.5;
                    text=tag == "NF" ? "NF" : string(tag, " pc"),
                    align=(:right, :bottom), color=:white, fontsize=32)
            end
            hidespines!(ax, :t, :r)
            col > 1 && hideydecorations!(ax, ticks=false, ticklabels=false, grid=false)
            row == 1 && hidexdecorations!(ax, ticks=false, ticklabels=false, grid=false)

            if col == length(labels) && row == 1
                Colorbar(fig[:, col + 1], hm;
                    label=L"P_{\max}\ [\mathrm{K}]",
                    labelsize=38,
                    ticklabelsize=30,
                    tickformat=InstrumentalEffect.cb_latex_ticks)
            end
        end
    end
    Label(fig[0, 1:5], "Zone-based test on depolarized regions"; fontsize=36)
    save(pdf_path, fig)
    save(png_path, fig)
end

function _outdir_for_los(los::AbstractString)
    return joinpath(@__DIR__, "..", "outputs", "zone_test", "LOS$(los)")
end

function measure_depolarized_zones()
    los = isempty(ARGS) ? nothing : String(ARGS[1])
    paths = _load_paths(CFG_PATH; los_override=los)
    outdir = _outdir_for_los(paths.los)
    mkpath(outdir)
    cfg = _build_cfg(paths)
    b1, b2, _ = InstrumentalEffect._project_bperp_maps(cfg)
    θB = InstrumentalEffect._orientation_map_bperp_from_b1_b2(b1, b2)

    n = size(b1, 1)
    m = size(b1, 2)
    x_pc = collect(range(paths.lbox_pc / (2n), paths.lbox_pc - paths.lbox_pc / (2n), length=n))
    y_pc = collect(range(paths.lbox_pc / (2m), paths.lbox_pc - paths.lbox_pc / (2m), length=m))

    labels = ["NF"]
    maps = Dict{String, Matrix{Float64}}()
    maps["NF"] = Float64.(InstrumentalEffect.read_FITS(cfg.Pmax_nofilter_path))
    for Llarge in reverse(sort(paths.llarge_list))
        tag = InstrumentalEffect._filter_tag(Llarge)
        ppath = joinpath(paths.base_out, tag, "RMSynthesis", "Pphi_max.fits")
        label = _map_label(paths.lbox_pc * Llarge / n)
        push!(labels, label)
        maps[label] = Float64.(InstrumentalEffect.read_FITS(ppath))
    end

    summaries = Dict{String, NamedTuple}()
    for tag in labels
        summaries[tag] = _zone_summary(maps[tag], θB, x_pc, y_pc)
    end

    out_csv = joinpath(outdir, "depolarized_zone_summary.csv")
    out_pdf = joinpath(outdir, "depolarized_zone_test.pdf")
    out_png = joinpath(outdir, "depolarized_zone_test.png")
    _write_csv(out_csv, labels, summaries)
    _plot_zone_test(out_pdf, out_png, labels, maps, summaries, paths.lbox_pc)

    @info "Saved zone-based test" pdf_path=out_pdf csv_path=out_csv
end

if abspath(PROGRAM_FILE) == @__FILE__
    measure_depolarized_zones()
end
