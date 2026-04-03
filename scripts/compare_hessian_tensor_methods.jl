using TOML
using CairoMakie
using LaTeXStrings
using Statistics
using Random
using DelimitedFiles

include(joinpath(@__DIR__, "..", "src", "code", "instrumental_effect", "InstrumentalEffect.jl"))
using .InstrumentalEffect

const CFG_PATH = joinpath(@__DIR__, "..", "config", "default.toml")
const OUTDIR = joinpath(@__DIR__, "..", "outputs", "method_comparison")

function _load_cfg()
    cfg = TOML.parsefile(CFG_PATH)
    sim_root = String(cfg["paths"]["simulations_root"])
    sim = String(cfg["simulation"]["name"])
    los = String(cfg["simulation"]["los"])
    llarge_list = Float64.(cfg["tasks"]["instrumental_effect"]["llarge_list"])

    ENV["SIMULATIONS_ROOT"] = sim_root
    ENV["SIMU_NAME"] = sim
    ENV["SIMU_LOS"] = los
    icfg = InstrumentalEffect.load_config_defaults()
    return cfg, icfg, sim_root, sim, los, llarge_list
end

function _filter_dir(llarge::Real)
    return "HardBandPass_remove_L0_to_1pc_and_$(round(Int, llarge))to50pc"
end

function _label_for_filter(llarge::Real, lbox_pc::Real, n::Int)
    return string(round(Int, lbox_pc * float(llarge) / n))
end

function _component_count(mask::AbstractMatrix{Bool})
    n, m = size(mask)
    seen = falses(n, m)
    count = 0
    neigh = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
    stack = CartesianIndex{2}[]
    @inbounds for j in 1:m, i in 1:n
        mask[i, j] || continue
        seen[i, j] && continue
        count += 1
        empty!(stack)
        push!(stack, CartesianIndex(i, j))
        seen[i, j] = true
        while !isempty(stack)
            I = pop!(stack)
            for (di, dj) in neigh
                ii = I[1] + di
                jj = I[2] + dj
                if 1 <= ii <= n && 1 <= jj <= m && mask[ii, jj] && !seen[ii, jj]
                    seen[ii, jj] = true
                    push!(stack, CartesianIndex(ii, jj))
                end
            end
        end
    end
    return count
end

function _random_filament_map(n::Int, m::Int; seed::Int=42, nxcell::Int=10, nycell::Int=10)
    rng = MersenneTwister(seed)
    A = ones(Float64, n, m)
    xcenters = collect(range(16.0, n - 16.0; length=nxcell))
    ycenters = collect(range(16.0, m - 16.0; length=nycell))
    dxcell = nxcell > 1 ? step(range(16.0, n - 16.0; length=nxcell)) : 0.0
    dycell = nycell > 1 ? step(range(16.0, m - 16.0; length=nycell)) : 0.0

    injected = NamedTuple{(:xc, :yc, :theta_deg)}[]
    for yc0 in ycenters, xc0 in xcenters
        xc = xc0 + (rand(rng) - 0.5) * 0.25 * dxcell
        yc = yc0 + (rand(rng) - 0.5) * 0.25 * dycell
        θ = rand(rng) * π
        push!(injected, (xc=xc, yc=yc, theta_deg=rad2deg(θ)))

        σpar = 9.5 + 2.0 * rand(rng)
        σperp = 1.7 + 0.4 * rand(rng)
        amp = 0.32 + 0.05 * rand(rng)
        cθ = cos(θ)
        sθ = sin(θ)

        @inbounds for j in 1:m, i in 1:n
            dx = i - xc
            dy = j - yc
            dpar = dx * cθ + dy * sθ
            dperp = -dx * sθ + dy * cθ
            A[i, j] -= amp * exp(-0.5 * ((dpar / σpar)^2 + (dperp / σperp)^2))
        end
    end

    A .+= 0.01 .* randn(rng, n, m)
    A = InstrumentalEffect._box_mean(A, 1)
    Amin = minimum(A)
    Amax = maximum(A)
    A = (A .- Amin) ./ (Amax - Amin)
    return A, injected
end

function _recover_angles_near_injected_filaments(θmap::AbstractMatrix, mask::AbstractMatrix{Bool}, injected;
        radius_pix::Int=8)
    recovered = Float64[]
    input = Float64[]
    for fil in injected
        i0 = round(Int, fil.xc)
        j0 = round(Int, fil.yc)
        vals = Float64[]
        for j in max(1, j0 - radius_pix):min(size(mask, 2), j0 + radius_pix),
            i in max(1, i0 - radius_pix):min(size(mask, 1), i0 + radius_pix)
            mask[i, j] || continue
            θ = θmap[i, j]
            isfinite(θ) && push!(vals, Float64(θ))
        end
        isempty(vals) && continue
        c2 = mean(cos.(2 .* vals))
        s2 = mean(sin.(2 .* vals))
        θrec = rad2deg(mod(0.5 * atan(s2, c2), π))
        push!(recovered, θrec)
        push!(input, fil.theta_deg)
    end
    return input, recovered
end

@inline function _delta_deg_pi_periodic(a_deg::Real, b_deg::Real)
    Δ = mod((float(a_deg) - float(b_deg)) + 90.0, 180.0) - 90.0
    return abs(Δ)
end

function _random_bias_metrics(input_angles_deg, recovered_angles_deg)
    isempty(recovered_angles_deg) && return (nrec=0, resultant=NaN, median_error_deg=NaN)
    θ = deg2rad.(Float64.(recovered_angles_deg))
    r = hypot(mean(cos.(2 .* θ)), mean(sin.(2 .* θ)))
    errs = [_delta_deg_pi_periodic(a, b) for (a, b) in zip(recovered_angles_deg, input_angles_deg)]
    return (nrec=length(recovered_angles_deg), resultant=r, median_error_deg=median(errs))
end

function _hessian_orientation_map(Pmap::AbstractMatrix, Δx::Real, Δy::Real)
    score_raw = InstrumentalEffect._canal_score_map(Pmap)
    score = InstrumentalEffect._box_mean(score_raw, InstrumentalEffect._CANAL_RIDGE_SMOOTH_RADIUS_PIX)
    thr = InstrumentalEffect._finite_quantile(score, InstrumentalEffect._CANAL_SCORE_QUANTILE)
    fxx, fxy, fyy, gyx, gyy = InstrumentalEffect._hessian_components(score, Δx, Δy)

    θcanal = fill(NaN, size(score))
    mask = falses(size(score))

    @inbounds for j in axes(score, 2), i in axes(score, 1)
        Dij = score[i, j]
        isfinite(Dij) || continue
        Dij >= thr || continue

        λ1, λ2, v1x, v1y, v2x, v2y = InstrumentalEffect._symmetric_2x2_eigensystem(
            fxx[i, j], fxy[i, j], fyy[i, j]
        )
        λ1 < 0 || continue
        anis = abs(λ1) / max(abs(λ2), 1e-6)
        anis >= InstrumentalEffect._CANAL_RIDGE_ANIS_MIN || continue

        sp = InstrumentalEffect._sample_bilinear(score, i + v1x, j + v1y)
        sm = InstrumentalEffect._sample_bilinear(score, i - v1x, j - v1y)
        isfinite(sp) && isfinite(sm) || continue
        (Dij >= sp && Dij >= sm) || continue

        θcanal[i, j] = InstrumentalEffect._wrap_orientation_pi(atan(v2y, v2x))
        mask[i, j] = true
    end

    if count(mask) == 0
        fallback_mask = isfinite.(score) .& (score .>= thr)
        @inbounds for j in axes(score, 2), i in axes(score, 1)
            fallback_mask[i, j] || continue
            tx = -gyy[i, j]
            ty = gyx[i, j]
            if isfinite(tx) && isfinite(ty) && (tx != 0 || ty != 0)
                θcanal[i, j] = InstrumentalEffect._wrap_orientation_pi(atan(ty, tx))
                mask[i, j] = true
            end
        end
    end
    return θcanal, mask, score
end

function _tensor_orientation_map(Pmap::AbstractMatrix, Δx::Real, Δy::Real)
    θcanal, mask, score, _ = InstrumentalEffect._canal_orientation_map(Pmap, Δx, Δy)
    return θcanal, mask, score
end

function _measure_method(Pmap::AbstractMatrix, θB::AbstractMatrix, Δx::Real, Δy::Real, method::Symbol)
    θcanal, mask, _ = method == :tensor ? _tensor_orientation_map(Pmap, Δx, Δy) : _hessian_orientation_map(Pmap, Δx, Δy)
    valid = mask .& isfinite.(θB)
    deltas = InstrumentalEffect._alignment_deltas_deg(θcanal, θB, valid)
    stats = InstrumentalEffect._alignment_stats_from_deltas(deltas)
    return (
        npix=count(mask),
        ncomp=_component_count(mask),
        median_abs_delta_deg=stats.median_abs_delta_deg,
        mean_abs_delta_deg=stats.mean_abs_delta_deg,
        frac_parallel_15deg=stats.frac_parallel_tol,
        θcanal=θcanal,
        mask=mask,
    )
end

function _plot_results(rows, random_payload, outbase::AbstractString)
    order = ["NF", "30", "20", "15", "10"]
    xpos = 1:length(order)
    xticklabels = [LaTeXString("\\mathrm{$tag}") for tag in order]
    method_styles = Dict(
        :tensor => (label=LaTeXString("\\mathrm{structure\\ tensor}"), color=:dodgerblue3),
        :hessian => (label=LaTeXString("\\mathrm{Hessian}"), color=:red3),
    )

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1600, 560), figure_padding=(18, 18, 10, 8))

        ax1 = Axis(fig[1, 1],
            xlabel=LaTeXString("L_{\\mathrm{large}}\\,[\\mathrm{pc}]"),
            ylabel=LaTeXString("\\mathrm{med}\\,|\\Delta\\theta(B_\\perp)|\\,[^\\circ]"),
            xlabelsize=26, ylabelsize=26, xticklabelsize=22, yticklabelsize=22,
            xminorticksvisible=false, yminorticksvisible=false,
            xminorgridvisible=false, yminorgridvisible=false,
        )
        ax2 = Axis(fig[1, 2],
            xlabel=LaTeXString("L_{\\mathrm{large}}\\,[\\mathrm{pc}]"),
            ylabel=LaTeXString("N_{\\mathrm{comp}}"),
            xlabelsize=26, ylabelsize=26, xticklabelsize=22, yticklabelsize=22,
            xminorticksvisible=false, yminorticksvisible=false,
            xminorgridvisible=false, yminorgridvisible=false,
        )
        ax3 = Axis(fig[1, 3],
            xlabel=LaTeXString("\\theta\\,[^\\circ]"),
            ylabel=LaTeXString("p(\\theta)"),
            xlabelsize=26, ylabelsize=26, xticklabelsize=22, yticklabelsize=22,
            xminorticksvisible=false, yminorticksvisible=false,
            xminorgridvisible=false, yminorgridvisible=false,
        )

        for method in (:tensor, :hessian)
            st = method_styles[method]
            y1 = [rows[(tag, method)].median_abs_delta_deg for tag in order]
            y2 = [rows[(tag, method)].ncomp for tag in order]
            lines!(ax1, xpos, y1; color=st.color, linewidth=3)
            scatter!(ax1, xpos, y1; color=st.color, markersize=16, strokecolor=:white, strokewidth=1.5, label=st.label)
            lines!(ax2, xpos, y2; color=st.color, linewidth=3)
            scatter!(ax2, xpos, y2; color=st.color, markersize=16, strokecolor=:white, strokewidth=1.5, label=st.label)
        end
        ax1.xticks = (collect(xpos), xticklabels)
        ax2.xticks = (collect(xpos), xticklabels)
        axislegend(ax1; position=:rt, framevisible=true, labelsize=20)
        axislegend(ax2; position=:lt, framevisible=true, labelsize=20)

        bins = 0:10:180
        hist!(ax3, random_payload.input; bins=bins, normalization=:pdf, color=(:lightgray, 0.7), strokecolor=:transparent)
        hist!(ax3, random_payload.tensor; bins=bins, normalization=:pdf, color=(:dodgerblue3, 0.45), strokecolor=:dodgerblue3, strokewidth=1.0, label=LaTeXString("\\mathrm{structure\\ tensor}"))
        hist!(ax3, random_payload.hessian; bins=bins, normalization=:pdf, color=(:red3, 0.35), strokecolor=:red3, strokewidth=1.0, label=LaTeXString("\\mathrm{Hessian}"))
        hlines!(ax3, [1 / 180]; color=:black, linestyle=:dash, linewidth=2)
        ax3.xticks = (0:30:180, [LaTeXString("\\mathrm{$v}") for v in 0:30:180])
        axislegend(ax3; position=:rt, framevisible=true, labelsize=20)

        text!(ax3, 0.04, 0.96;
            text="R_tensor = $(round(random_payload.tensor_resultant; digits=3))\nR_Hessian = $(round(random_payload.hessian_resultant; digits=3))",
            space=:relative,
            align=(:left, :top),
            fontsize=20,
            color=:black,
        )

        save(outbase * ".pdf", fig)
        save(outbase * ".png", fig)
    end
end

function main()
    cfg_raw, cfg, sim_root, sim, los, llarge_list = _load_cfg()
    Δx = cfg.Lbox_pc / cfg.n
    Δy = cfg.Lbox_pc / cfg.m

    p0 = Float64.(InstrumentalEffect.read_FITS(joinpath(sim_root, sim, los, "Synchrotron", "WithFaraday", "Pmax.fits")))
    base_out = joinpath(@__DIR__, "..", "outputs", "instrumental", sim, "LOS$(los)")

    filtered = Dict{String, Matrix{Float64}}()
    for L in llarge_list
        tag = _label_for_filter(L, cfg.Lbox_pc, cfg.n)
        path = joinpath(base_out, _filter_dir(L), "RMSynthesis", "Pphi_max.fits")
        filtered[tag] = Float64.(InstrumentalEffect.read_FITS(path))
    end

    b1, b2, _ = InstrumentalEffect._project_bperp_maps(cfg)
    θB = InstrumentalEffect._orientation_map_bperp_from_b1_b2(b1, b2)

    rows = Dict{Tuple{String, Symbol}, NamedTuple}()
    maps = Dict("NF" => p0)
    merge!(maps, filtered)
    for tag in ["NF", "30", "20", "15", "10"]
        Pmap = maps[tag]
        rows[(tag, :tensor)] = _measure_method(Pmap, θB, Δx, Δy, :tensor)
        rows[(tag, :hessian)] = _measure_method(Pmap, θB, Δx, Δy, :hessian)
    end

    randmap, injected = _random_filament_map(cfg.n, cfg.m)
    θt, mt, _ = _tensor_orientation_map(randmap, 1.0, 1.0)
    θh, mh, _ = _hessian_orientation_map(randmap, 1.0, 1.0)
    input_t, rec_t = _recover_angles_near_injected_filaments(θt, mt, injected)
    input_h, rec_h = _recover_angles_near_injected_filaments(θh, mh, injected)
    bias_t = _random_bias_metrics(input_t, rec_t)
    bias_h = _random_bias_metrics(input_h, rec_h)

    mkpath(OUTDIR)
    csv_path = joinpath(OUTDIR, "hessian_tensor_comparison.csv")
    open(csv_path, "w") do io
        println(io, "map_tag,method,npix,ncomp,median_abs_delta_deg,mean_abs_delta_deg,frac_parallel_15deg")
        for tag in ["NF", "30", "20", "15", "10"], method in (:tensor, :hessian)
            r = rows[(tag, method)]
            println(io, string(tag, ",", method, ",", r.npix, ",", r.ncomp, ",", r.median_abs_delta_deg, ",", r.mean_abs_delta_deg, ",", r.frac_parallel_15deg))
        end
        println(io, "random,tensor,$(bias_t.nrec),,,$(bias_t.resultant),$(bias_t.median_error_deg)")
        println(io, "random,hessian,$(bias_h.nrec),,,$(bias_h.resultant),$(bias_h.median_error_deg)")
    end

    fig_base = joinpath(OUTDIR, "hessian_tensor_comparison")
    _plot_results(rows, (
        input=input_t,
        tensor=rec_t,
        hessian=rec_h,
        tensor_resultant=bias_t.resultant,
        hessian_resultant=bias_h.resultant,
    ), fig_base)

    println("CSV: " * csv_path)
    println("FIG: " * fig_base * ".pdf")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
