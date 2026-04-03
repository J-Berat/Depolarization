using TOML
using Statistics
using Printf

include(joinpath(@__DIR__, "..", "src", "code", "instrumental_effect", "InstrumentalEffect.jl"))

using .InstrumentalEffect

const CFG_PATH = joinpath(@__DIR__, "..", "config", "default.toml")
const OUTDIR = joinpath(@__DIR__, "..", "outputs", "window_sweep_tensor")
const ANGLE_TOL_DEG = 15.0

function _load_paths(cfg_path::AbstractString)
    cfg = TOML.parsefile(cfg_path)
    sim_root = String(cfg["paths"]["simulations_root"])
    sim = String(cfg["simulation"]["name"])
    los = String(cfg["simulation"]["los"])
    lbox_pc = Float64(get(cfg["tasks"]["canal_metrics"], "lbox_pc", 50.0))
    llarge_list = Float64.(cfg["tasks"]["instrumental_effect"]["llarge_list"])
    base_out = joinpath(@__DIR__, "..", "outputs", "instrumental", sim, "LOS$(los)")
    return (; sim_root, sim, los, lbox_pc, llarge_list, base_out)
end

function _finite_quantile(A, q::Real)
    vals = vec(Float64.(A))
    vals = vals[isfinite.(vals)]
    isempty(vals) && return NaN
    return quantile(vals, q)
end

function _box_mean(A::AbstractMatrix, radius::Int)
    radius >= 0 || error("radius must be >= 0")
    n, m = size(A)
    S = zeros(Float64, n + 1, m + 1)
    @inbounds for j in 1:m
        rowsum = 0.0
        for i in 1:n
            rowsum += Float64(A[i, j])
            S[i + 1, j + 1] = S[i + 1, j] + rowsum
        end
    end
    out = Matrix{Float64}(undef, n, m)
    @inbounds for j in 1:m, i in 1:n
        i1 = max(1, i - radius)
        i2 = min(n, i + radius)
        j1 = max(1, j - radius)
        j2 = min(m, j + radius)
        area = (i2 - i1 + 1) * (j2 - j1 + 1)
        s = S[i2 + 1, j2 + 1] - S[i1, j2 + 1] - S[i2 + 1, j1] + S[i1, j1]
        out[i, j] = s / area
    end
    return out
end

function _gradients_central(A::AbstractMatrix, Δx::Real, Δy::Real)
    n, m = size(A)
    gx = Matrix{Float64}(undef, n, m)
    gy = Matrix{Float64}(undef, n, m)
    invdx = inv(float(Δx))
    invdy = inv(float(Δy))
    inv2dx = inv(2.0 * float(Δx))
    inv2dy = inv(2.0 * float(Δy))
    @inbounds for j in 1:m
        gx[1, j] = (A[2, j] - A[1, j]) * invdx
        for i in 2:(n - 1)
            gx[i, j] = (A[i + 1, j] - A[i - 1, j]) * inv2dx
        end
        gx[n, j] = (A[n, j] - A[n - 1, j]) * invdx
    end
    @inbounds for i in 1:n
        gy[i, 1] = (A[i, 2] - A[i, 1]) * invdy
        for j in 2:(m - 1)
            gy[i, j] = (A[i, j + 1] - A[i, j - 1]) * inv2dy
        end
        gy[i, m] = (A[i, m] - A[i, m - 1]) * invdy
    end
    return gx, gy
end

function _symmetric_2x2_eigensystem(a::Real, b::Real, c::Real)
    tr = float(a) + float(c)
    disc = sqrt(max((float(a) - float(c))^2 + 4.0 * float(b)^2, 0.0))
    λ1 = 0.5 * (tr - disc)
    λ2 = 0.5 * (tr + disc)
    v1x = float(b)
    v1y = λ1 - float(a)
    if abs(v1x) + abs(v1y) < 1e-12
        v1x = λ2 - float(c)
        v1y = float(b)
    end
    v2x = float(b)
    v2y = λ2 - float(a)
    if abs(v2x) + abs(v2y) < 1e-12
        v2x = λ1 - float(c)
        v2y = float(b)
    end
    return λ1, λ2, v1x, v1y, v2x, v2y
end

@inline _wrap_orientation_pi(θ::Real) = mod(float(θ), π)

@inline function _line_orientation_delta_deg(θa::Real, θb::Real)
    Δ = mod((float(θa) - float(θb)) + (π / 2), π) - (π / 2)
    return rad2deg(Δ)
end

function _orientation_map_from_components(vx::AbstractMatrix, vy::AbstractMatrix)
    θ = fill(NaN, size(vx))
    @inbounds for j in axes(vx, 2), i in axes(vx, 1)
        ax = vx[i, j]
        ay = vy[i, j]
        if isfinite(ax) && isfinite(ay) && (ax != 0 || ay != 0)
            θ[i, j] = _wrap_orientation_pi(atan(ay, ax))
        end
    end
    return θ
end

function _project_cube_mean(cube, los::String, target_size::Tuple{Int,Int})
    dim = los == "x" ? 1 : (los == "y" ? 2 : 3)
    proj = dropdims(mean(Float64.(cube); dims=dim), dims=dim)
    if size(proj) == target_size
        return Matrix{Float64}(proj)
    end
    proj_t = permutedims(proj, (2, 1))
    size(proj_t) == target_size || error("Projected map size mismatch")
    return Matrix{Float64}(proj_t)
end

function _project_hypot_mean(B1_cube, B2_cube, los::String, target_size::Tuple{Int,Int})
    nx, ny, nz = size(B1_cube)
    if los == "x"
        out = Matrix{Float64}(undef, ny, nz)
        @inbounds for j in 1:nz, i in 1:ny
            acc = 0.0
            for k in 1:nx
                acc += hypot(B1_cube[k, i, j], B2_cube[k, i, j])
            end
            out[i, j] = acc / nx
        end
    elseif los == "y"
        out = Matrix{Float64}(undef, nx, nz)
        @inbounds for j in 1:nz, i in 1:nx
            acc = 0.0
            for k in 1:ny
                acc += hypot(B1_cube[i, k, j], B2_cube[i, k, j])
            end
            out[i, j] = acc / ny
        end
    else
        out = Matrix{Float64}(undef, nx, ny)
        @inbounds for j in 1:ny, i in 1:nx
            acc = 0.0
            for k in 1:nz
                acc += hypot(B1_cube[i, j, k], B2_cube[i, j, k])
            end
            out[i, j] = acc / nz
        end
    end
    if size(out) == target_size
        return out
    end
    out_t = permutedims(out, (2, 1))
    size(out_t) == target_size || error("Projected hypot size mismatch")
    return out_t
end

function _project_bperp_maps(cfg)
    Bx = Float64.(InstrumentalEffect.read_FITS(cfg.Bx_in))
    By = Float64.(InstrumentalEffect.read_FITS(cfg.By_in))
    Bz = Float64.(InstrumentalEffect.read_FITS(cfg.Bz_in))
    target = (cfg.n, cfg.m)
    if cfg.los == "x"
        B1_cube = By
        B2_cube = Bz
    elseif cfg.los == "y"
        B1_cube = Bx
        B2_cube = Bz
    else
        B1_cube = Bx
        B2_cube = By
    end
    b1 = _project_cube_mean(B1_cube, cfg.los, target)
    b2 = _project_cube_mean(B2_cube, cfg.los, target)
    bperp = _project_hypot_mean(B1_cube, B2_cube, cfg.los, target)
    return b1, b2, bperp
end

function _canal_orientation_map_param(Pmap::AbstractMatrix, Δx::Real, Δy::Real;
        local_radius_pix::Int, ridge_smooth_radius_pix::Int,
        tensor_smooth_radius_pix::Int, score_quantile::Real=0.95, anis_min::Real=3.0)
    P = Float64.(Pmap)
    valid = isfinite.(P)
    fillval = _finite_quantile(P, 0.50)
    Pfill = copy(P)
    Pfill[.!valid] .= fillval
    Plocal = _box_mean(Pfill, local_radius_pix)
    score_raw = Plocal .- Pfill
    score_raw[.!valid] .= NaN
    score = _box_mean(score_raw, ridge_smooth_radius_pix)
    thr = _finite_quantile(score, score_quantile)
    gx, gy = _gradients_central(score, Δx, Δy)
    Jxx = _box_mean(gx .^ 2, tensor_smooth_radius_pix)
    Jxy = _box_mean(gx .* gy, tensor_smooth_radius_pix)
    Jyy = _box_mean(gy .^ 2, tensor_smooth_radius_pix)
    θcanal = fill(NaN, size(score))
    mask = falses(size(score))
    @inbounds for j in axes(score, 2), i in axes(score, 1)
        Dij = score[i, j]
        isfinite(Dij) || continue
        Dij >= thr || continue
        λtan, λnorm, tvecx, tvecy, nvecx, nvecy = _symmetric_2x2_eigensystem(Jxx[i, j], Jxy[i, j], Jyy[i, j])
        λnorm > 0 || continue
        anis = λnorm / max(λtan, 1e-6)
        anis >= anis_min || continue
        abs(tvecx) + abs(tvecy) > 1e-12 || continue
        abs(nvecx) + abs(nvecy) > 1e-12 || continue
        sx = round(Int, clamp(sign(nvecx), -1, 1))
        sy = round(Int, clamp(sign(nvecy), -1, 1))
        (sx != 0 || sy != 0) || continue
        i1 = clamp(i - sx, first(axes(score, 1)), last(axes(score, 1)))
        i2 = clamp(i + sx, first(axes(score, 1)), last(axes(score, 1)))
        j1 = clamp(j - sy, first(axes(score, 2)), last(axes(score, 2)))
        j2 = clamp(j + sy, first(axes(score, 2)), last(axes(score, 2)))
        (Dij >= score[i1, j1] && Dij >= score[i2, j2]) || continue
        θcanal[i, j] = _wrap_orientation_pi(atan(tvecy, tvecx))
        mask[i, j] = true
    end
    return θcanal, mask, score
end

function _alignment_deltas_deg(θcanal::AbstractMatrix, θref::AbstractMatrix, mask::AbstractMatrix)
    deltas = Float64[]
    @inbounds for j in axes(θcanal, 2), i in axes(θcanal, 1)
        if mask[i, j]
            a = θcanal[i, j]
            b = θref[i, j]
            if isfinite(a) && isfinite(b)
                push!(deltas, _line_orientation_delta_deg(a, b))
            end
        end
    end
    return deltas
end

function _alignment_stats_from_deltas(deltas::AbstractVector{<:Real})
    absd = abs.(Float64.(deltas))
    perp = abs.(absd .- 90.0)
    n = length(absd)
    n == 0 && return (npix=0, mean_abs_delta_deg=NaN, median_abs_delta_deg=NaN, frac_parallel_tol=NaN, frac_perp_tol=NaN, median_perp_offset_deg=NaN)
    return (
        npix=n,
        mean_abs_delta_deg=mean(absd),
        median_abs_delta_deg=median(absd),
        frac_parallel_tol=mean(absd .<= ANGLE_TOL_DEG),
        frac_perp_tol=mean(perp .<= ANGLE_TOL_DEG),
        median_perp_offset_deg=median(perp),
    )
end

function _count_true_runs(v::AbstractVector{Bool})
    nruns = 0
    in_run = false
    @inbounds for b in v
        if b
            if !in_run
                nruns += 1
                in_run = true
            end
        else
            in_run = false
        end
    end
    return nruns
end

function _rotate_mask_nearest(mask::AbstractMatrix{Bool}, α::Real)
    nx, ny = size(mask)
    cx = 0.5 * (nx + 1)
    cy = 0.5 * (ny + 1)
    cα = cos(float(α))
    sα = sin(float(α))
    out = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        x = i - cx
        y = j - cy
        xs = cα * x + sα * y
        ys = -sα * x + cα * y
        ii = round(Int, xs + cx)
        jj = round(Int, ys + cy)
        if 1 <= ii <= nx && 1 <= jj <= ny
            out[i, j] = mask[ii, jj]
        end
    end
    return out
end

function _crossings_per_normal_cut(mask::AbstractMatrix{Bool}, θdom::Real)
    rot = _rotate_mask_nearest(mask, -float(θdom))
    counts = Float64[]
    @inbounds for i in axes(rot, 1)
        push!(counts, _count_true_runs(view(rot, i, :)))
    end
    return counts
end

function _dominant_line_orientation(θvals::AbstractVector{<:Real})
    isempty(θvals) && return NaN
    c2 = mean(cos.(2.0 .* θvals))
    s2 = mean(sin.(2.0 .* θvals))
    θ = 0.5 * atan(s2, c2)
    return mod(θ, π)
end

function analyze_case(Pmap::AbstractMatrix, b1::AbstractMatrix, b2::AbstractMatrix, bperp::AbstractMatrix;
        local_radius_pix::Int, ridge_smooth_radius_pix::Int, tensor_smooth_radius_pix::Int, Δ::Real)
    θcanal, mask, _ = _canal_orientation_map_param(Pmap, Δ, Δ;
        local_radius_pix=local_radius_pix, ridge_smooth_radius_pix=ridge_smooth_radius_pix,
        tensor_smooth_radius_pix=tensor_smooth_radius_pix)
    θB = _orientation_map_from_components(b2, b1)
    mask_B = mask .& isfinite.(θB) .& isfinite.(bperp) .& (bperp .> 0)
    deltas = _alignment_deltas_deg(θcanal, θB, mask_B)
    bstats = _alignment_stats_from_deltas(deltas)
    θvals = Float64[]
    @inbounds for j in axes(mask, 2), i in axes(mask, 1)
        mask[i, j] || continue
        θ = θcanal[i, j]
        isfinite(θ) || continue
        push!(θvals, θ)
    end
    θdom = _dominant_line_orientation(θvals)
    delta_map = fill(NaN, size(mask))
    deltas_par = Float64[]
    @inbounds for j in axes(mask, 2), i in axes(mask, 1)
        mask[i, j] || continue
        θ = θcanal[i, j]
        isfinite(θ) || continue
        δ = _line_orientation_delta_deg(θ, θdom)
        push!(deltas_par, δ)
        delta_map[i, j] = δ
    end
    parallel_mask = mask .& isfinite.(delta_map) .& (abs.(delta_map) .<= ANGLE_TOL_DEG)
    counts = _crossings_per_normal_cut(parallel_mask, θdom)
    pstats = (
        npix=count(mask),
        mean_parallel_deg=mean(deltas_par),
        sigma_parallel_deg=std(deltas_par),
        median_abs_parallel_deg=median(abs.(deltas_par)),
        frac_parallel_family_15deg=mean(abs.(deltas_par) .<= ANGLE_TOL_DEG),
        mean_crossings=mean(counts),
        median_crossings=median(counts),
    )
    return merge(bstats, pstats, (theta_dom_deg=rad2deg(θdom),))
end

function main()
    paths = _load_paths(CFG_PATH)
    cfg = InstrumentalEffect.InstrumentalConfig(
        Q_in = joinpath(paths.sim_root, paths.sim, paths.los, "Synchrotron", "WithFaraday", "Qnu.fits"),
        U_in = joinpath(paths.sim_root, paths.sim, paths.los, "Synchrotron", "WithFaraday", "Unu.fits"),
        Pmax_nofilter_path = joinpath(paths.sim_root, paths.sim, paths.los, "Synchrotron", "WithFaraday", "Pmax.fits"),
        Q_in_phi = joinpath(paths.sim_root, paths.sim, paths.los, "Synchrotron", "WithFaraday", "realFDF.fits"),
        U_in_phi = joinpath(paths.sim_root, paths.sim, paths.los, "Synchrotron", "WithFaraday", "imagFDF.fits"),
        base_out = joinpath(paths.base_out),
        Bx_in = joinpath(paths.sim_root, paths.sim, "Bx.fits"),
        By_in = joinpath(paths.sim_root, paths.sim, "By.fits"),
        Bz_in = joinpath(paths.sim_root, paths.sim, "Bz.fits"),
        dens_in = joinpath(paths.sim_root, paths.sim, "density.fits"),
        los = paths.los,
    )

    b1, b2, bperp = _project_bperp_maps(cfg)
    maps = Dict{String, Matrix{Float64}}()
    maps["NF"] = Float64.(InstrumentalEffect.read_FITS(cfg.Pmax_nofilter_path))
    for Llarge in paths.llarge_list
        tag = InstrumentalEffect._filter_tag(Llarge)
        ppath = joinpath(cfg.base_out, tag, "RMSynthesis", "Pphi_max.fits")
        maps[string(round(Int, cfg.Lbox_pc * Llarge / cfg.n))] = Float64.(InstrumentalEffect.read_FITS(ppath))
    end

    local_radii = [1, 2, 3, 4]
    ridge_radii = [5, 7, 9]
    tensor_radii = [2, 3, 4, 5]

    mkpath(OUTDIR)
    out_csv = joinpath(OUTDIR, "canal_window_tensor_sweep.csv")
    open(out_csv, "w") do io
        println(io, "map_tag,local_window,ridge_window,tensor_window,mean_abs_delta_deg,median_abs_delta_deg,frac_parallel_15deg,npix,mean_crossings,median_crossings,sigma_parallel_deg,median_abs_parallel_deg,frac_parallel_family_15deg,theta_dom_deg")
        for (map_tag, Pmap) in maps
            for rl in local_radii
                for rr in ridge_radii
                    for rt in tensor_radii
                        stats = analyze_case(Pmap, b1, b2, bperp;
                            local_radius_pix=rl, ridge_smooth_radius_pix=rr,
                            tensor_smooth_radius_pix=rt, Δ=paths.lbox_pc / 256.0)
                        @printf(io, "%s,%d,%d,%d,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                            map_tag, 2 * rl + 1, 2 * rr + 1, 2 * rt + 1,
                            stats.mean_abs_delta_deg, stats.median_abs_delta_deg, stats.frac_parallel_tol,
                            stats.npix, stats.mean_crossings, stats.median_crossings,
                            stats.sigma_parallel_deg, stats.median_abs_parallel_deg,
                            stats.frac_parallel_family_15deg, stats.theta_dom_deg)
                    end
                end
            end
        end
    end

    println("Wrote ", out_csv)
end

main()
