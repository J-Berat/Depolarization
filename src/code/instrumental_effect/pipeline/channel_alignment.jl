const _ALIGNMENT_ANGLE_TOL_DEG = 15.0
const _ALIGNMENT_MAP_ENTRY = NamedTuple{
    (:map_tag, :llarge_eff_pc, :Pmap),
    Tuple{String, Float64, Matrix{Float64}},
}
const _ALIGNMENT_ROW = NamedTuple{
    (:map_tag, :llarge_eff_pc, :reference, :npix,
     :mean_abs_delta_deg, :median_abs_delta_deg,
     :frac_parallel_15deg, :frac_perp_15deg, :median_perp_offset_deg),
    Tuple{String, Float64, String, Int, Float64, Float64, Float64, Float64, Float64},
}
const _ALIGNMENT_DELTA_ENTRY = NamedTuple{
    (:map_tag, :llarge_eff_pc, :reference, :deltas),
    Tuple{String, Float64, String, Vector{Float64}},
}

"""
    _finite_quantile(...)

Returns a robust quantile over finite values only.
"""
function _finite_quantile(A, q::Real)
    0.0 <= q <= 1.0 || error("q must be in [0,1], got $q")
    vals = vec(Float64.(A))
    vals = vals[isfinite.(vals)]
    isempty(vals) && return NaN
    return quantile(vals, q)
end

@inline _wrap_orientation_pi(θ::Real) = mod(float(θ), π)

@inline function _line_orientation_delta_deg(θa::Real, θb::Real)
    Δ = mod((float(θa) - float(θb)) + (π / 2), π) - (π / 2)
    return rad2deg(Δ)
end

"""
    _gradients_central(...)

Central (interior) + one-sided (boundary) finite differences.
"""
function _gradients_central(A::AbstractMatrix, Δx::Real, Δy::Real)
    n, m = size(A)
    n >= 2 || error("Need n>=2 to compute gradients, got n=$n")
    m >= 2 || error("Need m>=2 to compute gradients, got m=$m")
    isfinite(Δx) && Δx > 0 || error("Δx must be positive finite, got $Δx")
    isfinite(Δy) && Δy > 0 || error("Δy must be positive finite, got $Δy")

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

"""
    _hessian_components(...)

Builds Hessian components from finite differences.
"""
function _hessian_components(A::AbstractMatrix, Δx::Real, Δy::Real)
    gx, gy = _gradients_central(A, Δx, Δy)
    fxx, fxy1 = _gradients_central(gx, Δx, Δy)
    fxy2, fyy = _gradients_central(gy, Δx, Δy)
    fxy = 0.5 .* (fxy1 .+ fxy2)
    return fxx, fxy, fyy, gx, gy
end

"""
    _orientation_map_from_components(...)

Converts vector components to line orientation in `[0,π)`.
"""
function _orientation_map_from_components(vx::AbstractMatrix, vy::AbstractMatrix)
    size(vx) == size(vy) || error("Component size mismatch: vx=$(size(vx)) vy=$(size(vy))")
    n, m = size(vx)
    θ = fill(NaN, n, m)
    @inbounds for j in 1:m, i in 1:n
        ax = vx[i, j]
        ay = vy[i, j]
        if isfinite(ax) && isfinite(ay) && (ax != 0 || ay != 0)
            θ[i, j] = _wrap_orientation_pi(atan(ay, ax))
        end
    end
    return θ
end

"""
    _orientation_map_bperp_from_b1_b2(...)

Converts `(B1,B2)` to line orientation using the requested form `atan(B1, B2)`.
"""
function _orientation_map_bperp_from_b1_b2(B1::AbstractMatrix, B2::AbstractMatrix)
    size(B1) == size(B2) || error("Component size mismatch: B1=$(size(B1)) B2=$(size(B2))")
    n, m = size(B1)
    θ = fill(NaN, n, m)
    @inbounds for j in 1:m, i in 1:n
        b1 = B1[i, j]
        b2 = B2[i, j]
        if isfinite(b1) && isfinite(b2) && (b1 != 0 || b2 != 0)
            θ[i, j] = _wrap_orientation_pi(atan(b1, b2))
        end
    end
    return θ
end

"""
    _canal_orientation_map(...)

Estimates canal tangent orientation on low-intensity canals (`P <= 1st decile`)
using only first-order gradients (no Hessian).
"""
function _canal_orientation_map(Pmap::AbstractMatrix, Δx::Real, Δy::Real)
    P = Float64.(Pmap)
    base_valid = isfinite.(P)
    p10 = _finite_quantile(P, 0.10)
    mask = base_valid .& (P .<= p10)
    count(mask) == 0 && (mask .= base_valid)

    Px, Py = _gradients_central(P, Δx, Δy)
    θcanal = fill(NaN, size(P))
    @inbounds for j in axes(P, 2), i in axes(P, 1)
        if mask[i, j]
            gx = Px[i, j]
            gy = Py[i, j]
            if isfinite(gx) && isfinite(gy) && (gx != 0 || gy != 0)
                # Canal direction is tangent to local iso-intensity contour.
                θcanal[i, j] = _wrap_orientation_pi(atan(gy, gx) + (π / 2))
            end
        end
    end

    return θcanal, mask, Px, Py
end

"""
    _alignment_deltas_deg(...)

Collects per-pixel `Δθ` values (degrees) on a selected validity mask.
"""
function _alignment_deltas_deg(θcanal::AbstractMatrix, θref::AbstractMatrix, mask::AbstractMatrix)
    size(θcanal) == size(θref) || error("Orientation size mismatch: canal=$(size(θcanal)) ref=$(size(θref))")
    size(θcanal) == size(mask) || error("Mask size mismatch: canal=$(size(θcanal)) mask=$(size(mask))")

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

"""
    _alignment_stats_from_deltas(...)

Computes alignment/perpendicularity statistics on selected pixels.
"""
function _alignment_stats_from_deltas(deltas::AbstractVector{<:Real}; tol_deg::Real=_ALIGNMENT_ANGLE_TOL_DEG)
    absd = abs.(Float64.(deltas))
    perp = abs.(absd .- 90.0)
    n = length(absd)

    n == 0 && return (
        npix=0,
        mean_abs_delta_deg=NaN,
        median_abs_delta_deg=NaN,
        frac_parallel_tol=NaN,
        frac_perp_tol=NaN,
        median_perp_offset_deg=NaN,
    )

    return (
        npix=n,
        mean_abs_delta_deg=mean(absd),
        median_abs_delta_deg=median(absd),
        frac_parallel_tol=mean(absd .<= tol_deg),
        frac_perp_tol=mean(perp .<= tol_deg),
        median_perp_offset_deg=median(perp),
    )
end

"""
    _alignment_stats(...)

Computes alignment/perpendicularity statistics on selected pixels.
"""
function _alignment_stats(θcanal::AbstractMatrix, θref::AbstractMatrix, mask::AbstractMatrix; tol_deg::Real=_ALIGNMENT_ANGLE_TOL_DEG)
    deltas = _alignment_deltas_deg(θcanal, θref, mask)
    return _alignment_stats_from_deltas(deltas; tol_deg=tol_deg)
end

"""
    _project_cube_mean(...)

LOS-average projection with map-shape reconciliation.
"""
function _project_cube_mean(cube, los::String, target_size::Tuple{Int,Int})
    ndims(cube) == 3 || error("Expected 3D cube for LOS projection, got ndims=$(ndims(cube)) size=$(size(cube))")
    los in ("x", "y", "z") || error("LOS must be x/y/z, got $los")

    dim = los == "x" ? 1 : (los == "y" ? 2 : 3)
    proj = dropdims(mean(Float64.(cube); dims=dim), dims=dim)

    if size(proj) == target_size
        return Matrix{Float64}(proj)
    end

    proj_t = permutedims(proj, (2, 1))
    if size(proj_t) == target_size
        return Matrix{Float64}(proj_t)
    end

    error("Projected map has size $(size(proj)) (or transposed $(size(proj_t))), expected $target_size")
end

"""
    _project_hypot_mean(...)

LOS-average projection of `hypot(B1, B2)` without allocating a full intermediate cube.
"""
function _project_hypot_mean(B1_cube, B2_cube, los::String, target_size::Tuple{Int,Int})
    size(B1_cube) == size(B2_cube) || error("B1/B2 size mismatch: B1=$(size(B1_cube)) B2=$(size(B2_cube))")
    ndims(B1_cube) == 3 || error("B1_cube must be 3D, got ndims=$(ndims(B1_cube))")

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
    elseif los == "z"
        out = Matrix{Float64}(undef, nx, ny)
        @inbounds for j in 1:ny, i in 1:nx
            acc = 0.0
            for k in 1:nz
                acc += hypot(B1_cube[i, j, k], B2_cube[i, j, k])
            end
            out[i, j] = acc / nz
        end
    else
        error("LOS must be x/y/z, got $los")
    end

    if size(out) == target_size
        return out
    end

    out_t = permutedims(out, (2, 1))
    if size(out_t) == target_size
        return out_t
    end
    error("Projected hypot map has size $(size(out)) (or transposed $(size(out_t))), expected $target_size")
end

"""
    _project_bperp_maps(...)

Projects magnetic field onto the sky plane according to LOS.
"""
function _project_bperp_maps(cfg::InstrumentalConfig)
    Bx = read_FITS(cfg.Bx_in)
    By = read_FITS(cfg.By_in)
    Bz = read_FITS(cfg.Bz_in)

    ndims(Bx) == 3 || error("Bx must be 3D, got size=$(size(Bx))")
    ndims(By) == 3 || error("By must be 3D, got size=$(size(By))")
    ndims(Bz) == 3 || error("Bz must be 3D, got size=$(size(Bz))")
    size(Bx) == size(By) || error("Bx/By shape mismatch: Bx=$(size(Bx)) By=$(size(By))")
    size(By) == size(Bz) || error("By/Bz shape mismatch: By=$(size(By)) Bz=$(size(Bz))")

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

    # Components in the sky plane.
    b1 = _project_cube_mean(B1_cube, cfg.los, target)
    b2 = _project_cube_mean(B2_cube, cfg.los, target)

    # Requested definition: B_perp = sqrt(B1^2 + B2^2), then LOS projection.
    bperp = _project_hypot_mean(B1_cube, B2_cube, cfg.los, target)

    return b1, b2, bperp
end

"""
    _phi_peak_map(...)

Computes per-pixel Faraday-depth peak location.
"""
function _phi_peak_map(Qphi_cube, Uphi_cube, PhiArray::AbstractVector, n::Int, m::Int)
    nphi = length(PhiArray)
    ndims(Qphi_cube) == 3 || error("realFDF must be 3D, got ndims=$(ndims(Qphi_cube)) size=$(size(Qphi_cube))")
    ndims(Uphi_cube) == 3 || error("imagFDF must be 3D, got ndims=$(ndims(Uphi_cube)) size=$(size(Uphi_cube))")

    qsz = size(Qphi_cube)
    usz = size(Uphi_cube)
    qsz == usz || error("realFDF/imagFDF shape mismatch: realFDF=$qsz imagFDF=$usz")

    layout = 0
    if qsz == (n, m, nphi)
        layout = 3
    elseif qsz == (nphi, n, m)
        layout = 1
    else
        error("FDF cubes must have shape (n,m,nphi)=($n,$m,$nphi) or (nphi,n,m)=($nphi,$n,$m), got size=$qsz")
    end

    phi_peak = Matrix{Float64}(undef, n, m)
    @inbounds for j in 1:m, i in 1:n
        kbest = 0
        vmax = -Inf
        for k in 1:nphi
            qv = layout == 3 ? Qphi_cube[i, j, k] : Qphi_cube[k, i, j]
            uv = layout == 3 ? Uphi_cube[i, j, k] : Uphi_cube[k, i, j]
            if isfinite(qv) && isfinite(uv)
                av = hypot(qv, uv)
                if av > vmax
                    vmax = av
                    kbest = k
                end
            end
        end
        phi_peak[i, j] = kbest == 0 ? NaN : Float64(PhiArray[kbest])
    end
    return phi_peak
end

"""
    _write_channel_alignment_csv(...)

Writes alignment summary table.
"""
function _write_channel_alignment_csv(path::AbstractString, rows::AbstractVector{<:NamedTuple})
    open(path, "w") do io
        println(io, "map_tag,llarge_eff_pc,reference,npix,mean_abs_delta_deg,median_abs_delta_deg,frac_parallel_15deg,frac_perp_15deg,median_perp_offset_deg")
        for r in rows
            println(io, @sprintf(
                "%s,%.6f,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f",
                r.map_tag,
                r.llarge_eff_pc,
                r.reference,
                r.npix,
                r.mean_abs_delta_deg,
                r.median_abs_delta_deg,
                r.frac_parallel_15deg,
                r.frac_perp_15deg,
                r.median_perp_offset_deg,
            ))
        end
    end
    return path
end

"""
    _kde_bandwidth(...)

Robust Gaussian KDE bandwidth (Silverman-like).
"""
function _kde_bandwidth(vals::AbstractVector{<:Real})
    n = length(vals)
    n <= 1 && return 5.0

    v = Float64.(vals)
    σ = std(v)
    q25 = quantile(v, 0.25)
    q75 = quantile(v, 0.75)
    iqr = q75 - q25
    s = min(σ, iqr / 1.34)
    h = 0.9 * s * n^(-0.2)

    if !(isfinite(h) && h > 0)
        h = isfinite(σ) && σ > 0 ? σ : 5.0
    end
    return clamp(h, 0.5, 20.0)
end

"""
    _kde_curve_density(...)

Builds a smooth KDE density curve `p(Δθ)`.
"""
function _kde_curve_density(vals::AbstractVector{<:Real};
                            xmin::Real=-90.0, xmax::Real=90.0,
                            ngrid::Int=721)
    ngrid >= 2 || error("ngrid must be >= 2, got $ngrid")

    xminf = float(xmin)
    xmaxf = float(xmax)
    xmaxf > xminf || error("xmax must be > xmin, got xmin=$xminf xmax=$xmaxf")

    x = collect(range(xminf, xmaxf; length=ngrid))
    y = zeros(Float64, ngrid)
    isempty(vals) && return x, y

    v = Float64.(vals)
    h = _kde_bandwidth(v)
    inv_norm = inv(length(v) * h * sqrt(2π))
    inv_h = inv(h)
    @inbounds for vk in v
        for i in eachindex(x)
            z = (x[i] - vk) * inv_h
            y[i] += exp(-0.5 * z^2)
        end
    end
    y .*= inv_norm
    return x, y
end

"""
    _plot_delta_theta_histograms(...)

Builds smooth `Δθ` KDE curves for `B_perp`,
with optional slider UI (controlled by `cfg.channel_alignment_pdf_plain`).
"""
function _plot_delta_theta_histograms(cfg::InstrumentalConfig, delta_entries::AbstractVector{<:NamedTuple})
    refs = ("B_perp",)
    ref_labels = Dict(
        "B_perp" => LaTeXString("B_{\\perp}"),
    )
    ref_colors = Dict("B_perp" => :dodgerblue3)

    tag_to_Leff = Dict{String, Float64}()
    deltas_lookup = Dict{Tuple{String, String}, Vector{Float64}}()
    for d in delta_entries
        tag_to_Leff[d.map_tag] = d.llarge_eff_pc
        deltas_lookup[(d.map_tag, d.reference)] = d.deltas
    end
    tags = collect(keys(tag_to_Leff))
    sort!(tags; by=t -> t == "nofilter" ? -Inf : tag_to_Leff[t])
    cols = _cycled_curve_colors(length(tags))
    tag_colors = Dict(tags[i] => cols[i] for i in eachindex(tags))

    with_theme(theme_latexfonts()) do
        set_theme_spectra!()

        ngrid = 721
        xkde, _ = _kde_curve_density(Float64[]; xmin=-90.0, xmax=90.0, ngrid=ngrid)
        zero_kde = zeros(Float64, length(xkde))
        kde_lookup = Dict{Tuple{String, String}, Vector{Float64}}()
        ymax_shared = 0.0
        for tag in tags, ref in refs
            vals = get(deltas_lookup, (tag, ref), Float64[])
            _, yk = _kde_curve_density(vals; xmin=-90.0, xmax=90.0, ngrid=ngrid)
            kde_lookup[(tag, ref)] = yk
            ymax_shared = max(ymax_shared, maximum(yk))
        end
        ymax_shared = ymax_shared > 0 ? 1.15 * ymax_shared : 1.0

        if cfg.channel_alignment_pdf_plain
            fig = Figure(size=(1400, 800), figure_padding=30)
            for (k, ref) in enumerate(refs)
                ax = Axis(fig[k, 1],
                    xlabel=LaTeXString("\\Delta\\theta\\,[^\\circ]"),
                    ylabel=LaTeXString("p(\\Delta\\theta)"),
                    title=ref_labels[ref],
                )
                ax.xminorticksvisible = false
                ax.yminorticksvisible = false
                ax.xticks = -90:45:90
                if k != length(refs)
                    ax.xlabelvisible = false
                    ax.xticklabelsvisible = false
                end

                ncurves = 0
                for tag in tags
                    vals = get(deltas_lookup, (tag, ref), Float64[])
                    isempty(vals) && continue
                    yk = get(kde_lookup, (tag, ref), zero_kde)
                    band!(ax, xkde, zero_kde, yk; color=(tag_colors[tag], 0.18))
                    lines!(ax, xkde, yk; color=tag_colors[tag], linewidth=4,
                        label=_filter_display_label_latex(tag, tag_to_Leff[tag]))
                    ncurves += 1
                end
                ncurves > 0 && axislegend(ax; position=:rt, framevisible=true, labelsize=34)

                # Visual guides: parallel band around 0 deg and perpendicular bands near +/-90 deg.
                vspan!(ax, -_ALIGNMENT_ANGLE_TOL_DEG, _ALIGNMENT_ANGLE_TOL_DEG; color=(:gray60, 0.16))
                vspan!(ax, 90.0 - _ALIGNMENT_ANGLE_TOL_DEG, 90.0; color=(:seagreen4, 0.12))
                vspan!(ax, -90.0, -90.0 + _ALIGNMENT_ANGLE_TOL_DEG; color=(:seagreen4, 0.12))

                xlims!(ax, -90.0, 90.0)
                ylims!(ax, 0.0, ymax_shared)
                ax.xtickformat = (vs) -> latex_linear_tickformat(vs; digits=0)
                ax.ytickformat = (vs) -> latex_linear_tickformat(vs; digits=3)
            end

            out = _save_figure(cfg, fig, "channel_alignment_delta_theta_hist.pdf")
            display(fig)
            return out
        else
            fig = Figure(size=(1400, 1000), figure_padding=30)
            selected_idx = Observable(1)
            kde_obs = Dict(ref => Observable(copy(get(kde_lookup, (tags[1], ref), zero_kde))) for ref in refs)

            for ref in refs
                kde_obs[ref][] = copy(get(kde_lookup, (tags[1], ref), zero_kde))
            end

            for (k, ref) in enumerate(refs)
                ax = Axis(fig[k, 1],
                    xlabel=LaTeXString("\\Delta\\theta\\,[^\\circ]"),
                    ylabel=LaTeXString("p(\\Delta\\theta)"),
                    title=ref_labels[ref],
                )
                ax.xminorticksvisible = false
                ax.yminorticksvisible = false
                ax.xticks = -90:45:90
                if k != length(refs)
                    ax.xlabelvisible = false
                    ax.xticklabelsvisible = false
                end

                band!(ax, xkde, zero_kde, kde_obs[ref]; color=(ref_colors[ref], 0.18))
                lines!(ax, xkde, kde_obs[ref]; color=ref_colors[ref], linewidth=4)

                # Visual guides: parallel band around 0 deg and perpendicular bands near +/-90 deg.
                vspan!(ax, -_ALIGNMENT_ANGLE_TOL_DEG, _ALIGNMENT_ANGLE_TOL_DEG; color=(:gray60, 0.16))
                vspan!(ax, 90.0 - _ALIGNMENT_ANGLE_TOL_DEG, 90.0; color=(:seagreen4, 0.12))
                vspan!(ax, -90.0, -90.0 + _ALIGNMENT_ANGLE_TOL_DEG; color=(:seagreen4, 0.12))

                xlims!(ax, -90.0, 90.0)
                ylims!(ax, 0.0, ymax_shared)
                ax.xtickformat = (vs) -> latex_linear_tickformat(vs; digits=0)
                ax.ytickformat = (vs) -> latex_linear_tickformat(vs; digits=3)
            end

            slider = Slider(fig[2, 1];
                range=1:length(tags),
                startvalue=1,
                horizontal=true,
            )
            filter_label = lift(selected_idx) do idx
                tag = tags[idx]
                return @sprintf("Filter: %s", _filter_display_label(tag, tag_to_Leff[tag]))
            end
            Label(fig[3, 1], filter_label; fontsize=28)

            on(slider.value) do idx
                selected_idx[] = idx
                for ref in refs
                    kde_obs[ref][] = copy(get(kde_lookup, (tags[idx], ref), zero_kde))
                end
            end

            out = _save_figure(cfg, fig, "channel_alignment_delta_theta_hist.pdf")
            display(fig)
            return out
        end
    end
end

"""
    _run_channel_b_alignment(...)

Measures `Δθ = θ_canal - θ_{B_perp}` over valid image pixels.
"""
function _run_channel_b_alignment(cfg::InstrumentalConfig, data::FilterPassResult)
    b1, b2, bperp = _project_bperp_maps(cfg)
    θB = _orientation_map_bperp_from_b1_b2(b1, b2)

    maps = _ALIGNMENT_MAP_ENTRY[(map_tag="nofilter", llarge_eff_pc=0.0, Pmap=Float64.(data.Pmax0))]
    for Llarge in sort(data.L_ok)
        Leff = (cfg.Lbox_pc / cfg.n) * Llarge
        push!(maps, (
            map_tag=_filter_tag(Llarge),
            llarge_eff_pc=Leff,
            Pmap=Float64.(data.Pmax_filt[Llarge]),
        ))
    end

    rows = _ALIGNMENT_ROW[]
    delta_entries = _ALIGNMENT_DELTA_ENTRY[]
    for entry in maps
        θcanal, canal_mask, _, _ = _canal_orientation_map(entry.Pmap, data.scales.Δx, data.scales.Δy)
        mask_B = canal_mask .& isfinite.(θB) .& isfinite.(bperp) .& (bperp .> 0)
        deltas = _alignment_deltas_deg(θcanal, θB, mask_B)
        push!(delta_entries, (
            map_tag=entry.map_tag,
            llarge_eff_pc=entry.llarge_eff_pc,
            reference="B_perp",
            deltas=deltas,
        ))
        stats = _alignment_stats_from_deltas(deltas; tol_deg=_ALIGNMENT_ANGLE_TOL_DEG)
        push!(rows, (
            map_tag=entry.map_tag,
            llarge_eff_pc=entry.llarge_eff_pc,
            reference="B_perp",
            npix=stats.npix,
            mean_abs_delta_deg=stats.mean_abs_delta_deg,
            median_abs_delta_deg=stats.median_abs_delta_deg,
            frac_parallel_15deg=stats.frac_parallel_tol,
            frac_perp_15deg=stats.frac_perp_tol,
            median_perp_offset_deg=stats.median_perp_offset_deg,
        ))
    end

    summary_path = _write_channel_alignment_csv(joinpath(cfg.base_out, "channel_alignment_summary.csv"), rows)
    figure_path = _plot_delta_theta_histograms(cfg, delta_entries)
    @info "Channel-angle alignment analysis complete" summary_path=summary_path figure_path=figure_path los=cfg.los

    return (summary_path=summary_path, figure_path=figure_path, rows=rows)
end
