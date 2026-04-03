const _CANAL_MORPH_NEIGHBORS = (
    CartesianIndex(-1, -1), CartesianIndex(-1, 0), CartesianIndex(-1, 1),
    CartesianIndex(0, -1),                         CartesianIndex(0, 1),
    CartesianIndex(1, -1),  CartesianIndex(1, 0), CartesianIndex(1, 1),
)
const _CANAL_MORPH_MIN_BRANCH_PIX = 4
const _CANAL_MORPH_CURVATURE_SMOOTH_RADIUS = 1
const _CANAL_MORPH_PROFILE_HALF_RANGE_PIX = 12.0
const _CANAL_MORPH_PROFILE_STEP_PIX = 0.25
const _CANAL_MORPH_RELATIVE_WIDTH_LEVEL = 0.2
const _CANAL_MORPH_RESOLVED_WIDTH_PIX = 2.5
const _CANAL_MORPH_LOG_BINS = 18
const _CANAL_MORPH_AXIS_LABELSIZE = 36
const _CANAL_MORPH_AXIS_TICKSIZE = 28
const _CANAL_MORPH_LEGEND_LABELSIZE = 24

@inline _branch_sort_key(idx::CartesianIndex{2}) = Tuple(idx)
@inline _edge_key(a::CartesianIndex{2}, b::CartesianIndex{2}) = isless(_branch_sort_key(a), _branch_sort_key(b)) ? (a, b) : (b, a)

@inline function _branch_step_length_pc(a::CartesianIndex{2}, b::CartesianIndex{2}, Δx::Real, Δy::Real)
    di = (b[1] - a[1]) * float(Δx)
    dj = (b[2] - a[2]) * float(Δy)
    return hypot(di, dj)
end

function _skeleton_component_graph(comp::Vector{CartesianIndex{2}})
    comp_set = Set(comp)
    neighbors = Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}}()
    for idx in comp
        nbs = CartesianIndex{2}[]
        for δ in _CANAL_MORPH_NEIGHBORS
            nb = idx + δ
            nb in comp_set && push!(nbs, nb)
        end
        sort!(nbs; by=_branch_sort_key)
        neighbors[idx] = nbs
    end
    return neighbors
end

function _trace_skeleton_chain(start::CartesianIndex{2}, next::CartesianIndex{2},
                               neighbors::Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}},
                               degrees::Dict{CartesianIndex{2}, Int},
                               used::Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}})
    seq = CartesianIndex{2}[start]
    prev = start
    curr = next
    push!(used, _edge_key(start, next))

    while true
        push!(seq, curr)
        get(degrees, curr, 0) == 2 || break
        nbs = get(neighbors, curr, CartesianIndex{2}[])
        nexts = [nb for nb in nbs if nb != prev]
        length(nexts) == 1 || break
        nxt = only(nexts)
        ekey = _edge_key(curr, nxt)
        ekey in used && break
        push!(used, ekey)
        prev, curr = curr, nxt
    end

    return seq
end

function _trace_skeleton_loop(start::CartesianIndex{2},
                              neighbors::Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}},
                              used::Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}})
    seq = CartesianIndex{2}[start]
    nbs0 = get(neighbors, start, CartesianIndex{2}[])
    isempty(nbs0) && return seq

    prev = start
    curr = nbs0[1]
    push!(used, _edge_key(start, curr))

    while true
        push!(seq, curr)
        nbs = get(neighbors, curr, CartesianIndex{2}[])
        nexts = [nb for nb in nbs if nb != prev]
        isempty(nexts) && break
        nxt = nexts[1]
        nxt == start && break
        ekey = _edge_key(curr, nxt)
        ekey in used && break
        push!(used, ekey)
        prev, curr = curr, nxt
    end

    return seq
end

"""
    _extract_skeleton_branches(...)

Splits a one-pixel skeleton into ordered branches. Junctions are cut at pixels
with degree different from two, so geodesic length and curvature remain
well-defined on each returned branch.
"""
function _extract_skeleton_branches(mask::AbstractMatrix{Bool})
    comps = _zone_connected_components(mask)
    branches = Vector{NamedTuple}()

    for (component_id, comp) in enumerate(comps)
        isempty(comp) && continue
        neighbors = _skeleton_component_graph(comp)
        degrees = Dict(idx => length(neighbors[idx]) for idx in keys(neighbors))
        used = Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}}()
        local_branch_id = 0

        nodes = sort(collect(keys(neighbors)); by=_branch_sort_key)
        branch_nodes = [idx for idx in nodes if degrees[idx] != 2]

        if isempty(branch_nodes)
            seq = _trace_skeleton_loop(first(nodes), neighbors, used)
            if length(seq) >= _CANAL_MORPH_MIN_BRANCH_PIX
                local_branch_id += 1
                push!(branches, (
                    component_id=component_id,
                    branch_id_local=local_branch_id,
                    pixels=seq,
                    is_loop=true,
                ))
            end
            continue
        end

        for node in branch_nodes
            for nb in neighbors[node]
                ekey = _edge_key(node, nb)
                ekey in used && continue
                seq = _trace_skeleton_chain(node, nb, neighbors, degrees, used)
                length(seq) >= _CANAL_MORPH_MIN_BRANCH_PIX || continue
                local_branch_id += 1
                push!(branches, (
                    component_id=component_id,
                    branch_id_local=local_branch_id,
                    pixels=seq,
                    is_loop=false,
                ))
            end
        end
    end

    return branches
end

function _fallback_branches_from_mask(mask::AbstractMatrix{Bool}, score::AbstractMatrix, θcanal::AbstractMatrix)
    comps = _zone_connected_components(mask)
    branches = Vector{NamedTuple}()

    for (component_id, comp) in enumerate(comps)
        isempty(comp) && continue
        θvals = Float64[θcanal[I] for I in comp if isfinite(θcanal[I])]
        isempty(θvals) && continue

        # Line orientations are π-periodic, so average doubled angles.
        c2 = mean(cos.(2 .* θvals))
        s2 = mean(sin.(2 .* θvals))
        θdom = _wrap_orientation_pi(0.5 * atan(s2, c2))
        group_by_i = abs(cos(θdom)) >= abs(sin(θdom))

        best = Dict{Int, CartesianIndex{2}}()
        best_score = Dict{Int, Float64}()
        for idx in comp
            key = group_by_i ? idx[1] : idx[2]
            s = float(score[idx])
            if !haskey(best, key) || s > best_score[key]
                best[key] = idx
                best_score[key] = s
            end
        end

        keys_sorted = sort(collect(keys(best)))
        seq = CartesianIndex{2}[best[k] for k in keys_sorted]
        length(seq) >= _CANAL_MORPH_MIN_BRANCH_PIX || continue
        push!(branches, (
            component_id=component_id,
            branch_id_local=1,
            pixels=seq,
            is_loop=false,
        ))
    end

    return branches
end

function _branch_arclengths(pixels::Vector{CartesianIndex{2}}, Δx::Real, Δy::Real)
    n = length(pixels)
    s = zeros(Float64, n)
    for k in 2:n
        s[k] = s[k - 1] + _branch_step_length_pc(pixels[k - 1], pixels[k], Δx, Δy)
    end
    return s
end

function _unwrap_line_angles(θvals::AbstractVector{<:Real})
    n = length(θvals)
    n == 0 && return Float64[]
    out = zeros(Float64, n)
    out[1] = float(θvals[1])
    for i in 2:n
        θ = float(θvals[i])
        while θ - out[i - 1] > π / 2
            θ -= π
        end
        while θ - out[i - 1] < -π / 2
            θ += π
        end
        out[i] = θ
    end
    return out
end

function _smooth_vector(vals::AbstractVector{<:Real}; radius::Int=_CANAL_MORPH_CURVATURE_SMOOTH_RADIUS)
    n = length(vals)
    out = fill(NaN, n)
    for i in 1:n
        lo = max(1, i - radius)
        hi = min(n, i + radius)
        acc = 0.0
        count = 0
        for j in lo:hi
            v = vals[j]
            if isfinite(v)
                acc += float(v)
                count += 1
            end
        end
        count > 0 && (out[i] = acc / count)
    end
    return out
end

function _curvature_from_branch(θvals::AbstractVector{<:Real}, s::AbstractVector{<:Real})
    n = length(θvals)
    n == length(s) || error("Angle/arc-length size mismatch: $(n) vs $(length(s))")
    n >= 3 || return Float64[]

    θu = _unwrap_line_angles(θvals)
    θs = _smooth_vector(θu)
    κ = fill(NaN, n)

    if isfinite(θs[1]) && isfinite(θs[2]) && s[2] > s[1]
        κ[1] = abs((θs[2] - θs[1]) / (s[2] - s[1]))
    end
    for i in 2:(n - 1)
        ds = s[i + 1] - s[i - 1]
        if ds > 0 && isfinite(θs[i - 1]) && isfinite(θs[i + 1])
            κ[i] = abs((θs[i + 1] - θs[i - 1]) / ds)
        end
    end
    if isfinite(θs[n]) && isfinite(θs[n - 1]) && s[n] > s[n - 1]
        κ[n] = abs((θs[n] - θs[n - 1]) / (s[n] - s[n - 1]))
    end

    return κ[isfinite.(κ)]
end

function _score_interpolator(score::AbstractMatrix)
    return interpolate(Float64.(score), BSpline(Linear()))
end

@inline function _sample_itp(itp, x::Real, y::Real, nx::Int, ny::Int)
    (1 <= x <= nx && 1 <= y <= ny) || return NaN
    return itp(x, y)
end

function _profile_side_crossing(us::AbstractVector{<:Real}, vals::AbstractVector{<:Real}, level::Real)
    length(us) == length(vals) || error("Profile size mismatch")
    length(us) >= 2 || return NaN
    for k in 2:length(us)
        v1 = vals[k - 1]
        v2 = vals[k]
        isfinite(v1) && isfinite(v2) || return NaN
        if (v1 >= level && v2 <= level) || (v1 <= level && v2 >= level)
            dv = v2 - v1
            abs(dv) <= eps(Float64) && return float(us[k])
            t = clamp((float(level) - float(v1)) / dv, 0.0, 1.0)
            return float(us[k - 1]) + t * (float(us[k]) - float(us[k - 1]))
        end
    end
    return NaN
end

function _normal_width_metrics(score_itp, score::AbstractMatrix, i::Int, j::Int, θ::Real, Δx::Real, Δy::Real;
                               half_range_pix::Real=_CANAL_MORPH_PROFILE_HALF_RANGE_PIX,
                               step_pix::Real=_CANAL_MORPH_PROFILE_STEP_PIX,
                               rel_level::Real=_CANAL_MORPH_RELATIVE_WIDTH_LEVEL)
    center = float(score[i, j])
    isfinite(center) && center > 0 || return (width_fwhm_pc=NaN, width_rel_pc=NaN)

    nx, ny = size(score)
    half_range_pc = half_range_pix * min(float(Δx), float(Δy))
    step_pc = step_pix * min(float(Δx), float(Δy))
    us = collect(0.0:step_pc:half_range_pc)
    isempty(us) && return (width_fwhm_pc=NaN, width_rel_pc=NaN)

    nvecx = -sin(float(θ))
    nvecy = cos(float(θ))
    vals_pos = similar(us)
    vals_neg = similar(us)

    for k in eachindex(us)
        u = us[k]
        di = (u * nvecx) / float(Δx)
        dj = (u * nvecy) / float(Δy)
        vals_pos[k] = _sample_itp(score_itp, i + di, j + dj, nx, ny)
        vals_neg[k] = _sample_itp(score_itp, i - di, j - dj, nx, ny)
    end

    w_half_pos = _profile_side_crossing(us, vals_pos, 0.5 * center)
    w_half_neg = _profile_side_crossing(us, vals_neg, 0.5 * center)
    w_rel_pos = _profile_side_crossing(us, vals_pos, rel_level * center)
    w_rel_neg = _profile_side_crossing(us, vals_neg, rel_level * center)

    width_fwhm_pc = (isfinite(w_half_pos) && isfinite(w_half_neg)) ? (w_half_pos + w_half_neg) : NaN
    width_rel_pc = (isfinite(w_rel_pos) && isfinite(w_rel_neg)) ? (w_rel_pos + w_rel_neg) : NaN
    return (width_fwhm_pc=width_fwhm_pc, width_rel_pc=width_rel_pc)
end

function _branch_morphology_row(entry, branch, branch_id::Int, θcanal::AbstractMatrix, score::AbstractMatrix, score_itp,
                                Δx::Real, Δy::Real)
    pixels = branch.pixels
    s = _branch_arclengths(pixels, Δx, Δy)
    length_pc = isempty(s) ? NaN : s[end] + (branch.is_loop && length(pixels) > 1 ? _branch_step_length_pc(pixels[end], pixels[1], Δx, Δy) : 0.0)

    widths_fwhm = Float64[]
    widths_rel = Float64[]
    θvals = Float64[]
    svals = Float64[]
    score_vals = Float64[]

    for (k, idx) in enumerate(pixels)
        i, j = Tuple(idx)
        θ = θcanal[i, j]
        isfinite(θ) || continue
        widths = _normal_width_metrics(score_itp, score, i, j, θ, Δx, Δy)
        isfinite(widths.width_fwhm_pc) && push!(widths_fwhm, widths.width_fwhm_pc)
        isfinite(widths.width_rel_pc) && push!(widths_rel, widths.width_rel_pc)
        push!(θvals, float(θ))
        push!(svals, s[k])
        push!(score_vals, float(score[i, j]))
    end

    curvature_vals = _curvature_from_branch(θvals, svals)
    return (
        map_tag=entry.map_tag,
        llarge_eff_pc=float(entry.llarge_eff_pc),
        component_id=branch.component_id,
        branch_id=branch_id,
        branch_id_local=branch.branch_id_local,
        is_loop=branch.is_loop,
        npix=length(pixels),
        length_pc=length_pc,
        width_fwhm_median_pc=_finite_quantile_or_nan(widths_fwhm, 0.50),
        width_fwhm_q16_pc=_finite_quantile_or_nan(widths_fwhm, 0.16),
        width_fwhm_q84_pc=_finite_quantile_or_nan(widths_fwhm, 0.84),
        width_rel_median_pc=_finite_quantile_or_nan(widths_rel, 0.50),
        width_rel_q16_pc=_finite_quantile_or_nan(widths_rel, 0.16),
        width_rel_q84_pc=_finite_quantile_or_nan(widths_rel, 0.84),
        width_n_samples=length(widths_fwhm),
        curvature_median_pc_inv=_finite_quantile_or_nan(curvature_vals, 0.50),
        curvature_p90_pc_inv=_finite_quantile_or_nan(curvature_vals, 0.90),
        curvature_max_pc_inv=isempty(curvature_vals) ? NaN : maximum(curvature_vals),
        mean_score=isempty(score_vals) ? NaN : mean(score_vals),
    )
end

function _analyze_canal_morphology_map(Pmap::AbstractMatrix, Δx::Real, Δy::Real; map_tag::AbstractString="nofilter", llarge_eff_pc::Real=0.0)
    θcanal, mask, score, _ = _canal_orientation_map(Pmap, Δx, Δy)
    branches = _extract_skeleton_branches(mask)
    isempty(branches) && count(mask) > 0 && (branches = _fallback_branches_from_mask(mask, score, θcanal))
    score_itp = _score_interpolator(score)
    rows = Vector{NamedTuple}()
    branch_id = 0

    for branch in branches
        branch_id += 1
        push!(rows, _branch_morphology_row((map_tag=map_tag, llarge_eff_pc=llarge_eff_pc), branch, branch_id, θcanal, score, score_itp, Δx, Δy))
    end

    return (rows=rows, θcanal=θcanal, mask=mask, score=score)
end

function _positive_finite(vals::AbstractVector{<:Real}; xmin::Real=0.0)
    return Float64[v for v in vals if isfinite(v) && v > max(float(xmin), 0.0)]
end

function _power_law_alpha(vals::AbstractVector{<:Real}; xmin::Real)
    xv = _positive_finite(vals; xmin=xmin)
    n = length(xv)
    n >= 2 || return (alpha=NaN, n=n)
    denom = sum(log.(xv ./ float(xmin)))
    denom > 0 || return (alpha=NaN, n=n)
    return (alpha=1.0 + n / denom, n=n)
end

function _loglog_relation(widths::AbstractVector{<:Real}, lengths::AbstractVector{<:Real}; xmin::Real, ymin::Real)
    good = isfinite.(widths) .& isfinite.(lengths) .& (widths .>= xmin) .& (lengths .>= ymin) .& (widths .> 0) .& (lengths .> 0)
    x = log10.(Float64.(widths[good]))
    y = log10.(Float64.(lengths[good]))
    length(x) >= 2 || return (slope=NaN, intercept=NaN, corr=NaN, n=length(x))

    xmean = mean(x)
    ymean = mean(y)
    varx = sum((x .- xmean) .^ 2)
    varx > 0 || return (slope=NaN, intercept=NaN, corr=NaN, n=length(x))

    covxy = sum((x .- xmean) .* (y .- ymean))
    slope = covxy / varx
    intercept = ymean - slope * xmean
    corr = cor(x, y)
    return (slope=slope, intercept=intercept, corr=corr, n=length(x))
end

function _write_canal_morphology_branch_csv(path::AbstractString, rows::AbstractVector{<:NamedTuple})
    open(path, "w") do io
        println(io, "map_tag,llarge_eff_pc,component_id,branch_id,branch_id_local,is_loop,npix,length_pc,width_fwhm_median_pc,width_fwhm_q16_pc,width_fwhm_q84_pc,width_rel_median_pc,width_rel_q16_pc,width_rel_q84_pc,width_n_samples,curvature_median_pc_inv,curvature_p90_pc_inv,curvature_max_pc_inv,mean_score")
        for row in rows
            println(io, string(
                row.map_tag, ",", row.llarge_eff_pc, ",", row.component_id, ",", row.branch_id, ",",
                row.branch_id_local, ",", row.is_loop, ",", row.npix, ",", row.length_pc, ",",
                row.width_fwhm_median_pc, ",", row.width_fwhm_q16_pc, ",", row.width_fwhm_q84_pc, ",",
                row.width_rel_median_pc, ",", row.width_rel_q16_pc, ",", row.width_rel_q84_pc, ",",
                row.width_n_samples, ",", row.curvature_median_pc_inv, ",", row.curvature_p90_pc_inv, ",",
                row.curvature_max_pc_inv, ",", row.mean_score
            ))
        end
    end
    return path
end

function _write_canal_morphology_summary_csv(path::AbstractString, rows::AbstractVector{<:NamedTuple})
    open(path, "w") do io
        println(io, "map_tag,llarge_eff_pc,n_branches,n_resolved_width,n_resolved_length,width_xmin_pc,length_xmin_pc,alpha_width,alpha_length,loglog_slope_L_of_w,loglog_intercept_L_of_w,loglog_corr")
        for row in rows
            println(io, string(
                row.map_tag, ",", row.llarge_eff_pc, ",", row.n_branches, ",", row.n_resolved_width, ",",
                row.n_resolved_length, ",", row.width_xmin_pc, ",", row.length_xmin_pc, ",",
                row.alpha_width, ",", row.alpha_length, ",", row.loglog_slope_L_of_w, ",",
                row.loglog_intercept_L_of_w, ",", row.loglog_corr
            ))
        end
    end
    return path
end

function _histogram_density(vals::AbstractVector{<:Real}, edges::AbstractVector{<:Real}; geometric_centers::Bool=false)
    counts = zeros(Float64, length(edges) - 1)
    for v in vals
        isfinite(v) && v > 0 || continue
        idx = searchsortedlast(edges, float(v)) - 1
        idx < 1 && continue
        idx > length(counts) && (idx = length(counts))
        counts[idx] += 1
    end
    widths = diff(Float64.(edges))
    total = sum(counts)
    density = total > 0 ? (counts ./ total) ./ widths : fill(0.0, length(counts))
    centers = geometric_centers ?
        sqrt.(Float64.(edges[1:end-1]) .* Float64.(edges[2:end])) :
        0.5 .* (Float64.(edges[1:end-1]) .+ Float64.(edges[2:end]))
    return centers, density
end

function _metric_groups(rows::AbstractVector{<:NamedTuple}, metric::Symbol)
    groups = Dict{String, Vector{Float64}}()
    leffs = Dict{String, Float64}()
    for row in rows
        haskey(groups, row.map_tag) || (groups[row.map_tag] = Float64[])
        isfinite(row.llarge_eff_pc) && !haskey(leffs, row.map_tag) && (leffs[row.map_tag] = row.llarge_eff_pc)
        value = getproperty(row, metric)
        isfinite(value) && value > 0 && push!(groups[row.map_tag], float(value))
    end
    return groups, leffs
end

function _plot_morphology_distribution(cfg::InstrumentalConfig, rows::AbstractVector{<:NamedTuple}, metric::Symbol;
                                       xlabel::AbstractString, ylabel::AbstractString, filename::AbstractString)
    groups, leffs = _metric_groups(rows, metric)
    allvals = reduce(vcat, values(groups); init=Float64[])
    allvals = _positive_finite(allvals)
    isempty(allvals) && return nothing

    lo = minimum(allvals)
    hi = maximum(allvals)
    if !(lo < hi)
        lo = max(0.8 * lo, eps(Float64))
        hi = max(1.2 * hi, 1.25 * lo)
    end
    edges = exp.(range(log(lo), log(hi); length=_CANAL_MORPH_LOG_BINS + 1))
    tags = sort(collect(keys(groups)); by=tag -> tag == "nofilter" ? (-Inf) : get(leffs, tag, Inf))
    cols = _cycled_curve_colors(length(tags))

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1200, 860))
        ax = Axis(
            fig[1, 1],
            xlabel=LaTeXString(xlabel),
            ylabel=LaTeXString(ylabel),
            xgridvisible=false,
            ygridvisible=false,
            xlabelsize=_CANAL_MORPH_AXIS_LABELSIZE,
            ylabelsize=_CANAL_MORPH_AXIS_LABELSIZE,
            xticklabelsize=_CANAL_MORPH_AXIS_TICKSIZE,
            yticklabelsize=_CANAL_MORPH_AXIS_TICKSIZE,
        )
        yvals = Float64[]

        for (idx, tag) in enumerate(tags)
            vals = groups[tag]
            isempty(vals) && continue
            centers, density = _histogram_density(vals, edges; geometric_centers=true)
            ok = density .> 0
            any(ok) || continue
            append!(yvals, density[ok])
            lines!(ax, centers[ok], density[ok]; color=cols[idx], linewidth=3,
                   label=_filter_display_label_latex(tag, get(leffs, tag, 0.0)))
            scatter!(ax, centers[ok], density[ok]; color=cols[idx], markersize=12)
        end

        isempty(yvals) && return nothing
        configure_axis_style!(ax; xdata=allvals, ydata=yvals, xscale=:log10, yscale=:log10, xminor_subdiv=0, yminor_subdiv=0)
        axislegend(ax; position=:rt, labelsize=_CANAL_MORPH_LEGEND_LABELSIZE)
        return _save_figure(cfg, fig, filename)
    end
end

function _plot_length_width_relation(cfg::InstrumentalConfig, rows::AbstractVector{<:NamedTuple}, summary_rows::AbstractVector{<:NamedTuple})
    tags = sort(unique(row.map_tag for row in rows); by=tag -> tag == "nofilter" ? (-Inf) : minimum(r.llarge_eff_pc for r in rows if r.map_tag == tag))
    isempty(tags) && return nothing
    cols = _cycled_curve_colors(length(tags))

    widths_all = Float64[row.width_fwhm_median_pc for row in rows if isfinite(row.width_fwhm_median_pc) && row.width_fwhm_median_pc > 0]
    lengths_all = Float64[row.length_pc for row in rows if isfinite(row.length_pc) && row.length_pc > 0]
    isempty(widths_all) && return nothing
    isempty(lengths_all) && return nothing

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1200, 860))
        ax = Axis(
            fig[1, 1],
            xlabel=LaTeXString("w_{\\mathrm{FWHM}}\\,[\\mathrm{pc}]"),
            ylabel=LaTeXString("L\\,[\\mathrm{pc}]"),
            xgridvisible=false,
            ygridvisible=false,
            xlabelsize=_CANAL_MORPH_AXIS_LABELSIZE,
            ylabelsize=_CANAL_MORPH_AXIS_LABELSIZE,
            xticklabelsize=_CANAL_MORPH_AXIS_TICKSIZE,
            yticklabelsize=_CANAL_MORPH_AXIS_TICKSIZE,
        )

        for (idx, tag) in enumerate(tags)
            vals = [row for row in rows if row.map_tag == tag && isfinite(row.width_fwhm_median_pc) && row.width_fwhm_median_pc > 0 && isfinite(row.length_pc) && row.length_pc > 0]
            isempty(vals) && continue
            ws = Float64[row.width_fwhm_median_pc for row in vals]
            ls = Float64[row.length_pc for row in vals]
            scatter!(ax, ws, ls; color=(cols[idx], 0.55), markersize=11,
                     label=_filter_display_label_latex(tag, minimum(row.llarge_eff_pc for row in vals)))

            summary = findfirst(r -> r.map_tag == tag, summary_rows)
            summary === nothing && continue
            fit = summary_rows[summary]
            if isfinite(fit.loglog_slope_L_of_w) && isfinite(fit.loglog_intercept_L_of_w)
                xmin = max(fit.width_xmin_pc, minimum(ws))
                xmax = maximum(ws)
                if xmin < xmax
                    xfit = exp.(range(log(xmin), log(xmax); length=100))
                    yfit = 10.0 .^ (fit.loglog_intercept_L_of_w .+ fit.loglog_slope_L_of_w .* log10.(xfit))
                    lines!(ax, xfit, yfit; color=cols[idx], linewidth=3, linestyle=:dash)
                end
            end
        end

        configure_axis_style!(ax; xdata=widths_all, ydata=lengths_all, xscale=:log10, yscale=:log10, xminor_subdiv=0, yminor_subdiv=0)
        axislegend(ax; position=:rt, labelsize=_CANAL_MORPH_LEGEND_LABELSIZE)
        return _save_figure(cfg, fig, "canal_length_width_relation.pdf")
    end
end

function _build_canal_morphology_summary(rows::AbstractVector{<:NamedTuple}, Δx::Real, Δy::Real)
    tags = sort(unique(row.map_tag for row in rows); by=tag -> tag == "nofilter" ? (-Inf) : minimum(r.llarge_eff_pc for r in rows if r.map_tag == tag))
    width_xmin = _CANAL_MORPH_RESOLVED_WIDTH_PIX * min(float(Δx), float(Δy))
    length_xmin = _CANAL_MORPH_MIN_BRANCH_PIX * min(float(Δx), float(Δy))
    summary_rows = Vector{NamedTuple}()

    for tag in tags
        sub = [row for row in rows if row.map_tag == tag]
        widths = Float64[row.width_fwhm_median_pc for row in sub]
        lengths = Float64[row.length_pc for row in sub]
        alpha_w = _power_law_alpha(widths; xmin=width_xmin)
        alpha_l = _power_law_alpha(lengths; xmin=length_xmin)
        fit = _loglog_relation(widths, lengths; xmin=width_xmin, ymin=length_xmin)
        push!(summary_rows, (
            map_tag=tag,
            llarge_eff_pc=minimum(Float64[row.llarge_eff_pc for row in sub]),
            n_branches=length(sub),
            n_resolved_width=alpha_w.n,
            n_resolved_length=alpha_l.n,
            width_xmin_pc=width_xmin,
            length_xmin_pc=length_xmin,
            alpha_width=alpha_w.alpha,
            alpha_length=alpha_l.alpha,
            loglog_slope_L_of_w=fit.slope,
            loglog_intercept_L_of_w=fit.intercept,
            loglog_corr=fit.corr,
        ))
    end

    return summary_rows
end

function _run_canal_morphology(cfg::InstrumentalConfig, data::FilterPassResult)
    maps = [(map_tag="nofilter", llarge_eff_pc=0.0, Pmap=Float64.(data.Pmax0))]
    for Llarge in sort(data.L_ok)
        Leff = (cfg.Lbox_pc / cfg.n) * Llarge
        push!(maps, (
            map_tag=_filter_tag(Llarge),
            llarge_eff_pc=Leff,
            Pmap=Float64.(data.Pmax_filt[Llarge]),
        ))
    end

    rows = Vector{NamedTuple}()
    for entry in maps
        result = _analyze_canal_morphology_map(entry.Pmap, data.scales.Δx, data.scales.Δy;
            map_tag=entry.map_tag, llarge_eff_pc=entry.llarge_eff_pc)
        append!(rows, result.rows)
    end

    branch_csv = joinpath(cfg.base_out, "canal_morphology_branches.csv")
    summary_rows = _build_canal_morphology_summary(rows, data.scales.Δx, data.scales.Δy)
    summary_csv = joinpath(cfg.base_out, "canal_morphology_scaling_summary.csv")
    _write_canal_morphology_branch_csv(branch_csv, rows)
    _write_canal_morphology_summary_csv(summary_csv, summary_rows)

    width_fig = _plot_morphology_distribution(cfg, rows, :width_fwhm_median_pc;
        xlabel="w_{\\mathrm{FWHM}}\\,[\\mathrm{pc}]",
        ylabel="\\mathrm{d}N/\\mathrm{d}w",
        filename="canal_width_distribution.pdf")
    length_fig = _plot_morphology_distribution(cfg, rows, :length_pc;
        xlabel="L\\,[\\mathrm{pc}]",
        ylabel="\\mathrm{d}N/\\mathrm{d}L",
        filename="canal_length_distribution.pdf")
    curvature_fig = _plot_morphology_distribution(cfg, rows, :curvature_median_pc_inv;
        xlabel="\\kappa\\,[\\mathrm{pc}^{-1}]",
        ylabel="\\mathrm{d}N/\\mathrm{d}\\kappa",
        filename="canal_curvature_distribution.pdf")
    relation_fig = _plot_length_width_relation(cfg, rows, summary_rows)

    @info "Canal morphology analysis complete" branch_csv=branch_csv summary_csv=summary_csv
    return (
        branch_csv=branch_csv,
        summary_csv=summary_csv,
        width_figure=width_fig,
        length_figure=length_fig,
        curvature_figure=curvature_fig,
        relation_figure=relation_fig,
        rows=rows,
        summary_rows=summary_rows,
    )
end
