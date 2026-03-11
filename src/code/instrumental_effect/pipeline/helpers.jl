"""
    _filter_tag(...)

    Builds filter tag text.
"""
function _filter_tag(Llarge::Real)
    @sprintf("HardBandPass_remove_L0_to_1pc_and_%dto50pc", Int(round(Llarge)))
end

@inline function _filter_display_label(map_tag::AbstractString, llarge_eff_pc::Real)
    map_tag == "nofilter" && return "No filter"
    return @sprintf("L_eff = %.2f pc", float(llarge_eff_pc))
end

@inline function _filter_display_label_latex(map_tag::AbstractString, llarge_eff_pc::Real)
    map_tag == "nofilter" && return LaTeXString("\\mathrm{No\\ filter}")
    return LaTeXString(@sprintf("L_{\\mathrm{eff}}=%.2f\\,\\mathrm{pc}", float(llarge_eff_pc)))
end

"""
    _save_figure(...)

Saves a figure under `<base_out>/figures`.
"""
function _save_figure(cfg::InstrumentalConfig, fig, filename::AbstractString)
    fig_dir = joinpath(cfg.base_out, "figures")
    mkpath(fig_dir)
    out = joinpath(fig_dir, filename)
    save(out, fig)
    @info "Saved figure" path=out
    return out
end

"""
    _robust_colorrange(...)

    Computes robust map color bounds.
"""
function _robust_colorrange(Pmax0, Pmax_filt::AbstractDict{<:Real, <:AbstractMatrix}, L_ok::Vector{Float64})
    vals = vec(float.(Pmax0))
    for Llarge in L_ok
        append!(vals, vec(float.(Pmax_filt[Llarge])))
    end
    vals = vals[isfinite.(vals)]
    isempty(vals) && error("No finite values available to compute color range")
    sort!(vals)
    lo = vals[clamp(round(Int, 0.01 * length(vals)), 1, length(vals))]
    hi = vals[clamp(round(Int, 0.99 * length(vals)), 1, length(vals))]
    return (lo, hi)
end

"""
    _collect_series(...)

    Assembles spectra series descriptors.
"""
function _collect_series(no_filter_img, filtered_imgs::AbstractDict{<:Real, <:AbstractMatrix},
                         kx::AbstractVector, Lsorted::Vector{Float64};
                         nofilter_label::LaTeXString=LaTeXString("\\mathrm{no\\ filter}"))
    cols = curve_colors(1 + length(Lsorted))
    series = Vector{NamedTuple}()

    kx0, P0 = psd1d_x_mean_over_y(no_filter_img, kx; remove_mean=true)
    ok0 = isfinite.(P0) .& (kx0 .> 0) .& (P0 .> 0)
    push!(series, (label=nofilter_label, color=cols[1], kx=kx0, Pk=P0, ok=ok0))

    for (i, Llarge) in enumerate(Lsorted)
        haskey(filtered_imgs, Llarge) || continue
        kxv, Pk = psd1d_x_mean_over_y(filtered_imgs[Llarge], kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)
        push!(series, (
            label=LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"),
            color=cols[i + 1],
            kx=kxv,
            Pk=Pk,
            ok=ok,
        ))
    end
    return series
end

"""
    _select_display_filters(...)

Selects at most `n_display` representative filter scales for the main spectrum plot.
"""
function _select_display_filters(Lsorted::AbstractVector{<:Real}; n_display::Int=4)
    n = length(Lsorted)
    n == 0 && return Float64[]
    n <= n_display && return Float64.(Lsorted)

    idx = unique(clamp.(round.(Int, range(1, n; length=n_display)), 1, n))
    return Float64.(Lsorted[idx])
end

"""
    _cycled_curve_colors(...)

Returns `n` colors by cycling through the base spectrum palette.
"""
function _cycled_curve_colors(n::Int)
    n <= 0 && return Symbol[]
    base = curve_colors(min(n, 9))
    n <= length(base) && return base
    return [base[mod1(i, length(base))] for i in 1:n]
end

"""
    _rm_log_progress_enabled(...)

Controls RM synthesis progress logging. Progress logging is enabled by default
outside CI, with optional explicit override via
`DEPOL_RM_LOG_PROGRESS=1|0`.
"""
function _rm_log_progress_enabled()
    raw = lowercase(strip(get(ENV, "DEPOL_RM_LOG_PROGRESS", "")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    return isempty(get(ENV, "CI", ""))
end
