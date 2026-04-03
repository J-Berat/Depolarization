"""
    _plot_pmax_maps(...)

    Plots no-filter vs filtered `Pmax` maps.
"""
function _sky_plane_labels(cfg::InstrumentalConfig)
    if cfg.los == "x"
        return (LaTeXString("y\\,[\\mathrm{pc}]"), LaTeXString("z\\,[\\mathrm{pc}]"))
    elseif cfg.los == "y"
        return (LaTeXString("z\\,[\\mathrm{pc}]"), LaTeXString("x\\,[\\mathrm{pc}]"))
    else
        return (LaTeXString("x\\,[\\mathrm{pc}]"), LaTeXString("y\\,[\\mathrm{pc}]"))
    end
end

function _cell_edges_for_plot(vals::AbstractVector)
    n = length(vals)
    n >= 2 || error("Need at least two sample locations to build cell edges, got n=$n")
    Δ = vals[2] - vals[1]
    return collect(range(vals[1] - 0.5 * Δ, vals[end] + 0.5 * Δ; length=n + 1))
end

function _finite_mean_map(A::AbstractMatrix)
    acc = 0.0
    n = 0
    @inbounds for val in A
        if isfinite(val)
            acc += Float64(val)
            n += 1
        end
    end
    n == 0 && return NaN
    return acc / n
end

function _bperp_arrow_unit_vector(B1::AbstractMatrix, B2::AbstractMatrix)
    # Plot coordinates use horizontal = second sky-plane component and
    # vertical = first sky-plane component.
    vx = _finite_mean_map(B2)
    vy = _finite_mean_map(B1)
    nrm = hypot(vx, vy)
    if !isfinite(nrm) || nrm <= 0
        return (1.0, 0.0)
    end
    return (vx / nrm, vy / nrm)
end

function _overlay_b_arrow!(ax::Axis, los::String, Lbox_pc::Real, ux::Real, uy::Real)
    margin_pc = 3.8
    arrow_len_pc = 7.2
    head_len_pc = 1.9
    head_open_deg = 28.0

    x0 = margin_pc
    y0 = 2.2

    if los == "x"
        text!(ax, x0 + 0.2, y0 + 0.3;
            text=LaTeXString("\\odot\\;\\langle B \\rangle"),
            align=(:left, :bottom),
            color=:white,
            fontsize=30,
        )
        return nothing
    end

    x1 = x0 + arrow_len_pc * ux
    y1 = y0 + arrow_len_pc * uy

    α = atan(uy, ux)
    β1 = α + π - deg2rad(head_open_deg)
    β2 = α + π + deg2rad(head_open_deg)
    xh1 = x1 + head_len_pc * cos(β1)
    yh1 = y1 + head_len_pc * sin(β1)
    xh2 = x1 + head_len_pc * cos(β2)
    yh2 = y1 + head_len_pc * sin(β2)

    for (lw, col) in ((6.0, (:black, 0.85)), (3.2, (:white, 0.98)))
        lines!(ax, Point2f[(x0, y0), (x1, y1)]; color=col, linewidth=lw)
        lines!(ax, Point2f[(x1, y1), (xh1, yh1)]; color=col, linewidth=lw)
        lines!(ax, Point2f[(x1, y1), (xh2, yh2)]; color=col, linewidth=lw)
    end
    text!(ax, x0 + 0.6, y0 + 0.3;
        text=LaTeXString("\\langle B \\rangle"),
        align=(:left, :bottom),
        color=:white,
        fontsize=30,
    )
    return nothing
end

function _b_field_segments(B1::AbstractMatrix, B2::AbstractMatrix, axes; step::Integer=20, segment_length_pc::Real=1.0)
    x_pc = collect(Float64.(axes.x_pc))
    y_pc = collect(Float64.(axes.y_pc))
    nx, ny = size(B1)
    size(B2) == (nx, ny) || error("B1/B2 size mismatch: B1=$(size(B1)) B2=$(size(B2))")

    segments = Point2f[]
    sizehint!(segments, 2 * cld(nx, step) * cld(ny, step))
    half_length = segment_length_pc / 2

    for i0 in 1:step:nx
        i1 = min(i0 + step - 1, nx)
        for j0 in 1:step:ny
            j1 = min(j0 + step - 1, ny)
            vx = 0.0
            vy = 0.0
            nacc = 0
            @inbounds for j in j0:j1, i in i0:i1
                b1 = Float64(B1[i, j])
                b2 = Float64(B2[i, j])
                isfinite(b1) && isfinite(b2) || continue
                vx += b2
                vy += b1
                nacc += 1
            end
            nacc == 0 && continue
            nrm = hypot(vx, vy)
            nrm > 0 || continue

            ux = vx / nrm
            uy = vy / nrm
            ic = clamp(round(Int, 0.5 * (i0 + i1)), 1, nx)
            jc = clamp(round(Int, 0.5 * (j0 + j1)), 1, ny)
            xc = x_pc[ic]
            yc = y_pc[jc]

            push!(segments, Point2f(xc - half_length * ux, yc - half_length * uy))
            push!(segments, Point2f(xc + half_length * ux, yc + half_length * uy))
        end
    end
    return segments
end

function _overlay_b_field_segments!(ax::Axis, B1::AbstractMatrix, B2::AbstractMatrix, axes; step::Integer=20, segment_length_pc::Real=1.0)
    segments = _b_field_segments(B1, B2, axes; step=step, segment_length_pc=segment_length_pc)
    isempty(segments) && return nothing
    linesegments!(ax, segments; color=(:white, 0.72), linewidth=1.5)
    return nothing
end

function _pmax_canal_segments(Pmap::AbstractMatrix, axes; step::Integer=1, segment_length_pc::Real=1.15)
    x_pc = collect(Float64.(axes.x_pc))
    y_pc = collect(Float64.(axes.y_pc))
    Δx = x_pc[2] - x_pc[1]
    Δy = y_pc[2] - y_pc[1]

    θcanal, mask, score, _ = _canal_orientation_map(Pmap, Δx, Δy)
    nx, ny = size(Pmap)

    segments = Point2f[]
    sizehint!(segments, 2 * cld(nx, step) * cld(ny, step))
    half_length = segment_length_pc / 2

    for i0 in 1:step:nx
        i1 = min(i0 + step - 1, nx)
        for j0 in 1:step:ny
            j1 = min(j0 + step - 1, ny)
            best_score = -Inf
            best_i = 0
            best_j = 0
            @inbounds for j in j0:j1, i in i0:i1
                mask[i, j] || continue
                sij = score[i, j]
                if sij > best_score
                    best_score = sij
                    best_i = i
                    best_j = j
                end
            end
            best_i == 0 && continue

            θ = θcanal[best_i, best_j]
            isfinite(θ) || continue

            ux = cos(θ)
            uy = sin(θ)
            xc = x_pc[best_i]
            yc = y_pc[best_j]

            push!(segments, Point2f(xc - half_length * ux, yc - half_length * uy))
            push!(segments, Point2f(xc + half_length * ux, yc + half_length * uy))
        end
    end

    return segments
end

function _pmax_zone_mask(Pmap::AbstractMatrix; min_area_pix::Integer=12)
    score = _canal_score_map(Float64.(Pmap))
    Ds = _box_mean(score, _CANAL_RIDGE_SMOOTH_RADIUS_PIX)
    thr = _finite_quantile(Ds, _CANAL_SCORE_QUANTILE)
    mask = isfinite.(Ds) .& (Ds .>= thr)
    comps = _zone_connected_components(mask)
    filtered_mask = falses(size(mask))
    for comp in comps
        length(comp) >= min_area_pix || continue
        for I in comp
            filtered_mask[I] = true
        end
    end
    return filtered_mask
end

function _pmax_zone_display_alpha(Pmap::AbstractMatrix)
    base_mask = _pmax_zone_mask(Pmap)
    any(base_mask) || return fill(NaN, size(Pmap))

    score = _canal_score_map(Float64.(Pmap))
    Ds = _box_mean(score, _CANAL_RIDGE_SMOOTH_RADIUS_PIX)
    display_thr = _finite_quantile(Ds, 0.82)

    # Keep the measured zones, but expand the display by about one pixel so
    # the white patches stay well defined on the valleys.
    support = _box_mean(Float64.(base_mask), 1)
    candidate = isfinite.(Ds) .& (support .> 0.0) .& ((Ds .>= display_thr) .| base_mask)

    alpha = Float64.(candidate)
    alpha = _box_mean(alpha, 1)
    alpha .*= _box_mean(Float64.(candidate), 1)
    maxα = maximum(alpha)
    maxα > 0 && (alpha ./= maxα)
    alpha[alpha .< 0.18] .= NaN
    return alpha
end

function _overlay_zone_contours!(ax::Axis, Pmap::AbstractMatrix, axes)
    zalpha = _pmax_zone_display_alpha(Pmap)
    any(isfinite, zalpha) || return nothing
    x_edges = _cell_edges_for_plot(collect(Float64.(axes.x_pc)))
    y_edges = _cell_edges_for_plot(collect(Float64.(axes.y_pc)))
    zone_cmap = cgrad([
        RGBAf(1, 1, 1, 0.0),
        RGBAf(1, 1, 1, 0.16),
        RGBAf(1, 1, 1, 0.30),
        RGBAf(1, 1, 1, 0.42),
    ])
    heatmap!(ax, x_edges, y_edges, zalpha; colormap=zone_cmap, colorrange=(0, 1))
    return nothing
end

function _overlay_canal_segments!(ax::Axis, Pmap::AbstractMatrix, axes)
    segments = _pmax_canal_segments(Pmap, axes)
    isempty(segments) && return nothing
    linesegments!(ax, segments; color=(:black, 0.82), linewidth=5.0)
    linesegments!(ax, segments; color=(:white, 0.97), linewidth=2.6)
    return nothing
end

function _plot_pmax_maps(cfg::InstrumentalConfig, data::FilterPassResult)
    L_ok = data.L_ok
    axes = data.axes
    xlabel_txt, ylabel_txt = _sky_plane_labels(cfg)
    x_edges = _cell_edges_for_plot(collect(axes.x_pc))
    y_edges = _cell_edges_for_plot(collect(axes.y_pc))
    pmax_colorrange = _robust_colorrange(data.Pmax0, data.Pmax_filt, L_ok)
    b1, b2, _ = _project_bperp_maps(cfg)
    arrow_ux, arrow_uy = _bperp_arrow_unit_vector(b1, b2)
    seg_step = cfg.los == "x" ? 28 : 24
    seg_length_pc = cfg.los == "x" ? 0.9 : 1.0

    with_theme(theme_latexfonts()) do
        set_theme_heatmaps!()

        fig = Figure(size=(3000, 1360), figure_padding=(18, 18, 10, 10))

        ticklabel_size = 40
        axis_label_size = 44
        function _panel_axis(row::Int, col::Int; show_xlabel::Bool=false, show_ylabel::Bool=false,
                show_xticklabels::Bool=false, show_yticklabels::Bool=false)
            return Axis(fig[row, col],
                xlabel=show_xlabel ? xlabel_txt : "",
                ylabel=show_ylabel ? ylabel_txt : "",
                xlabelsize=axis_label_size,
                ylabelsize=axis_label_size,
                aspect=DataAspect(),
                xaxisposition=:bottom,
                xticklabelsize=ticklabel_size,
                yticklabelsize=ticklabel_size,
                xticklabelsvisible=show_xticklabels,
                yticklabelsvisible=show_yticklabels,
                limits=(0.0, cfg.Lbox_pc, 0.0, cfg.Lbox_pc),
                xautolimitmargin=(0.0f0, 0.0f0),
                yautolimitmargin=(0.0f0, 0.0f0),
            )
        end

        ax0_top = _panel_axis(1, 1; show_ylabel=false, show_yticklabels=false)
        hm0 = heatmap!(ax0_top, x_edges, y_edges, data.Pmax0; colorrange=pmax_colorrange, colormap=:magma)
        _overlay_zone_contours!(ax0_top, data.Pmax0, axes)
        _overlay_b_arrow!(ax0_top, cfg.los, cfg.Lbox_pc, arrow_ux, arrow_uy)
        set_latex_linear_ticks_pc!(ax0_top; step_pc=10.0, Lbox_pc=cfg.Lbox_pc)
        ax0_top.xticklabelsvisible = false
        ax0_top.yticklabelsvisible = false

        text!(ax0_top, 0.95, 0.08;
            text=LaTeXString("\\mathrm{no\\ filter}"),
            space=:relative,
            align=(:right, :bottom),
            color=:white,
            fontsize=36,
        )

        ax0_bot = _panel_axis(2, 1; show_xlabel=true, show_ylabel=true, show_xticklabels=true, show_yticklabels=true)
        heatmap!(ax0_bot, x_edges, y_edges, data.Pmax0; colorrange=pmax_colorrange, colormap=:magma)
        _overlay_b_field_segments!(ax0_bot, b1, b2, axes; step=seg_step, segment_length_pc=seg_length_pc)
        set_latex_linear_ticks_pc!(ax0_bot; step_pc=10.0, Lbox_pc=cfg.Lbox_pc)

        for (j, Llarge) in enumerate(L_ok)
            ax_top = _panel_axis(1, j + 1)
            heatmap!(ax_top, x_edges, y_edges, data.Pmax_filt[Llarge]; colorrange=pmax_colorrange, colormap=:magma)
            _overlay_zone_contours!(ax_top, data.Pmax_filt[Llarge], axes)
            _overlay_b_arrow!(ax_top, cfg.los, cfg.Lbox_pc, arrow_ux, arrow_uy)
            set_latex_linear_ticks_pc!(ax_top; step_pc=10.0, Lbox_pc=cfg.Lbox_pc)
            ax_top.xticklabelsvisible = false
            ax_top.yticklabelsvisible = false

            text!(ax_top, 0.95, 0.08;
                text=LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"),
                space=:relative,
                align=(:right, :bottom),
                color=:white,
                fontsize=36,
            )

            ax_bot = _panel_axis(2, j + 1)
            heatmap!(ax_bot, x_edges, y_edges, data.Pmax_filt[Llarge]; colorrange=pmax_colorrange, colormap=:magma)
            _overlay_b_field_segments!(ax_bot, b1, b2, axes; step=seg_step, segment_length_pc=seg_length_pc)
            set_latex_linear_ticks_pc!(ax_bot; step_pc=10.0, Lbox_pc=cfg.Lbox_pc)
        end

        Colorbar(fig[1:2, length(L_ok) + 2], hm0;
            label=LaTeXString("P_{\\max}\\,[\\mathrm{K}]"),
            tellheight=true,
            tickformat=cb_latex_ticks,
            ticklabelsize=ticklabel_size,
            labelsize=axis_label_size,
        )

        for j in 1:(length(L_ok) + 1)
            colsize!(fig.layout, j, Fixed(360))
        end
        for r in 1:2
            rowsize!(fig.layout, r, Fixed(360))
        end
        colsize!(fig.layout, length(L_ok) + 2, Fixed(54))
        rowgap!(fig.layout, 8)
        colgap!(fig.layout, length(L_ok) + 1, 5)
        resize_to_layout!(fig)
        _save_figure(cfg, fig, "pmax_maps.pdf")
        display(fig)
    end
end

"""
    _plot_psd_panels(...)

    Plots isotropic PSD panels.
"""
function _plot_psd_panels(cfg::InstrumentalConfig, data::FilterPassResult)
    L_ok = data.L_ok
    axes = data.axes

    with_theme(theme_latexfonts()) do
        set_theme_spectra!()

        figPSD = Figure(size=(3000, 950), figure_padding=30)

        nbins_psd = 90
        kmin_plot = 1 / cfg.Lbox_pc
        kmax_plot = maximum(sqrt.(repeat(axes.kx, 1, cfg.m).^2 .+ repeat(axes.ky', cfg.n, 1).^2))

        kcen0, Pk0 = psd1d_isotropic(data.Pmax0, axes.kx, axes.ky; nbins=nbins_psd, kmin=kmin_plot, kmax=kmax_plot)
        ok0 = isfinite.(Pk0) .& (kcen0 .> 0) .& (Pk0 .> 0)

        axpsd0 = Axis(figPSD[1, 1],
            xlabel=LaTeXString("k\\,[\\mathrm{pc^{-1}}]"),
            ylabel=LaTeXString("S_{P}(k)"),
        )
        lines!(axpsd0, kcen0[ok0], Pk0[ok0])
        configure_axis_style!(axpsd0;
            xdata=kcen0[ok0],
            ydata=Pk0[ok0],
            xscale=:log10,
            yscale=:log10,
            xminor_subdiv=9,
            yminor_subdiv=9,
        )

        for (j, Llarge) in enumerate(L_ok)
            kcen, Pk = psd1d_isotropic(data.Pmax_filt[Llarge], axes.kx, axes.ky; nbins=nbins_psd, kmin=kmin_plot, kmax=kmax_plot)
            ok = isfinite.(Pk) .& (kcen .> 0) .& (Pk .> 0)

            axpsd = Axis(figPSD[1, j + 1],
                xlabel=LaTeXString("k\\,[\\mathrm{pc^{-1}}]"),
                ylabelvisible=false,
                yticklabelsvisible=false,
            )
            lines!(axpsd, kcen[ok], Pk[ok])
            configure_axis_style!(axpsd;
                xdata=kcen[ok],
                ydata=Pk[ok],
                xscale=:log10,
                yscale=:log10,
                xminor_subdiv=9,
                yminor_subdiv=9,
            )
        end

        _save_figure(cfg, figPSD, "psd_panels.pdf")
        display(figPSD)
    end
end

function _powerlaw_slope_fit(k::AbstractVector, P::AbstractVector; kmax::Real)
    mask = isfinite.(k) .& isfinite.(P) .& (k .> 0) .& (P .> 0) .& (k .< float(kmax))
    count(mask) >= 3 || return nothing

    x = log10.(Float64.(k[mask]))
    y = log10.(Float64.(P[mask]))
    n = length(x)

    x̄ = sum(x) / n
    ȳ = sum(y) / n
    sxx = sum((x .- x̄) .^ 2)
    sxx > 0 || return nothing
    sxy = sum((x .- x̄) .* (y .- ȳ))
    α = sxy / sxx
    β = ȳ - α * x̄

    kvals = Float64.(k[mask])
    xlo = quantile(kvals, 0.20)
    xhi = quantile(kvals, 0.80)
    xhi > xlo || return nothing
    ylo = 10.0 ^ (β + α * log10(xlo))
    yhi = 10.0 ^ (β + α * log10(xhi))
    xmid = sqrt(xlo * xhi)
    ymid = 10.0 ^ (β + α * log10(xmid))

    return (alpha=α, beta=β, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, xmid=xmid, ymid=ymid)
end

function _add_lowk_slope_guides!(ax::Axis, series; kmax::Real)
    isempty(series) && return nothing
    s = series[1]
    fit = _powerlaw_slope_fit(s.kx, s.Pk; kmax=kmax)
    fit === nothing && return nothing

    yshift = 3.2
    lines!(ax, [fit.xlo, fit.xhi], [fit.ylo * yshift, fit.yhi * yshift];
        color=(s.color, 0.95),
        linewidth=4.5,
        linestyle=:dash,
    )

    text!(ax, fit.xmid * 1.04, fit.ymid * 5.1;
        text=LaTeXString("\\mathrm{slope}=$(round(fit.alpha; digits=2))"),
        color=:black,
        fontsize=32,
        align=(:left, :bottom),
    )
    return nothing
end

"""
    _plot_pmax_kx(...)

    Plots `Pmax` spectra along `kx` with peak inset.
"""
function _plot_pmax_kx(cfg::InstrumentalConfig, data::FilterPassResult)
    Lall = sort(data.L_ok)
    Ldisplay = _select_display_filters(Lall; n_display=4)
    series = _collect_series(data.Pmax0, Dict(k => data.Pmax_filt[k] for k in Ldisplay), data.axes.kx, Ldisplay)
    kmin_win = 0.0191
    kmax_win = Inf

    with_theme(theme_latexfonts()) do
        set_theme_spectra!()

        figS = Figure(size=(1500, 1200), figure_padding=30)
        axS = Axis(figS[1, 1],
            xlabel=LaTeXString("k_{\\parallel}\\,[\\mathrm{pc^{-1}}]"),
            ylabel=LaTeXString("S_{P}(k_{\\parallel})"),
        )

        plot_multi_spectrum!(axS, series;
            add_verticals=false,
            Llist_pc=cfg.Llarge_list,
            peak_window=(kmin_win, kmax_win),
        )
        _add_lowk_slope_guides!(axS, series; kmax=4.8e-3)

        let
            @info "Inset Pmax window configured" kmin_rad_pc_inv=kmin_win kmax_rad_pc_inv=kmax_win

            Leff = Float64[]
            kpeak = Float64[]
            cpeak = Symbol[]

            cols = _cycled_curve_colors(1 + length(Lall))
            for (i, Llarge) in enumerate(Lall)
                kxv, Pk = psd1d_x_mean_over_y(data.Pmax_filt[Llarge], data.axes.kx; remove_mean=false)
                kp = kpeak_in_window(kxv, Pk; kmin=kmin_win, kmax=kmax_win, smooth_radius=2)
                push!(Leff, (50 / 256) * Llarge)
                push!(kpeak, kp)
                push!(cpeak, cols[i + 1])
            end

            ok = isfinite.(Leff) .& isfinite.(kpeak) .& (kpeak .> 0)

            axInset = Axis(figS[1, 1],
                width=Relative(0.32),
                height=Relative(0.28),
                halign=0.17,
                valign=0.18,
                xlabel=LaTeXString("L_{\\mathrm{large}}\\,[\\mathrm{pc}]"),
                ylabel=LaTeXString("k_{\\parallel,\\mathrm{peak}}"),
                xlabelsize=36,
                ylabelsize=36,
                xticksize=8,
                yticksize=8,
                xtickwidth=2,
                ytickwidth=2,
                xticklabelsize=34,
                yticklabelsize=34,
                xtickalign=1,
                ytickalign=1,
            )

            lines!(axInset, Leff[ok], kpeak[ok], color=:black, linewidth=2)
            scatter!(axInset, Leff[ok], kpeak[ok];
                color=cpeak[ok],
                markersize=18,
                strokecolor=:black,
                strokewidth=1.5,
            )

            xlo = minimum(Leff[ok])
            xhi = maximum(Leff[ok])
            xpad = max(0.05 * (xhi - xlo), 0.5)
            xlims!(axInset, max(0.0, xlo - xpad), xhi + xpad)

            ymin = max(minimum(kpeak[ok]) * 0.85, eps(Float64))
            ymax = maximum(kpeak[ok]) * 1.15
            configure_axis_style!(axInset;
                xdata=Leff[ok],
                ydata=kpeak[ok],
                xscale=:linear,
                yscale=:linear,
                xminor_subdiv=0,
                yminor_subdiv=0,
                xaxisposition=:bottom,
                yaxisposition=:left,
                x_linear_digits=1,
                y_linear_digits=3,
            )
            ylims!(axInset, ymin, ymax)
            axInset.backgroundcolor = (:white, 0.94)
        end

        _save_figure(cfg, figS, "pmax_kx.pdf")
        display(figS)
    end
end

"""
    _plot_component_spectrum(...)

    Generic component-spectrum plotting routine.
"""
function _plot_component_spectrum(cfg::InstrumentalConfig, kx::AbstractVector,
                                  no_filter_img::AbstractMatrix,
                                  filtered_imgs::AbstractDict{<:Real, <:AbstractMatrix},
                                  ylabel::LaTeXString;
                                  add_verticals::Bool=true,
                                  peak_window::Union{Nothing, Tuple{Real, Real}}=nothing,
                                  inset_title::Union{Nothing, Symbol}=nothing,
                                  save_tag::AbstractString="spectrum")
    Lall = sort(collect(keys(filtered_imgs)))
    Ldisplay = _select_display_filters(Lall; n_display=4)
    series = _collect_series(no_filter_img, filtered_imgs, kx, Ldisplay)

    with_theme(theme_latexfonts()) do
        set_theme_spectra!()

        fig = Figure(size=(1500, 1200), figure_padding=30)
        ax = Axis(fig[1, 1],
            xlabel=LaTeXString("k_{\\parallel}\\,[\\mathrm{rad\\,pc^{-1}}]"),
            ylabel=ylabel,
        )

        plot_multi_spectrum!(ax, series;
            add_verticals=add_verticals,
            Llist_pc=cfg.Llarge_list,
            peak_window=peak_window,
        )

        if inset_title !== nothing && peak_window !== nothing
            kmin_win, kmax_win = peak_window
            Leff = Float64[]
            kpeak = Float64[]
            cpeak = Symbol[]
            cols = _cycled_curve_colors(1 + length(Lall))

            for (i, Llarge) in enumerate(Lall)
                kxv, Pk = psd1d_x_mean_over_y(filtered_imgs[Llarge], kx; remove_mean=false)
                kp = kpeak_in_window(kxv, Pk; kmin=kmin_win, kmax=kmax_win, smooth_radius=2)
                push!(Leff, (50 / 256) * Llarge)
                push!(kpeak, kp)
                push!(cpeak, cols[i + 1])
                @info "Inset peak located" inset=String(inset_title) Llarge_pc=Llarge kx_peak_rad_pc_inv=kp
            end

            ok = isfinite.(Leff) .& isfinite.(kpeak) .& (kpeak .> 0)
            axInset = Axis(fig[1, 1],
                width=Relative(0.38),
                height=Relative(0.34),
                halign=:left,
                valign=:bottom,
                xlabel=LaTeXString("L_{\\mathrm{large}}\\,[\\mathrm{pc}]"),
                ylabel=LaTeXString("k_{\\parallel,\\mathrm{peak}}"),
                xlabelsize=36,
                ylabelsize=36,
                xticksize=8,
                yticksize=8,
                xtickwidth=2,
                ytickwidth=2,
                xticklabelsize=34,
                yticklabelsize=34,
                xtickalign=1,
                ytickalign=1,
            )

            lines!(axInset, Leff[ok], kpeak[ok], color=:black, linewidth=2)
            scatter!(axInset, Leff[ok], kpeak[ok];
                color=cpeak[ok],
                markersize=18,
                strokecolor=:black,
                strokewidth=1.5,
            )

            xlo = minimum(Leff[ok])
            xhi = maximum(Leff[ok])
            xpad = max(0.05 * (xhi - xlo), 0.5)
            xlims!(axInset, max(0.0, xlo - xpad), xhi + xpad)

            ymin = max(minimum(kpeak[ok]) * 0.85, eps(Float64))
            ymax = maximum(kpeak[ok]) * 1.15
            configure_axis_style!(axInset;
                xdata=Leff[ok],
                ydata=kpeak[ok],
                xscale=:linear,
                yscale=:linear,
                xminor_subdiv=0,
                yminor_subdiv=0,
                xaxisposition=:top,
                yaxisposition=:right,
                x_linear_digits=1,
                y_linear_digits=3,
            )
            ylims!(axInset, ymin, ymax)
            axInset.backgroundcolor = (:white, 0.85)
        end

        _save_figure(cfg, fig, "$(save_tag).pdf")
        display(fig)
    end
end

"""
    _build_filtered_dict_from_slice(...)

    Builds per-filter dictionary from a 2D slice.
"""
function _build_filtered_dict_from_slice(slice0::AbstractMatrix, cfg::InstrumentalConfig, scales)
    size(slice0) == (cfg.n, cfg.m) || error("Input slice must have size ($(cfg.n),$(cfg.m)), got $(size(slice0))")

    out = Dict{Float64, Matrix{Float64}}()
    for Llarge in Float64.(cfg.Llarge_list)
        H, _ = instrument_bandpass_L(cfg.n, cfg.m;
                                     Δx=scales.Δx,
                                     Δy=scales.Δy,
                                     Lcut_small=cfg.Lcut_small,
                                     Llarge=Llarge,
                                     fNy=scales.fNy)
        out[Llarge] = apply_instrument_2d(slice0, H)
    end
    return out
end

"""
    _validate_phi_cubes(...)

    Validates `Qphi/Uphi` cubes against phi axis and configured grid.
"""
function _validate_phi_cubes(cfg::InstrumentalConfig, data::FilterPassResult, Qphi_cube, Uphi_cube)
    nϕ = length(data.PhiArray)
    nϕ > 0 || error("Phi axis is empty")
    require_channel_index(cfg.iphi, nϕ, "iphi")

    q_layout = require_stokes_cube_layout(Qphi_cube, "realFDF", cfg, nϕ)
    u_layout = require_stokes_cube_layout(Uphi_cube, "imagFDF", cfg, nϕ)
    q_layout == u_layout || error("realFDF/imagFDF layout mismatch: realFDF frequency axis dim=$q_layout, imagFDF frequency axis dim=$u_layout")

    return true
end
