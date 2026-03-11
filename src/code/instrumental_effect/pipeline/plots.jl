"""
    _plot_pmax_maps(...)

    Plots no-filter vs filtered `Pmax` maps.
"""
function _plot_pmax_maps(cfg::InstrumentalConfig, data::FilterPassResult)
    L_ok = data.L_ok
    axes = data.axes
    pmax_colorrange = _robust_colorrange(data.Pmax0, data.Pmax_filt, L_ok)

    with_theme(theme_latexfonts()) do
        set_theme_heatmaps!()

        fig = Figure(size=(3000, 1050), figure_padding=30)

        ax0 = Axis(fig[1, 1],
            xlabel=LaTeXString("x\\,[\\mathrm{pc}]"),
            ylabel=LaTeXString("y\\,[\\mathrm{pc}]"),
            aspect=DataAspect(),
        )
        hm0 = heatmap!(ax0, axes.x_pc, axes.y_pc, data.Pmax0; colorrange=pmax_colorrange, colormap=:magma)
        set_latex_linear_ticks_pc!(ax0; step_pc=10.0, Lbox_pc=cfg.Lbox_pc)

        text!(ax0, 0.97, 0.97;
            text=LaTeXString("\\mathrm{no\\ filter}"),
            space=:relative,
            align=(:right, :top),
            color=:white,
            fontsize=50,
        )

        for (j, Llarge) in enumerate(L_ok)
            ax = Axis(fig[1, j + 1],
                xlabel=LaTeXString("x\\,[\\mathrm{pc}]"),
                aspect=DataAspect(),
                ylabelvisible=false,
                yticklabelsvisible=false,
            )
            heatmap!(ax, axes.x_pc, axes.y_pc, data.Pmax_filt[Llarge]; colorrange=pmax_colorrange, colormap=:magma)
            set_latex_linear_ticks_pc!(ax; step_pc=10.0, Lbox_pc=cfg.Lbox_pc)

            text!(ax, 0.97, 0.97;
                text=LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"),
                space=:relative,
                align=(:right, :top),
                color=:white,
                fontsize=50,
            )
        end

        Colorbar(fig[1, length(L_ok) + 2], hm0;
            label=LaTeXString("P_{\\max}\\,[\\mathrm{K}]"),
            tellheight=true,
            tickformat=cb_latex_ticks,
        )

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
        kmin_plot = 2π * (1 / cfg.Lbox_pc)
        kmax_plot = maximum(sqrt.(repeat(axes.kx, 1, cfg.m).^2 .+ repeat(axes.ky', cfg.n, 1).^2))

        kcen0, Pk0 = psd1d_isotropic(data.Pmax0, axes.kx, axes.ky; nbins=nbins_psd, kmin=kmin_plot, kmax=kmax_plot)
        ok0 = isfinite.(Pk0) .& (kcen0 .> 0) .& (Pk0 .> 0)

        axpsd0 = Axis(figPSD[1, 1],
            xlabel=LaTeXString("k\\,[\\mathrm{rad\\,pc^{-1}}]"),
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
                xlabel=LaTeXString("k\\,[\\mathrm{rad\\,pc^{-1}}]"),
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

"""
    _plot_pmax_kx(...)

    Plots `Pmax` spectra along `kx` with peak inset.
"""
function _plot_pmax_kx(cfg::InstrumentalConfig, data::FilterPassResult)
    Lall = sort(data.L_ok)
    Ldisplay = _select_display_filters(Lall; n_display=4)
    series = _collect_series(data.Pmax0, Dict(k => data.Pmax_filt[k] for k in Ldisplay), data.axes.kx, Ldisplay)

    with_theme(theme_latexfonts()) do
        set_theme_spectra!()

        figS = Figure(size=(1500, 1200), figure_padding=30)
        axS = Axis(figS[1, 1],
            xlabel=LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
            ylabel=LaTeXString("S_{P}(k_x)"),
        )

        plot_multi_spectrum!(axS, series; add_verticals=true, Llist_pc=cfg.Llarge_list)

        let
            kmin_win = 0.12
            kmax_win = Inf
            @info "Inset Pmax window configured" kmin_rad_pc_inv=kmin_win kmax_rad_pc_inv=kmax_win

            Leff = Float64[]
            kpeak = Float64[]
            cpeak = Symbol[]

            cols = _cycled_curve_colors(1 + length(Lall))
            for (i, Llarge) in enumerate(Lall)
                kxv, Pk = psd1d_x_mean_over_y(data.Pmax_filt[Llarge], data.axes.kx; remove_mean=true)
                kp = kpeak_in_window(kxv, Pk; kmin=kmin_win, kmax=kmax_win)
                push!(Leff, (50 / 256) * Llarge)
                push!(kpeak, kp)
                push!(cpeak, cols[i + 1])
            end

            ok = isfinite.(Leff) .& isfinite.(kpeak) .& (kpeak .> 0)

            axInset = Axis(figS[1, 1],
                width=Relative(0.38),
                height=Relative(0.34),
                halign=:left,
                valign=:bottom,
                xlabel=LaTeXString("L_{\\mathrm{large}}\\,[\\mathrm{pc}]"),
                ylabel=LaTeXString("k_{x,\\mathrm{peak}}"),
                xlabelsize=36,
                ylabelsize=36,
                xticksize=8,
                yticksize=8,
                xtickwidth=2,
                ytickwidth=2,
                xticklabelsize=22,
                yticklabelsize=22,
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
            xlabel=LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
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
                kxv, Pk = psd1d_x_mean_over_y(filtered_imgs[Llarge], kx; remove_mean=true)
                kp = kpeak_in_window(kxv, Pk; kmin=kmin_win, kmax=kmax_win)
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
                ylabel=LaTeXString("k_{x,\\mathrm{peak}}"),
                xlabelsize=36,
                ylabelsize=36,
                xticksize=8,
                yticksize=8,
                xtickwidth=2,
                ytickwidth=2,
                xticklabelsize=22,
                yticklabelsize=22,
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
