function _filter_tag(Llarge::Real)
    @sprintf("HardBandPass_remove_L0_to_1pc_and_%dto50pc", Int(round(Llarge)))
end

function _robust_colorrange(Pmax0, Pmax_filt::Dict{Float64, Matrix}, L_ok::Vector{Float64})
    vals = vec(float.(Pmax0))
    for Llarge in L_ok
        append!(vals, vec(float.(Pmax_filt[Llarge])))
    end
    vals = vals[isfinite.(vals)]
    sort!(vals)
    lo = vals[clamp(round(Int, 0.01 * length(vals)), 1, length(vals))]
    hi = vals[clamp(round(Int, 0.99 * length(vals)), 1, length(vals))]
    return (lo, hi)
end

function _collect_series(no_filter_img, filtered_imgs::Dict{Float64, Matrix},
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

function run_filter_pass(cfg::InstrumentalConfig)
    mkpath(cfg.base_out)

    scales = grid_scales(cfg)
    axes = spectral_axes(cfg, scales.Δx, scales.Δy)

    @info "Grid" cfg.n cfg.m scales.Δx scales.Δy scales.fNy

    nuArray, PhiArray = rm_synthesis_axes(cfg)

    Qdata = read_FITS(cfg.Q_in)
    Udata = read_FITS(cfg.U_in)

    Qslice_filt = Dict{Float64, Matrix{Float32}}()
    Uslice_filt = Dict{Float64, Matrix{Float32}}()

    for Llarge in cfg.Llarge_list
        tag = _filter_tag(Llarge)
        subdir = joinpath(cfg.base_out, tag)
        mkpath(subdir)

        H, _ = instrument_bandpass_L(cfg.n, cfg.m;
                                     Δx=scales.Δx, Δy=scales.Δy,
                                     Lcut_small=cfg.Lcut_small,
                                     Llarge=Llarge,
                                     fNy=scales.fNy)

        Qf = apply_to_array_xy(Qdata, H; n=cfg.n, m=cfg.m)
        Uf = apply_to_array_xy(Udata, H; n=cfg.n, m=cfg.m)

        Qslice_filt[Llarge] = Float32.(get_chan_xy(Qf, cfg.ichan, cfg.n, cfg.m))
        Uslice_filt[Llarge] = Float32.(get_chan_xy(Uf, cfg.ichan, cfg.n, cfg.m))

        write_FITS(joinpath(subdir, "Qnu_filtered.fits"), Qf)
        write_FITS(joinpath(subdir, "Unu_filtered.fits"), Uf)

        absF, _, _ = RMSynthesis(Qf, Uf, nuArray, PhiArray; log_progress=true)
        Pmaxf = Pphi_max_map(absF)

        mkpath(joinpath(subdir, "RMSynthesis"))
        write_FITS(joinpath(subdir, "RMSynthesis", "Pphi_max.fits"), Pmaxf)
    end

    Pmax0 = read_FITS(cfg.Pmax_nofilter_path)

    Pmax_filt = Dict{Float64, Matrix}()
    L_ok = Float64[]
    for Llarge in cfg.Llarge_list
        p = joinpath(cfg.base_out, _filter_tag(Llarge), "RMSynthesis", "Pphi_max.fits")
        isfile(p) || continue
        Pmax_filt[Llarge] = read_FITS(p)
        push!(L_ok, Llarge)
    end

    @assert !isempty(L_ok) "No filtered Pphi_max.fits found under: $(cfg.base_out)"

    return (
        Qdata=Qdata,
        Udata=Udata,
        Qslice_filt=Qslice_filt,
        Uslice_filt=Uslice_filt,
        Pmax0=Pmax0,
        Pmax_filt=Pmax_filt,
        L_ok=L_ok,
        scales=scales,
        axes=axes,
        nuArray=nuArray,
        PhiArray=PhiArray,
    )
end

function _plot_pmax_maps(cfg::InstrumentalConfig, data)
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

        display(fig)
    end
end

function _plot_psd_panels(cfg::InstrumentalConfig, data)
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
            ylabel=LaTeXString("\\langle |\\tilde{P}|^2 \\rangle(k)"),
        )
        lines!(axpsd0, kcen0[ok0], Pk0[ok0])
        set_log_safe!(axpsd0; xdata=kcen0[ok0], ydata=Pk0[ok0])
        autolimits!(axpsd0)
        set_pow10_ticks!(axpsd0; which=:both)

        for (j, Llarge) in enumerate(L_ok)
            kcen, Pk = psd1d_isotropic(data.Pmax_filt[Llarge], axes.kx, axes.ky; nbins=nbins_psd, kmin=kmin_plot, kmax=kmax_plot)
            ok = isfinite.(Pk) .& (kcen .> 0) .& (Pk .> 0)

            axpsd = Axis(figPSD[1, j + 1],
                xlabel=LaTeXString("k\\,[\\mathrm{rad\\,pc^{-1}}]"),
                ylabelvisible=false,
                yticklabelsvisible=false,
            )
            lines!(axpsd, kcen[ok], Pk[ok])
            set_log_safe!(axpsd; xdata=kcen[ok], ydata=Pk[ok])
            autolimits!(axpsd)
            set_pow10_ticks!(axpsd; which=:both)
        end

        display(figPSD)
    end
end

function _plot_pmax_kx(cfg::InstrumentalConfig, data)
    Lsorted = sort(data.L_ok)
    series = _collect_series(data.Pmax0, Dict(k => data.Pmax_filt[k] for k in Lsorted), data.axes.kx, Lsorted)

    with_theme(theme_latexfonts()) do
        set_theme_spectra!()

        figS = Figure(size=(1500, 1200), figure_padding=30)
        axS = Axis(figS[1, 1],
            xlabel=LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
            ylabel=LaTeXString("\\langle |\\tilde{P}(k_x, y)|^2 \\rangle_y"),
        )

        plot_multi_spectrum!(axS, series; add_verticals=true, Llist_pc=cfg.Llarge_list)

        let
            kmin_win = 0.12
            kmax_win = Inf
            @info "Inset Pmax: k-window" kmin_win kmax_win

            Leff = Float64[]
            kpeak = Float64[]
            cpeak = Symbol[]

            for (i, Llarge) in enumerate(Lsorted)
                kxv, Pk = psd1d_x_mean_over_y(data.Pmax_filt[Llarge], data.axes.kx; remove_mean=true)
                kp = kpeak_in_window(kxv, Pk; kmin=kmin_win, kmax=kmax_win)
                push!(Leff, (50 / 256) * Llarge)
                push!(kpeak, kp)
                push!(cpeak, curve_colors(1 + length(Lsorted))[i + 1])
            end

            ok = isfinite.(Leff) .& isfinite.(kpeak) .& (kpeak .> 0)

            axInset = Axis(figS[1, 1],
                width=Relative(0.38),
                height=Relative(0.34),
                halign=:left,
                valign=:bottom,
                xlabel=LaTeXString("L_{\\mathrm{large}}"),
                ylabel=LaTeXString("k_{\\mathrm{x,peak}}"),
                xminorticksvisible=false,
                yminorticksvisible=true,
                yminorticks=IntervalsBetween(9),
                yminorticksize=5,
                yminortickwidth=1.5,
                xlabelsize=25,
                ylabelsize=25,
                xticksize=8,
                yticksize=8,
                xtickwidth=2,
                ytickwidth=2,
                xticklabelsize=22,
                yticklabelsize=22,
                xtickalign=1,
                ytickalign=1,
                xaxisposition=:top,
                yaxisposition=:right,
            )

            lines!(axInset, Leff[ok], kpeak[ok], color=:black, linewidth=2)
            scatter!(axInset, Leff[ok], kpeak[ok]; color=cpeak[ok], markersize=14)
            xlims!(axInset, 0, 40)

            ymin = max(minimum(kpeak[ok]) * 0.9, eps(Float64))
            ymax = maximum(kpeak[ok]) * 1.1
            axInset.yscale = linear
            ylims!(axInset, ymin, ymax)

            axInset.xtickformat = (vs) -> [LaTeXString("\\mathrm{$(round(v; digits=1))}") for v in vs]
            set_pow10_ticks!(axInset; which=:y)
            axInset.backgroundcolor = (:white, 0.85)
        end

        display(figS)
    end
end

function _plot_component_spectrum(cfg::InstrumentalConfig, kx::AbstractVector,
                                  no_filter_img::AbstractMatrix,
                                  filtered_imgs::Dict{Float64, Matrix},
                                  ylabel::LaTeXString;
                                  add_verticals::Bool=true,
                                  peak_window::Union{Nothing, Tuple{Real, Real}}=nothing,
                                  inset_title::Union{Nothing, Symbol}=nothing)
    Lsorted = sort(collect(keys(filtered_imgs)))
    series = _collect_series(no_filter_img, filtered_imgs, kx, Lsorted)

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
            cols = curve_colors(1 + length(Lsorted))

            for (i, Llarge) in enumerate(Lsorted)
                s = series[i + 1]
                kp = kpeak_in_window(s.kx, s.Pk; kmin=kmin_win, kmax=kmax_win)
                push!(Leff, (50 / 256) * Llarge)
                push!(kpeak, kp)
                push!(cpeak, cols[i + 1])
                @info "Inset peak" inset_title Llarge kp
            end

            ok = isfinite.(Leff) .& isfinite.(kpeak) .& (kpeak .> 0)
            axInset = Axis(fig[1, 1],
                width=Relative(0.38),
                height=Relative(0.34),
                halign=:left,
                valign=:bottom,
                xlabel=LaTeXString("L_{\\mathrm{large}}"),
                ylabel=LaTeXString("k_{\\mathrm{x,peak}}"),
                xminorticksvisible=false,
                yminorticksvisible=true,
                yminorticks=IntervalsBetween(9),
                yminorticksize=5,
                yminortickwidth=1.5,
                xlabelsize=25,
                ylabelsize=25,
                xticksize=8,
                yticksize=8,
                xtickwidth=2,
                ytickwidth=2,
                xticklabelsize=22,
                yticklabelsize=22,
                xtickalign=1,
                ytickalign=1,
                xaxisposition=:top,
                yaxisposition=:right,
            )

            lines!(axInset, Leff[ok], kpeak[ok], color=:black, linewidth=2)
            scatter!(axInset, Leff[ok], kpeak[ok]; color=cpeak[ok], markersize=14)
            xlims!(axInset, 0, 40)
            ymin = max(minimum(kpeak[ok]) * 0.9, eps(Float64))
            ymax = maximum(kpeak[ok]) * 1.1
            ylims!(axInset, ymin, ymax)
            axInset.xtickformat = (vs) -> [LaTeXString("\\mathrm{$(round(v; digits=1))}") for v in vs]
            set_pow10_ticks!(axInset; which=:y)
            axInset.backgroundcolor = (:white, 0.85)
        end

        display(fig)
    end
end

function _build_filtered_dict_from_slice(slice0::AbstractMatrix, cfg::InstrumentalConfig, scales)
    out = Dict{Float64, Matrix}()
    for Llarge in cfg.Llarge_list
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

function run_pipeline(cfg::InstrumentalConfig; flags::RunFlags=RunFlags())
    data = run_filter_pass(cfg)

    if flags.run_pmax_maps
        _plot_pmax_maps(cfg, data)
    end

    if flags.run_psd
        _plot_psd_panels(cfg, data)
        _plot_pmax_kx(cfg, data)
    end

    if flags.run_q_u_p_q2
        Qslice0 = get_chan_xy(data.Qdata, cfg.ichan, cfg.n, cfg.m)
        Uslice0 = get_chan_xy(data.Udata, cfg.ichan, cfg.n, cfg.m)

        q_filt = Dict{Float64, Matrix}(k => Matrix{Float64}(data.Qslice_filt[k]) for k in keys(data.Qslice_filt))
        u_filt = Dict{Float64, Matrix}(k => Matrix{Float64}(data.Uslice_filt[k]) for k in keys(data.Uslice_filt))

        _plot_component_spectrum(cfg, data.axes.kx, Qslice0, q_filt,
            LaTeXString("\\langle |\\tilde{Q}(k_x, y;\\nu_{50})|^2 \\rangle_y"); add_verticals=true)

        _plot_component_spectrum(cfg, data.axes.kx, Uslice0, u_filt,
            LaTeXString("\\langle |\\tilde{U}(k_x, y;\\nu_{50})|^2 \\rangle_y"); add_verticals=true)

        p_no = sqrt.(Qslice0.^2 .+ Uslice0.^2)
        p_filt = Dict{Float64, Matrix}()
        for Llarge in sort(data.L_ok)
            p_filt[Llarge] = sqrt.(Float64.(data.Qslice_filt[Llarge]).^2 .+ Float64.(data.Uslice_filt[Llarge]).^2)
        end
        _plot_component_spectrum(cfg, data.axes.kx, p_no, p_filt,
            LaTeXString("\\langle |\\tilde{P}(k_x, y;\\nu_{50})|^2 \\rangle_y");
            add_verticals=false,
            peak_window=(0.12, Inf),
            inset_title=:Pnu)

        q2_no = Qslice0 .^ 2
        q2_filt = Dict{Float64, Matrix}(k => Float64.(data.Qslice_filt[k]).^2 for k in keys(data.Qslice_filt))
        _plot_component_spectrum(cfg, data.axes.kx, q2_no, q2_filt,
            LaTeXString("\\left\\langle \\left|\\tilde{Q^2}(k_x,y;\\nu_{50})\\right|^2 \\right\\rangle_y"); add_verticals=true)
    end

    if flags.run_phi_q_u_p
        Qphi_cube = read_FITS(cfg.Q_in_phi)
        Uphi_cube = read_FITS(cfg.U_in_phi)

        ϕval = data.PhiArray[cfg.iphi]
        @info "Phi channel" cfg.iphi ϕval

        Qφ0 = get_chan_xy(Qphi_cube, cfg.iphi, cfg.n, cfg.m)
        Uφ0 = get_chan_xy(Uphi_cube, cfg.iphi, cfg.n, cfg.m)
        Pφ0 = sqrt.(Qφ0.^2 .+ Uφ0.^2)

        Qφf = _build_filtered_dict_from_slice(Qφ0, cfg, data.scales)
        Uφf = _build_filtered_dict_from_slice(Uφ0, cfg, data.scales)
        Pφf = Dict{Float64, Matrix}(k => sqrt.(Qφf[k].^2 .+ Uφf[k].^2) for k in keys(Qφf))

        _plot_component_spectrum(cfg, data.axes.kx, Qφ0, Qφf,
            LaTeXString("\\langle |\\tilde{Q}_{\\phi}(k_x,y)|^2 \\rangle_y"); add_verticals=true)

        _plot_component_spectrum(cfg, data.axes.kx, Uφ0, Uφf,
            LaTeXString("\\langle |\\tilde{U}_{\\phi}(k_x,y)|^2 \\rangle_y"); add_verticals=true)

        _plot_component_spectrum(cfg, data.axes.kx, Pφ0, Pφf,
            LaTeXString("\\langle |\\tilde{P}_{\\phi}(k_x,y)|^2 \\rangle_y");
            add_verticals=false,
            peak_window=(0.12, Inf),
            inset_title=:Pphi)
    end

    if flags.run_lic
        @info "LIC section requested, but default pipeline keeps it disabled as in original script comments."
    end

    return data
end
