"""
    _filter_tag(...)

    Builds filter tag text.
"""
function _filter_tag(Llarge::Real)
    @sprintf("HardBandPass_remove_L0_to_1pc_and_%dto50pc", Int(round(Llarge)))
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

"""
    _integrity_push_check!(...)

Adds one integrity check record to the report lists.
"""
function _integrity_push_check!(checks::Vector{Dict{String,Any}}, critical_failures::Vector{String},
                                warnings::Vector{String}, name::String, passed::Bool;
                                critical::Bool=true, details::AbstractString="")
    level = critical ? "critical" : "warning"
    msg = passed ? "ok" : (isempty(details) ? "failed" : details)
    push!(checks, Dict{String,Any}(
        "name" => name,
        "level" => level,
        "passed" => passed,
        "message" => msg,
    ))
    if !passed
        if critical
            push!(critical_failures, name)
        else
            push!(warnings, name)
        end
    end
    return nothing
end

"""
    _write_integrity_report(...)

Writes a compact integrity report next to run outputs.
"""
function _write_integrity_report(path::AbstractString, report::Dict{String,Any})
    lines = String[
        "integrity_status=$(report["status"])",
        "critical_failures=$(report["critical_failures"])",
        "warnings=$(report["warnings"])",
        "checks_total=$(report["checks_total"])",
        "",
        "[checks]",
    ]

    for item in report["checks"]
        push!(lines, string(
            "- ", item["name"],
            " | level=", item["level"],
            " | passed=", item["passed"],
            " | message=", item["message"],
        ))
    end

    open(path, "w") do io
        write(io, join(lines, "\n"))
        write(io, "\n")
    end
    return path
end

"""
    _validate_filter_inputs(...)

Validates dimensions, NaN/Inf, units consistency and expected physical ranges.
Critical failures abort the run.
"""
function _validate_filter_inputs(cfg::InstrumentalConfig, Qdata, Udata, Pmax0, nuArray, PhiArray)
    checks = Dict{String,Any}[]
    critical_failures = String[]
    warnings = String[]

    try
        validate_instrumental_config!(cfg)
        _integrity_push_check!(checks, critical_failures, warnings, "config.basics", true; critical=true)
    catch err
        _integrity_push_check!(checks, critical_failures, warnings, "config.basics", false;
                               critical=true, details=sprint(showerror, err))
    end

    nν = length(nuArray)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.nu.non_empty", nν > 0; critical=true)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.phi.non_empty", length(PhiArray) > 0; critical=true)

    nu_finite = !isempty(nuArray) && all(isfinite, nuArray)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.nu.finite", nu_finite; critical=true)
    nu_positive = !isempty(nuArray) && all(>(0), nuArray)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.nu.positive_hz", nu_positive; critical=true)
    nu_sorted = length(nuArray) <= 1 || all(diff(nuArray) .> 0)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.nu.strictly_increasing", nu_sorted; critical=true)

    phi_finite = !isempty(PhiArray) && all(isfinite, PhiArray)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.phi.finite", phi_finite; critical=true)
    phi_sorted = length(PhiArray) <= 1 || all(diff(PhiArray) .> 0)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.phi.strictly_increasing", phi_sorted; critical=true)

    if !isempty(nuArray)
        expected_min_hz = cfg.νmin_MHz * 1e6
        expected_max_hz = cfg.νmax_MHz * 1e6
        expected_step_hz = cfg.Δν_MHz * 1e6

        first_ok = isapprox(first(nuArray), expected_min_hz; rtol=1e-10, atol=1e-6)
        step_ok = length(nuArray) <= 1 || all(isapprox.(diff(nuArray), expected_step_hz; rtol=1e-10, atol=1e-6))
        # Range syntax (νmin:Δν:νmax) may stop before νmax when Δν does not divide the span.
        last_not_above_max = last(nuArray) <= expected_max_hz + 1e-6

        nν = length(nuArray)
        expected_last_from_step = expected_min_hz + (nν - 1) * expected_step_hz
        last_consistent_with_step = isapprox(last(nuArray), expected_last_from_step; rtol=1e-10, atol=1e-6)

        units_ok = first_ok && step_ok && last_not_above_max && last_consistent_with_step
        units_details = units_ok ? "" : string(
            "nu-axis mismatch (first=", first(nuArray),
            ", last=", last(nuArray),
            ", step_expected=", expected_step_hz,
            ", max_cfg=", expected_max_hz, ")"
        )
        _integrity_push_check!(checks, critical_failures, warnings, "units.frequency.MHz_to_Hz", units_ok;
                               critical=true, details=units_details)
    else
        _integrity_push_check!(checks, critical_failures, warnings, "units.frequency.MHz_to_Hz", false;
                               critical=true, details="empty nuArray")
    end

    try
        q_layout = require_stokes_cube_layout(Qdata, "Qnu", cfg, nν)
        u_layout = require_stokes_cube_layout(Udata, "Unu", cfg, nν)
        layouts_match = q_layout == u_layout
        _integrity_push_check!(checks, critical_failures, warnings, "dims.stokes_layout", layouts_match;
                               critical=true, details=layouts_match ? "" : "Qnu/Unu layout mismatch")
    catch err
        _integrity_push_check!(checks, critical_failures, warnings, "dims.stokes_layout", false;
                               critical=true, details=sprint(showerror, err))
    end

    try
        require_channel_index(cfg.ichan, nν, "ichan")
        _integrity_push_check!(checks, critical_failures, warnings, "index.ichan.in_bounds", true; critical=true)
    catch err
        _integrity_push_check!(checks, critical_failures, warnings, "index.ichan.in_bounds", false;
                               critical=true, details=sprint(showerror, err))
    end

    try
        require_map_layout(Pmax0, "Pmax", cfg)
        _integrity_push_check!(checks, critical_failures, warnings, "dims.pmax_layout", true; critical=true)
    catch err
        _integrity_push_check!(checks, critical_failures, warnings, "dims.pmax_layout", false;
                               critical=true, details=sprint(showerror, err))
    end

    q_finite = all(isfinite, Qdata)
    u_finite = all(isfinite, Udata)
    pmax_finite = all(isfinite, Pmax0)
    _integrity_push_check!(checks, critical_failures, warnings, "nan_inf.Qnu", q_finite; critical=true)
    _integrity_push_check!(checks, critical_failures, warnings, "nan_inf.Unu", u_finite; critical=true)
    _integrity_push_check!(checks, critical_failures, warnings, "nan_inf.Pmax", pmax_finite; critical=true)

    pmax_non_negative = all(Pmax0 .>= 0)
    _integrity_push_check!(checks, critical_failures, warnings, "physical.Pmax.non_negative", pmax_non_negative;
                           critical=true, details="Pmax must be >= 0")

    freq_band_plausible = (1.0 <= cfg.νmin_MHz <= cfg.νmax_MHz <= 5.0e4)
    _integrity_push_check!(checks, critical_failures, warnings, "physical.frequency_band_plausible", freq_band_plausible;
                           critical=false, details="Expected ~[1, 50000] MHz")

    phi_span_plausible = !isempty(PhiArray) && (maximum(abs.(PhiArray)) <= 1.0e5)
    _integrity_push_check!(checks, critical_failures, warnings, "physical.phi_range_plausible", phi_span_plausible;
                           critical=false, details="Expected |phi| <= 1e5 rad m^-2")

    report = Dict{String,Any}(
        "status" => isempty(critical_failures) ? "pass" : "fail",
        "critical_failures" => length(critical_failures),
        "warnings" => length(warnings),
        "checks_total" => length(checks),
        "critical_failed_checks" => copy(critical_failures),
        "warning_checks" => copy(warnings),
        "checks" => checks,
    )

    report_path = _write_integrity_report(joinpath(cfg.base_out, "integrity_report.txt"), report)
    report["report_path"] = report_path

    if isempty(critical_failures)
        @info "Integrity checks passed" checks=length(checks) warnings=length(warnings) report=report_path
    else
        @error "Integrity critical checks failed; aborting run" failed=critical_failures report=report_path
        error("Critical integrity checks failed: $(join(critical_failures, ", "))")
    end

    return report
end

"""
    _process_single_filter(...)

    Computes one filter level and writes filtered outputs.
"""
function _process_single_filter(cfg::InstrumentalConfig, scales, nuArray, PhiArray,
                                Qdata, Udata, Llarge::Float64; log_progress::Bool)
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

    Qslice = Float32.(get_chan_xy(Qf, cfg.ichan, cfg.n, cfg.m))
    Uslice = Float32.(get_chan_xy(Uf, cfg.ichan, cfg.n, cfg.m))

    write_FITS(joinpath(subdir, "Qnu_filtered.fits"), Qf)
    write_FITS(joinpath(subdir, "Unu_filtered.fits"), Uf)

    absF, _, _ = RMSynthesis(Qf, Uf, nuArray, PhiArray; log_progress=log_progress)
    Pmaxf = Float32.(Pphi_max_map(absF))

    rm_dir = joinpath(subdir, "RMSynthesis")
    mkpath(rm_dir)
    write_FITS(joinpath(rm_dir, "Pphi_max.fits"), Pmaxf)

    return (Llarge=Llarge, Qslice=Qslice, Uslice=Uslice, Pmaxf=Pmaxf)
end

"""
    run_filter_pass(...)

    Executes filtering pass and writes filtered products.
"""
function run_filter_pass(cfg::InstrumentalConfig)
    mkpath(cfg.base_out)

    scales = grid_scales(cfg)
    axes = spectral_axes(cfg, scales.Δx, scales.Δy)

    @debug "Grid" cfg.n cfg.m scales.Δx scales.Δy scales.fNy

    nuArray, PhiArray = rm_synthesis_axes(cfg)

    Qdata = read_FITS(cfg.Q_in)
    Udata = read_FITS(cfg.U_in)
    Pmax0 = Float32.(read_FITS(cfg.Pmax_nofilter_path))

    integrity = _validate_filter_inputs(cfg, Qdata, Udata, Pmax0, nuArray, PhiArray)

    Qslice_filt = Dict{Float64, Matrix{Float32}}()
    Uslice_filt = Dict{Float64, Matrix{Float32}}()
    Pmax_filt = Dict{Float64, Matrix{Float32}}()
    L_ok = Float64[]

    log_progress = _rm_log_progress_enabled()

    for Llarge in Float64.(cfg.Llarge_list)
        out = _process_single_filter(cfg, scales, nuArray, PhiArray, Qdata, Udata, Llarge; log_progress=log_progress)
        Qslice_filt[Llarge] = out.Qslice
        Uslice_filt[Llarge] = out.Uslice
        Pmax_filt[Llarge] = out.Pmaxf
        push!(L_ok, Llarge)
    end

    isempty(L_ok) && error("No filtered Pphi_max.fits products generated under: $(cfg.base_out)")

    return FilterPassResult(
        Qdata,
        Udata,
        Qslice_filt,
        Uslice_filt,
        Pmax0,
        Pmax_filt,
        L_ok,
        scales,
        axes,
        Vector{Float64}(nuArray),
        Vector{Float64}(PhiArray),
        integrity,
    )
end

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
            @info "Inset Pmax: k-window" kmin_win kmax_win

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
            scatter!(axInset, Leff[ok], kpeak[ok]; color=cpeak[ok], markersize=14)
            xlims!(axInset, 0, 40)

            ymin = max(minimum(kpeak[ok]) * 0.9, eps(Float64))
            ymax = maximum(kpeak[ok]) * 1.1
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
                                  inset_title::Union{Nothing, Symbol}=nothing)
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
                @info "Inset peak" inset_title Llarge kp
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
            scatter!(axInset, Leff[ok], kpeak[ok]; color=cpeak[ok], markersize=14)
            xlims!(axInset, 0, 40)
            ymin = max(minimum(kpeak[ok]) * 0.9, eps(Float64))
            ymax = maximum(kpeak[ok]) * 1.1
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

"""
    run_pipeline(...)

    Orchestrates full instrumental-effect pipeline.
"""
function run_pipeline(cfg::InstrumentalConfig; flags::RunFlags=RunFlags(), on_step::Function=(_ -> nothing))
    on_step(:filter_pass)
    data = run_filter_pass(cfg)

    if flags.run_pmax_maps
        on_step(:run_pmax_maps)
        _plot_pmax_maps(cfg, data)
    end

    if flags.run_psd
        on_step(:run_psd)
        _plot_psd_panels(cfg, data)
        _plot_pmax_kx(cfg, data)
    end

    if flags.run_q_u_p_q2
        on_step(:run_q_u_p_q2)
        Qslice0 = get_chan_xy(data.Qdata, cfg.ichan, cfg.n, cfg.m)
        Uslice0 = get_chan_xy(data.Udata, cfg.ichan, cfg.n, cfg.m)

        q_filt = copy(data.Qslice_filt)
        u_filt = copy(data.Uslice_filt)

        _plot_component_spectrum(cfg, data.axes.kx, Qslice0, q_filt,
            LaTeXString("S_{Q}(k_x;\\nu_{50})"); add_verticals=true)

        _plot_component_spectrum(cfg, data.axes.kx, Uslice0, u_filt,
            LaTeXString("S_{U}(k_x;\\nu_{50})"); add_verticals=true)

        p_no = sqrt.(Qslice0.^2 .+ Uslice0.^2)
        p_filt = Dict{Float64, Matrix{Float64}}()
        for Llarge in sort(data.L_ok)
            qL = data.Qslice_filt[Llarge]
            uL = data.Uslice_filt[Llarge]
            p_filt[Llarge] = Float64.(sqrt.(qL.^2 .+ uL.^2))
        end
        _plot_component_spectrum(cfg, data.axes.kx, p_no, p_filt,
            LaTeXString("S_{P}(k_x;\\nu_{50})");
            add_verticals=false,
            peak_window=(0.12, Inf),
            inset_title=:Pnu)

        q2_no = Qslice0 .^ 2
        q2_filt = Dict{Float64, Matrix{Float64}}(k => Float64.(data.Qslice_filt[k].^2) for k in keys(data.Qslice_filt))
        _plot_component_spectrum(cfg, data.axes.kx, q2_no, q2_filt,
            LaTeXString("S_{Q^{2}}(k_x;\\nu_{50})"); add_verticals=true)
    end

    if flags.run_phi_q_u_p
        on_step(:run_phi_q_u_p)
        Qphi_cube = read_FITS(cfg.Q_in_phi)
        Uphi_cube = read_FITS(cfg.U_in_phi)
        _validate_phi_cubes(cfg, data, Qphi_cube, Uphi_cube)

        ϕval = data.PhiArray[cfg.iphi]
        @info "Phi channel" cfg.iphi ϕval

        Qφ0 = get_chan_xy(Qphi_cube, cfg.iphi, cfg.n, cfg.m)
        Uφ0 = get_chan_xy(Uphi_cube, cfg.iphi, cfg.n, cfg.m)
        Pφ0 = sqrt.(Qφ0.^2 .+ Uφ0.^2)

        Qφf = _build_filtered_dict_from_slice(Qφ0, cfg, data.scales)
        Uφf = _build_filtered_dict_from_slice(Uφ0, cfg, data.scales)
        Pφf = Dict{Float64, Matrix{Float64}}(k => sqrt.(Qφf[k].^2 .+ Uφf[k].^2) for k in keys(Qφf))

        _plot_component_spectrum(cfg, data.axes.kx, Qφ0, Qφf,
            LaTeXString("S_{Q_{\\phi}}(k_x)"); add_verticals=true)

        _plot_component_spectrum(cfg, data.axes.kx, Uφ0, Uφf,
            LaTeXString("S_{U_{\\phi}}(k_x)"); add_verticals=true)

        _plot_component_spectrum(cfg, data.axes.kx, Pφ0, Pφf,
            LaTeXString("S_{P_{\\phi}}(k_x)");
            add_verticals=false,
            peak_window=(0.12, Inf),
            inset_title=:Pphi)
    end

    if flags.run_lic
        on_step(:run_lic)
        @debug "LIC section requested, but default pipeline keeps it disabled as in original script comments."
    end

    return data
end
