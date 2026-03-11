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

        q_filt = data.Qslice_filt
        u_filt = data.Uslice_filt

        _plot_component_spectrum(cfg, data.axes.kx, Qslice0, q_filt,
            LaTeXString("S_{Q}(k_x;\\nu_{50})"); add_verticals=true, save_tag="q_nu50")

        _plot_component_spectrum(cfg, data.axes.kx, Uslice0, u_filt,
            LaTeXString("S_{U}(k_x;\\nu_{50})"); add_verticals=true, save_tag="u_nu50")

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
            inset_title=:Pnu,
            save_tag="p_nu50")

        q2_no = Qslice0 .^ 2
        q2_filt = Dict{Float64, Matrix{Float64}}(k => Float64.(data.Qslice_filt[k].^2) for k in keys(data.Qslice_filt))
        _plot_component_spectrum(cfg, data.axes.kx, q2_no, q2_filt,
            LaTeXString("S_{Q^{2}}(k_x;\\nu_{50})"); add_verticals=true, save_tag="q2_nu50")
    end

    if flags.run_phi_q_u_p
        on_step(:run_phi_q_u_p)
        Qphi_cube = read_FITS(cfg.Q_in_phi)
        Uphi_cube = read_FITS(cfg.U_in_phi)
        _validate_phi_cubes(cfg, data, Qphi_cube, Uphi_cube)

        ϕval = data.PhiArray[cfg.iphi]
        @info "Phi channel selected" iphi=cfg.iphi phi_rad_m2=ϕval

        Qφ0 = get_chan_xy(Qphi_cube, cfg.iphi, cfg.n, cfg.m)
        Uφ0 = get_chan_xy(Uphi_cube, cfg.iphi, cfg.n, cfg.m)
        Pφ0 = sqrt.(Qφ0.^2 .+ Uφ0.^2)

        Qφf = _build_filtered_dict_from_slice(Qφ0, cfg, data.scales)
        Uφf = _build_filtered_dict_from_slice(Uφ0, cfg, data.scales)
        Pφf = Dict{Float64, Matrix{Float64}}(k => sqrt.(Qφf[k].^2 .+ Uφf[k].^2) for k in keys(Qφf))

        _plot_component_spectrum(cfg, data.axes.kx, Qφ0, Qφf,
            LaTeXString("S_{Q_{\\phi}}(k_x)"); add_verticals=true, save_tag="q_phi")

        _plot_component_spectrum(cfg, data.axes.kx, Uφ0, Uφf,
            LaTeXString("S_{U_{\\phi}}(k_x)"); add_verticals=true, save_tag="u_phi")

        _plot_component_spectrum(cfg, data.axes.kx, Pφ0, Pφf,
            LaTeXString("S_{P_{\\phi}}(k_x)");
            add_verticals=false,
            peak_window=(0.12, Inf),
            inset_title=:Pphi,
            save_tag="p_phi")
    end

    if flags.run_channel_b_alignment
        on_step(:run_channel_b_alignment)
        _run_channel_b_alignment(cfg, data)
    end

    if flags.run_lic
        on_step(:run_lic)
        @debug "LIC section requested, but default pipeline keeps it disabled as in original script comments."
    end

    return data
end
