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

    Pmaxf = RMSynthesis_pmax_map(Qf, Uf, nuArray, PhiArray; log_progress=log_progress)

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

    @debug "Grid geometry" nx=cfg.n ny=cfg.m dx_pc=scales.Δx dy_pc=scales.Δy f_nyquist_rad_pc_inv=scales.fNy

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
