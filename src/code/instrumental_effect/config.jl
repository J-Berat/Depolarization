Base.@kwdef struct InstrumentalConfig
    Q_in::String
    U_in::String
    Pmax_nofilter_path::String
    Q_in_phi::String
    U_in_phi::String
    base_out::String

    Bx_in::String
    By_in::String
    Bz_in::String
    dens_in::String
    los::String = "y"

    n::Int = 256
    m::Int = 256
    Lbox_pc::Float64 = 50.0
    Lcut_small::Float64 = 1.0
    Llarge_list::Vector{Float64} = [200.0, 180.0, 165.0, 150.0, 135.0, 120.0, 105.0, 90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0]

    νmin_MHz::Float64 = 120.0
    νmax_MHz::Float64 = 167.0
    Δν_MHz::Float64 = 0.098
    PhiArray::Vector{Float64} = collect(-10.0:0.25:10.0)

    ichan::Int = 50
    iphi::Int = 30
end

Base.@kwdef struct RunFlags
    run_pmax_maps::Bool = true
    run_psd::Bool = true
    run_q_u_p_q2::Bool = true
    run_phi_q_u_p::Bool = true
    run_lic::Bool = false
    run_channel_b_alignment::Bool = true
end

# Backward-compatible positional constructor kept for existing call sites/tests.
RunFlags(run_pmax_maps::Bool, run_psd::Bool, run_q_u_p_q2::Bool, run_phi_q_u_p::Bool, run_lic::Bool) =
    RunFlags(
        run_pmax_maps=run_pmax_maps,
        run_psd=run_psd,
        run_q_u_p_q2=run_q_u_p_q2,
        run_phi_q_u_p=run_phi_q_u_p,
        run_lic=run_lic,
        run_channel_b_alignment=false,
    )

"""
    load_config_defaults(...)

    Builds default instrumental-pipeline config.
"""
function load_config_defaults()
    sim_root = get(ENV, "SIMULATIONS_ROOT", "./data/simu_RAMSES")
    simu_name = get(ENV, "SIMU_NAME", "d1cf05bx10rms18000nograv1024")
    los = get(ENV, "SIMU_LOS", "y")

    root = joinpath(sim_root, simu_name)
    base = joinpath(root, los, "Synchrotron", "WithFaraday")

    InstrumentalConfig(
        Q_in = joinpath(base, "Qnu.fits"),
        U_in = joinpath(base, "Unu.fits"),
        Pmax_nofilter_path = joinpath(base, "Pmax.fits"),
        Q_in_phi = joinpath(base, "realFDF.fits"),
        U_in_phi = joinpath(base, "imagFDF.fits"),
        base_out = joinpath(base, "WithInstrument"),
        Bx_in = joinpath(root, "Bx.fits"),
        By_in = joinpath(root, "By.fits"),
        Bz_in = joinpath(root, "Bz.fits"),
        dens_in = joinpath(root, "density.fits"),
        los = los,
    )
end

"""
    rm_synthesis_axes(...)

    Returns `nu` and `phi` axes for RM synthesis.
"""
function rm_synthesis_axes(cfg::InstrumentalConfig)
    nuArray = (cfg.νmin_MHz:cfg.Δν_MHz:cfg.νmax_MHz) .* 1e6
    return nuArray, cfg.PhiArray
end
