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

    n::Int = 256
    m::Int = 256
    Lbox_pc::Float64 = 50.0
    Lcut_small::Float64 = 1.0
    Llarge_list::Vector{Float64} = [100.0, 50.0, 40.0, 20.0]

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
end

function load_config_defaults()
    root = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024"
    base = joinpath(root, "y", "Synchrotron", "WithFaraday")

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
    )
end

function rm_synthesis_axes(cfg::InstrumentalConfig)
    nuArray = (cfg.νmin_MHz:cfg.Δν_MHz:cfg.νmax_MHz) .* 1e6
    return nuArray, cfg.PhiArray
end
