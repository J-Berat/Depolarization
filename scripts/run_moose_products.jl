using TOML
using JSON
using Printf

function _arg_value(args::Vector{String}, key::String; default=nothing)
    idx = findfirst(==(key), args)
    idx === nothing && return default
    idx == length(args) && error("Missing value after $key")
    return args[idx + 1]
end

function _norm_yesno(s::AbstractString)
    v = uppercase(strip(String(s)))
    v in ("Y", "N") || error("Expected Y or N, got $s")
    return v
end

function _load_settings(cfg_path::AbstractString)
    cfg = TOML.parsefile(cfg_path)
    sim_root = String(cfg["paths"]["simulations_root"])
    moose_root = String(cfg["paths"]["moose_root"])
    emissivity = String(cfg["paths"]["emissivity_table"])
    sim_name = String(cfg["simulation"]["name"])
    los = String(cfg["simulation"]["los"])
    phi_min = Float64(cfg["tasks"]["instrumental_effect"]["phi_min"])
    phi_max = Float64(cfg["tasks"]["instrumental_effect"]["phi_max"])
    dphi = Float64(cfg["tasks"]["instrumental_effect"]["dphi"])
    zeta = Float64(cfg["tasks"]["ne_dm_em"]["zeta"])
    Geff = Float64(cfg["tasks"]["ne_dm_em"]["Geff"])
    phiPAH = Float64(cfg["tasks"]["ne_dm_em"]["phiPAH"])
    XC = Float64(cfg["tasks"]["ne_dm_em"]["XC"])
    return (
        sim_root=sim_root,
        moose_root=moose_root,
        emissivity=emissivity,
        sim_name=sim_name,
        los=los,
        phi_min=phi_min,
        phi_max=phi_max,
        dphi=dphi,
        zeta=zeta,
        Geff=Geff,
        phiPAH=phiPAH,
        XC=XC,
    )
end

function _moose_config(settings; faraday::String, los::String)
    sim_path = joinpath(settings.sim_root, settings.sim_name)
    return Dict(
        "base_dir" => settings.sim_root,
        "simulations" => [sim_path],
        "chosen_LOS" => [los],
        "conversionB" => 1000.0,
        "conversionn" => 1.0,
        "conversionT" => 1.0,
        "BoxLength_pc" => 50.0,
        "BoxLength_pix" => 256,
        "nustart" => 120.0,
        "nuend" => 167.0,
        "dnu" => 0.098,
        "FaradayRotation" => faraday,
        "phimin" => settings.phi_min,
        "phimax" => settings.phi_max,
        "dphi" => settings.dphi,
        "responseSynchrotron" => "N",
        "add_noise" => "N",
        "interpolation_file_path" => settings.emissivity,
        "ne_option" => "1",
        "zeta" => settings.zeta,
        "Geff" => settings.Geff,
        "phiPAH" => settings.phiPAH,
        "XC" => settings.XC,
        "log_progress" => true,
    )
end

function _run_moose(settings, cfg::Dict{String, Any}, tag::String)
    mkpath(joinpath(@__DIR__, "..", "outputs", "moose_configs"))
    cfg_path = joinpath(@__DIR__, "..", "outputs", "moose_configs", "moose_$(settings.sim_name)_LOS$(only(cfg["chosen_LOS"]))_$(tag).json")
    open(cfg_path, "w") do io
        write(io, JSON.json(cfg))
    end
    runner = joinpath(@__DIR__, "moose_from_config_runner.jl")
    cmd = `$(Base.julia_cmd()) --project=$(settings.moose_root) $runner $cfg_path`
    @info "Running MOOSE" tag cfg_path
    run(cmd)
    return cfg_path
end

function _make_nofaraday_pmax(settings; los::String)
    sim_path = joinpath(settings.sim_root, settings.sim_name)
    outdir = joinpath(sim_path, los, "Synchrotron", "noFaraday")
    q_path = joinpath(outdir, "Qnu.fits")
    u_path = joinpath(outdir, "Unu.fits")
    helper = joinpath(@__DIR__, "moose_make_pmax_nofaraday.jl")
    cmd = `$(Base.julia_cmd()) --project=$(settings.moose_root) $helper --q $q_path --u $u_path --outdir $outdir --nustart 120.0 --nuend 167.0 --dnu 0.098 --phimin $(settings.phi_min) --phimax $(settings.phi_max) --dphi $(settings.dphi)`
    @info "Computing noFaraday Pmax with MOOSE RM synthesis" outdir
    run(cmd)
end

function main(args::Vector{String})
    cfg_path = something(_arg_value(args, "--config"; default=joinpath(@__DIR__, "..", "config", "default.toml")))
    settings = _load_settings(cfg_path)
    los = String(something(_arg_value(args, "--los"; default=settings.los)))
    mode = lowercase(String(something(_arg_value(args, "--mode"; default="both"))))
    mode in ("both", "withfaraday", "nofaraday") || error("mode must be one of both, withfaraday, nofaraday")

    if mode in ("both", "withfaraday")
        cfgY = _moose_config(settings; faraday="Y", los=los)
        _run_moose(settings, cfgY, "withfaraday")
    end
    if mode in ("both", "nofaraday")
        cfgN = _moose_config(settings; faraday="N", los=los)
        _run_moose(settings, cfgN, "nofaraday")
        _make_nofaraday_pmax(settings; los=los)
    end

    sim_path = joinpath(settings.sim_root, settings.sim_name, los, "Synchrotron")
    @info "MOOSE products ready" los mode path=sim_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
