include(joinpath(@__DIR__, "..", "lib", "DepolLib.jl"))
include(joinpath(@__DIR__, "..", "instrumental_effect", "InstrumentalEffect.jl"))

using .DepolLib
using .InstrumentalEffect

function run_instrumental_effect_job(cfg)::Dict{String,Any}
    sim_root = resolve_simulations_root(cfg)
    sim = string(cfg_require(cfg, ["simulation", "name"]))
    los = require_los(string(cfg_require(cfg, ["simulation", "los"])))

    q_in = withfaraday_path(sim_root, sim, los, "Qnu.fits")
    u_in = withfaraday_path(sim_root, sim, los, "Unu.fits")
    pmax_in = withfaraday_path(sim_root, sim, los, "Pmax.fits")
    q_phi = withfaraday_path(sim_root, sim, los, "realFDF.fits")
    u_phi = withfaraday_path(sim_root, sim, los, "imagFDF.fits")
    bx_in = simulation_field_path(sim_root, sim, "Bx")
    by_in = simulation_field_path(sim_root, sim, "By")
    bz_in = simulation_field_path(sim_root, sim, "Bz")
    dens_in = simulation_field_path(sim_root, sim, "density")
    llarge_list = Float64.(cfg_get(cfg, ["tasks", "instrumental_effect", "llarge_list"];
                                   default=[200.0, 180.0, 165.0, 150.0, 135.0, 120.0, 105.0, 90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0]))

    try
        require_existing_files([
            q_in, u_in, pmax_in, q_phi, u_phi, bx_in, by_in, bz_in, dens_in
        ]; context="instrumental_effect")
    catch err
        emsg = sprint(showerror, err)
        configured_root = string(cfg_get(cfg, ["paths", "simulations_root"]; default="(missing)"))
        error(
            string(
                emsg, "\n",
                "Hint: check paths.simulations_root in your config.\n",
                "Configured value: ", configured_root, "\n",
                "Example override:\n",
                "  julia --startup-file=no src/code/jobs/run_instrumental_effect_job.jl ",
                "--config config/default.toml ",
                "--set paths.simulations_root=/path/to/simu_RAMSES\n"
            )
        )
    end

    run_cfg = InstrumentalConfig(
        Q_in = q_in,
        U_in = u_in,
        Pmax_nofilter_path = pmax_in,
        Q_in_phi = q_phi,
        U_in_phi = u_phi,
        base_out = standard_output_dir(cfg, "instrumental"; simu=sim, los=los),
        Bx_in = bx_in,
        By_in = by_in,
        Bz_in = bz_in,
        dens_in = dens_in,
        los = los,
        Llarge_list = llarge_list,
    )

    flags = RunFlags(
        run_pmax_maps = Bool(cfg_get(cfg, ["tasks", "instrumental_effect", "run_pmax_maps"]; default=true)),
        run_psd = Bool(cfg_get(cfg, ["tasks", "instrumental_effect", "run_psd"]; default=true)),
        run_q_u_p_q2 = Bool(cfg_get(cfg, ["tasks", "instrumental_effect", "run_q_u_p_q2"]; default=true)),
        run_phi_q_u_p = Bool(cfg_get(cfg, ["tasks", "instrumental_effect", "run_phi_q_u_p"]; default=true)),
        run_lic = Bool(cfg_get(cfg, ["tasks", "instrumental_effect", "run_lic"]; default=false)),
        run_channel_b_alignment = Bool(cfg_get(cfg, ["tasks", "instrumental_effect", "run_channel_b_alignment"]; default=true)),
    )

    data = run_pipeline(run_cfg; flags=flags)

    return Dict(
        "task" => "instrumental",
        "output_dir" => run_cfg.base_out,
        "integrity_report" => get(data.integrity, "report_path", nothing),
        "integrity_status" => get(data.integrity, "status", "unknown"),
        "integrity_critical_failures" => get(data.integrity, "critical_failures", -1),
        "integrity_warnings" => get(data.integrity, "warnings", -1),
        "flags" => Dict(
            "run_pmax_maps" => flags.run_pmax_maps,
            "run_psd" => flags.run_psd,
            "run_q_u_p_q2" => flags.run_q_u_p_q2,
            "run_phi_q_u_p" => flags.run_phi_q_u_p,
            "run_lic" => flags.run_lic,
            "run_channel_b_alignment" => flags.run_channel_b_alignment,
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_job_entrypoint("instrumental", run_instrumental_effect_job)
end
