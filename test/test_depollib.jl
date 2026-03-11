using Test

@testset "DepolLib config and helpers" begin
    @testset "runtime config parsing" begin
        mktempdir() do td
            cfg_path = joinpath(td, "cfg.toml")
            open(cfg_path, "w") do io
                println(io, "[paths]")
                println(io, "simulations_root = \"$td\"")
                println(io, "desktop_output_root = \"$td/out\"")
                println(io, "[simulation]")
                println(io, "name = \"simA\"")
                println(io, "los = \"y\"")
                println(io, "[tasks.demo]")
                println(io, "enabled = false")
            end

            cfg = load_runtime_config("unit_job"; args=[
                "--config", cfg_path,
                "--set", "simulation.los=z",
                "--set", "tasks.demo.enabled=true",
                "--set", "runtime_answer=42",
                "--set", "floaty=3.5",
            ])

            @test cfg["runtime"]["task"] == "unit_job"
            @test cfg_get(cfg, ["simulation", "los"]) == "z"
            @test cfg_get(cfg, ["tasks", "demo", "enabled"]) === true
            @test cfg_get(cfg, ["runtime_answer"]) == 42
            @test cfg_get(cfg, ["floaty"]) == 3.5

            @test cfg_require(cfg, ["simulation", "name"]) == "simA"
            @test_throws ErrorException cfg_require(cfg, ["missing", "key"])
        end
    end

    @testset "build_cfg overrides" begin
        mktempdir() do td
            cfg_path = joinpath(td, "cfg.toml")
            open(cfg_path, "w") do io
                println(io, "[paths]")
                println(io, "simulations_root = \"$td\"")
                println(io, "desktop_output_root = \"$td/out\"")
                println(io, "[simulation]")
                println(io, "name = \"simA\"")
                println(io, "los = \"y\"")
                println(io, "[tasks.demo]")
                println(io, "enabled = false")
            end

            cfg = build_cfg("unit_job";
                config_path=cfg_path,
                overrides=Dict{Any,Any}(
                    "simulation.los" => "z",
                    "tasks.demo.enabled" => true,
                    :runtime_answer => 42,
                    "floaty" => "3.5",
                ),
            )

            @test cfg["runtime"]["task"] == "unit_job"
            @test cfg["runtime"]["config_path"] == abspath(cfg_path)
            @test cfg_get(cfg, ["simulation", "los"]) == "z"
            @test cfg_get(cfg, ["tasks", "demo", "enabled"]) === true
            @test cfg_get(cfg, ["runtime_answer"]) == 42
            @test cfg_get(cfg, ["floaty"]) == 3.5

            @test length(cfg["runtime"]["set_args"]) == 4

            @test_throws ErrorException build_cfg("unit_job";
                config_path=cfg_path,
                overrides=Dict{Any,Any}(1 => "bad"))
        end
    end

    @testset "math and interval utilities" begin
        vals = [1.0, NaN, Inf, -1.0]
        @test finite_values(vals) == [1.0, -1.0]
        @test finite_minmax(vals) == (-1.0, 1.0)

        @test require_los("x") == "x"
        @test_throws ErrorException require_los("bad")

        merged = merge_intervals([(2, 4), (5, 7), (12, 13)]; gap_pix=1)
        @test merged == [(2, 7), (12, 13)]

        A = zeros(2, 2)
        B = ones(2, 2)
        @test require_ndims(A, 2, "A") === A
        @test require_same_size([A, B], ["A", "B"]) === true
        @test_throws ErrorException require_same_size([A, ones(3, 3)], ["A", "B"])
    end

    @testset "Wolfire ne and LOS DM/EM maps" begin
        T = fill(Float32(100.0), 2, 2, 2)
        n = fill(Float32(10.0), 2, 2, 2)

        ne = Wolfire_ne(2.5e-16, 1.0, 0.5, 1.4e-4, T, n)
        expected_ne = 2.4e-3 * sqrt(2.5e-16 / 1e-16) * (100.0 / 100.0)^0.25 * sqrt(1.0) / 0.5 + 10.0 * 1.4e-4
        @test all(isapprox.(ne, expected_ne; rtol=1e-6, atol=0.0))

        dm_x, em_x, dl_x = dm_em_maps(ne, "x"; lbox_pc=50.0)
        @test size(dm_x) == (2, 2)
        @test size(em_x) == (2, 2)
        @test dl_x == 25.0
        @test all(isapprox.(dm_x, 2 * expected_ne * 25.0; rtol=1e-6, atol=0.0))
        @test all(isapprox.(em_x, 2 * expected_ne^2 * 25.0; rtol=1e-6, atol=0.0))

        dm_z, em_z, dl_z = dm_em_maps(ne, "z"; lbox_pc=100.0)
        @test size(dm_z) == (2, 2)
        @test size(em_z) == (2, 2)
        @test dl_z == 50.0
        @test all(isapprox.(dm_z, 2 * expected_ne * 50.0; rtol=1e-6, atol=0.0))
        @test all(isapprox.(em_z, 2 * expected_ne^2 * 50.0; rtol=1e-6, atol=0.0))
    end

    @testset "path + orchestration helpers" begin
        root = "/tmp/simroot"
        @test simulation_dir(root, "sim01") == "/tmp/simroot/sim01"
        @test simulation_field_path(root, "sim01", "Bx") == "/tmp/simroot/sim01/Bx.fits"
        @test synchrotron_dir(root, "sim01", "y") == "/tmp/simroot/sim01/y/Synchrotron"
        @test withfaraday_path(root, "sim01", "y", "Pmax.fits") == "/tmp/simroot/sim01/y/Synchrotron/WithFaraday/Pmax.fits"

        cfg = Dict(
            "tasks" => Dict("demo" => Dict("enabled" => false)),
            "simulation" => Dict("los" => "x"),
        )
        @test task_enabled(cfg, ["tasks", "demo", "enabled"]; default=true) == false
        @test task_enabled(cfg, ["tasks", "missing", "enabled"]; default=true) == true

        skipped = skipped_job_result("demo_job", "disabled")
        @test skipped["task"] == "demo_job"
        @test skipped["status"] == "skipped"
        @test skipped["reason"] == "disabled"

        mktempdir() do td
            cfg_path = joinpath(td, "cfg.toml")
            open(cfg_path, "w") do io
                println(io, "[paths]")
                println(io, "simulations_root = \"$td\"")
                println(io, "desktop_output_root = \"$td/out\"")
                println(io, "[simulation]")
                println(io, "name = \"simA\"")
                println(io, "los = \"x\"")
            end

            job_result = Dict{String,Any}()
            redirect_stdout(devnull) do
                job_result = run_job_entrypoint("demo_task", cfg -> Dict(
                    "task" => cfg["runtime"]["task"],
                    "los" => cfg["simulation"]["los"],
                ); args=["--config", cfg_path, "--set", "simulation.los=z"])
            end

            @test job_result["task"] == "demo_task"
            @test job_result["los"] == "z"
        end

        mktempdir() do td
            preferred = joinpath(td, "preferred")
            legacy = joinpath(td, "legacy")
            cfg_preferred = Dict(
                "paths" => Dict(
                    "desktop_output_root" => preferred,
                    "outputs_root" => legacy,
                ),
                "simulation" => Dict("name" => "simA", "los" => "y"),
            )
            dir_preferred = standard_output_dir(cfg_preferred, "demo_job")
            @test startswith(dir_preferred, preferred)

            cfg_legacy = Dict(
                "paths" => Dict("outputs_root" => legacy),
                "simulation" => Dict("name" => "simA", "los" => "y"),
            )
            dir_legacy = standard_output_dir(cfg_legacy, "demo_job")
            @test startswith(dir_legacy, legacy)
        end

        missing_path = "/tmp/does-not-exist-file.fits"
        err = try
            require_existing_files([missing_path]; context="unit-check")
            nothing
        catch ex
            ex
        end
        @test err isa ErrorException
        @test occursin("unit-check", sprint(showerror, err))
        @test occursin(missing_path, sprint(showerror, err))
    end

    @testset "repl help" begin
        help_text = sprint(io -> depol_help(io))
        alias_text = sprint(io -> repl_help(io))
        jobs_text = sprint(io -> depol_help(io; topic="jobs"))
        one_job_text = sprint(io -> depol_help(io; topic="instrumental"))
        canal_alias_text = sprint(io -> depol_help(io; topic="canal"))
        cli_text = sprint(io -> depol_help(io; topic="cli"))
        order_text = sprint(io -> depol_help(io; topic="order"))

        @test occursin("Depolarization REPL quick help", help_text)
        @test occursin("Jobs: what they do + how to run in REPL", help_text)
        @test occursin("Recommended order (fresh simulation)", help_text)
        @test occursin("CLI cheatsheet (one command per job)", help_text)
        @test occursin("run_ne_dm_em_job(cfg)", help_text)
        @test occursin("run_mach_suite(cfg)", help_text)
        @test occursin("run_reversal_transition_job(cfg)", help_text)
        @test occursin("run_canal_metrics_job(cfg)", help_text)
        @test occursin("run_segmentation_pipeline_job(cfg)", help_text)
        @test occursin("run_instrumental_effect_job(cfg)", help_text)
        @test occursin("depol_help(\"canal\")", help_text)
        @test occursin("depol_help(\"instrumental\")", help_text)
        @test occursin("DepolLib.help(\"instrumental\")", help_text)
        @test occursin("run_instrumental_effect_job(cfg)", one_job_text)
        @test !occursin("run_mach_suite(cfg)", one_job_text)
        @test occursin("run_canal_metrics_job(cfg)", canal_alias_text)
        @test !occursin("run_instrumental_effect_job(cfg)", canal_alias_text)
        @test occursin("run_mach_suite(cfg)", jobs_text)
        @test occursin("run_canal_metrics_job(cfg)", jobs_text)
        @test occursin("src/code/jobs/run_canal_metrics_job.jl --config config/default.toml", cli_text)
        @test occursin("for job in \\", cli_text)
        @test occursin("1. ne_dm_em", order_text)
        @test occursin("6. instrumental", order_text)
        @test_throws ErrorException depol_help(IOBuffer(); topic="unknown_topic")
        @test jobs_text == sprint(io -> DepolLib.help(io; topic="jobs"))
        @test help_text == alias_text
    end
end
