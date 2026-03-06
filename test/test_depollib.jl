using Test

@testset "DepolLib config and helpers" begin
    @testset "runtime config parsing" begin
        mktempdir() do td
            cfg_path = joinpath(td, "cfg.toml")
            open(cfg_path, "w") do io
                println(io, "[paths]")
                println(io, "simulations_root = \"$td\"")
                println(io, "outputs_root = \"$td/out\"")
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
                println(io, "outputs_root = \"$td/out\"")
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
                println(io, "outputs_root = \"$td/out\"")
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
end
