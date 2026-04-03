using Test
using Logging

module InstrumentalEffectJobUnderTest
include(joinpath(@__DIR__, "..", "src", "code", "jobs", "run_instrumental_effect_job.jl"))
end

const IEJ = InstrumentalEffectJobUnderTest
const _JOB_PIPELINE_CAPTURE = Ref{Any}(nothing)

@eval IEJ.InstrumentalEffect begin
    function run_pipeline(cfg::InstrumentalConfig; flags::RunFlags=RunFlags(), on_step::Function=(_ -> nothing))
        Main._JOB_PIPELINE_CAPTURE[] = (cfg=cfg, flags=flags)
        if flags.run_phi_q_u_p
            @info "Phi channel selected" iphi=cfg.iphi phi_rad_m2=cfg.PhiArray[cfg.iphi]
        end
        return (integrity=Dict{String,Any}(
            "status" => "pass",
            "report_path" => joinpath(cfg.base_out, "integrity_report.txt"),
            "critical_failures" => 0,
            "warnings" => 0,
        ),)
    end
end

function _write_job_fits(path::String, data)
    mkpath(dirname(path))
    IEJ.write_FITS(path, data)
end

function _make_job_cfg(td::String; phi_min=-1.0, phi_max=1.0, dphi=0.5, phi_target=0.0,
        run_phi_q_u_p::Bool=false, run_channel_b_alignment::Bool=false)
    sim = "simA"
    los = "y"
    withfaraday = joinpath(td, sim, los, "Synchrotron", "WithFaraday")

    phi_array = collect(phi_min:dphi:phi_max)

    _write_job_fits(joinpath(withfaraday, "Qnu.fits"), reshape(Float32[1], 1, 1, 1))
    _write_job_fits(joinpath(withfaraday, "Unu.fits"), reshape(Float32[1], 1, 1, 1))
    _write_job_fits(joinpath(withfaraday, "Pmax.fits"), reshape(Float32[1], 1, 1))
    _write_job_fits(joinpath(withfaraday, "realFDF.fits"), reshape(Float32[1], 1, 1, 1))
    _write_job_fits(joinpath(withfaraday, "imagFDF.fits"), reshape(Float32[1], 1, 1, 1))
    _write_job_fits(joinpath(td, sim, "Bx.fits"), reshape(Float32[1], 1, 1, 1))
    _write_job_fits(joinpath(td, sim, "By.fits"), reshape(Float32[1], 1, 1, 1))
    _write_job_fits(joinpath(td, sim, "Bz.fits"), reshape(Float32[1], 1, 1, 1))
    _write_job_fits(joinpath(td, sim, "density.fits"), reshape(Float32[1], 1, 1, 1))

    return Dict(
        "paths" => Dict(
            "simulations_root" => td,
            "desktop_output_root" => joinpath(td, "outputs"),
        ),
        "simulation" => Dict(
            "name" => sim,
            "los" => los,
        ),
        "tasks" => Dict(
            "instrumental_effect" => Dict(
                "llarge_list" => [4.0],
                "phi_min" => phi_min,
                "phi_max" => phi_max,
                "dphi" => dphi,
                "phi_target" => phi_target,
                "run_pmax_maps" => false,
                "run_psd" => false,
                "run_q_u_p_q2" => false,
                "run_phi_q_u_p" => run_phi_q_u_p,
                "run_lic" => false,
                "run_channel_b_alignment" => run_channel_b_alignment,
                "channel_alignment_pdf_plain" => true,
            ),
        ),
    )
end

function _phi_selection_log(logger::Test.TestLogger)
    matches = filter(rec -> rec.message == "Phi channel selected", logger.logs)
    isempty(matches) && error("Missing `Phi channel selected` log record")
    return only(matches)
end

@testset "Instrumental job configuration and orchestration" begin
    @testset "phi target selection uses nearest index" begin
        mktempdir() do td
            cfg = _make_job_cfg(td; phi_min=-1.0, phi_max=1.0, dphi=0.5, phi_target=0.2, run_phi_q_u_p=true)
            logger = Test.TestLogger(min_level=Logging.Info)
            result = nothing
            with_logger(logger) do
                _JOB_PIPELINE_CAPTURE[] = nothing
                result = IEJ.run_instrumental_effect_job(cfg)
            end

            record = _phi_selection_log(logger)
            captured = _JOB_PIPELINE_CAPTURE[]
            @test result["task"] == "instrumental"
            @test result["flags"]["run_phi_q_u_p"] == true
            @test result["integrity_status"] == "pass"
            @test captured !== nothing
            @test captured.cfg.PhiArray == [-1.0, -0.5, 0.0, 0.5, 1.0]
            @test captured.cfg.iphi == 3
            @test get(record.kwargs, :iphi, nothing) == 3
            @test isapprox(get(record.kwargs, :phi_rad_m2, NaN), 0.0; atol=1e-12)
            @test result["output_dir"] == captured.cfg.base_out
        end
    end

    @testset "phi target outside range clamps to nearest edge" begin
        mktempdir() do td
            cfg = _make_job_cfg(td; phi_min=-1.0, phi_max=1.0, dphi=0.5, phi_target=9.0, run_phi_q_u_p=true)
            logger = Test.TestLogger(min_level=Logging.Info)
            with_logger(logger) do
                _JOB_PIPELINE_CAPTURE[] = nothing
                IEJ.run_instrumental_effect_job(cfg)
            end

            record = _phi_selection_log(logger)
            captured = _JOB_PIPELINE_CAPTURE[]
            @test captured !== nothing
            @test captured.cfg.iphi == 5
            @test get(record.kwargs, :iphi, nothing) == 5
            @test isapprox(get(record.kwargs, :phi_rad_m2, NaN), 1.0; atol=1e-12)
        end
    end

    @testset "empty phi range errors clearly" begin
        cfg = Dict(
            "paths" => Dict(
                "simulations_root" => "/tmp/unused",
                "desktop_output_root" => "/tmp/unused_out",
            ),
            "simulation" => Dict(
                "name" => "simA",
                "los" => "y",
            ),
            "tasks" => Dict(
                "instrumental_effect" => Dict(
                    "phi_min" => -1.0,
                    "phi_max" => 1.0,
                    "dphi" => -0.5,
                ),
            ),
        )

        err = try
            IEJ.run_instrumental_effect_job(cfg)
            nothing
        catch ex
            ex
        end
        @test err isa ErrorException
        @test occursin("Empty PhiArray built", sprint(showerror, err))
    end

    @testset "zero phi step is rejected" begin
        cfg = Dict(
            "paths" => Dict(
                "simulations_root" => "/tmp/unused",
                "desktop_output_root" => "/tmp/unused_out",
            ),
            "simulation" => Dict(
                "name" => "simA",
                "los" => "y",
            ),
            "tasks" => Dict(
                "instrumental_effect" => Dict(
                    "phi_min" => -1.0,
                    "phi_max" => 1.0,
                    "dphi" => 0.0,
                ),
            ),
        )

        @test_throws ArgumentError IEJ.run_instrumental_effect_job(cfg)
    end
end
