using Test

function _write_basic_fits(path::String, data)
    mkpath(dirname(path))
    write_FITS(path, data)
end

struct _StopPipeline <: Exception
    step::Symbol
end

function _with_rm_progress_disabled(f::Function)
    old = get(ENV, "DEPOL_RM_LOG_PROGRESS", nothing)
    ENV["DEPOL_RM_LOG_PROGRESS"] = "0"
    try
        return f()
    finally
        if old === nothing
            delete!(ENV, "DEPOL_RM_LOG_PROGRESS")
        else
            ENV["DEPOL_RM_LOG_PROGRESS"] = old
        end
    end
end

function _make_minimal_cfg(td::String; bad_layout::Bool=false, pmax_nan::Bool=false)
    q_path = joinpath(td, "Qnu.fits")
    u_path = joinpath(td, "Unu.fits")
    pmax_path = joinpath(td, "Pmax.fits")
    qphi_path = joinpath(td, "realFDF.fits")
    uphi_path = joinpath(td, "imagFDF.fits")
    bx_path = joinpath(td, "Bx.fits")
    by_path = joinpath(td, "By.fits")
    bz_path = joinpath(td, "Bz.fits")
    dens_path = joinpath(td, "density.fits")

    q_shape = bad_layout ? (4, 4, 3) : (4, 4, 4)
    _write_basic_fits(q_path, randn(Float32, q_shape...))
    _write_basic_fits(u_path, randn(Float32, q_shape...))
    pmax = rand(Float32, 4, 4)
    if pmax_nan
        pmax[1, 1] = NaN32
    end
    _write_basic_fits(pmax_path, pmax)
    _write_basic_fits(qphi_path, randn(Float32, 4, 4, 3))
    _write_basic_fits(uphi_path, randn(Float32, 4, 4, 3))
    _write_basic_fits(bx_path, randn(Float32, 4, 4, 4))
    _write_basic_fits(by_path, randn(Float32, 4, 4, 4))
    _write_basic_fits(bz_path, randn(Float32, 4, 4, 4))
    _write_basic_fits(dens_path, rand(Float32, 4, 4, 4))

    return InstrumentalConfig(
        Q_in=q_path,
        U_in=u_path,
        Pmax_nofilter_path=pmax_path,
        Q_in_phi=qphi_path,
        U_in_phi=uphi_path,
        base_out=joinpath(td, "out"),
        Bx_in=bx_path,
        By_in=by_path,
        Bz_in=bz_path,
        dens_in=dens_path,
        los="y",
        n=4,
        m=4,
        Lbox_pc=4.0,
        Lcut_small=1.0,
        Llarge_list=[4.0],
        νmin_MHz=100.0,
        νmax_MHz=103.0,
        Δν_MHz=1.0,
        PhiArray=[-1.0, 0.0, 1.0],
        ichan=2,
        iphi=2,
    )
end

@testset "Instrumental filter pass" begin
    @testset "smoke run with synthetic FITS" begin
        mktempdir() do td
            cfg = _make_minimal_cfg(td)

            result = _with_rm_progress_disabled() do
                run_filter_pass(cfg)
            end
            @test result isa InstrumentalEffect.FilterPassResult
            @test result.L_ok == [4.0]
            @test size(result.Pmax0) == (4, 4)
            @test result.integrity["status"] == "pass"
            @test isfile(result.integrity["report_path"])
            @test haskey(result.Qslice_filt, 4.0)
            @test size(result.Qslice_filt[4.0]) == (4, 4)
            @test isfile(joinpath(cfg.base_out, "HardBandPass_remove_L0_to_1pc_and_4to50pc", "Qnu_filtered.fits"))
            @test isfile(joinpath(cfg.base_out, "HardBandPass_remove_L0_to_1pc_and_4to50pc", "RMSynthesis", "Pphi_max.fits"))
        end
    end

    @testset "NaN in critical map aborts run" begin
        mktempdir() do td
            cfg = _make_minimal_cfg(td; pmax_nan=true)

            err = try
                _with_rm_progress_disabled() do
                    run_filter_pass(cfg)
                end
                nothing
            catch ex
                ex
            end

            report_path = joinpath(cfg.base_out, "integrity_report.txt")
            @test err isa ErrorException
            @test occursin("Critical integrity checks failed", sprint(showerror, err))
            @test isfile(report_path)
            @test occursin("nan_inf.Pmax", read(report_path, String))
        end
    end

    @testset "invalid frequency layout is rejected" begin
        mktempdir() do td
            cfg = _make_minimal_cfg(td; bad_layout=true)

            err = try
                _with_rm_progress_disabled() do
                    run_filter_pass(cfg)
                end
                nothing
            catch ex
                ex
            end

            @test err isa ErrorException
            @test occursin("dims.stokes_layout", sprint(showerror, err))
        end
    end

    @testset "run_pipeline flags and section gating" begin
        @testset "all optional sections disabled" begin
            mktempdir() do td
                cfg = _make_minimal_cfg(td)
                steps = Symbol[]
                flags = RunFlags(false, false, false, false, false)
                result = _with_rm_progress_disabled() do
                    run_pipeline(cfg; flags=flags, on_step=step -> push!(steps, step))
                end
                @test result isa InstrumentalEffect.FilterPassResult
                @test steps == [:filter_pass]
            end
        end

        @testset "LIC-only section enabled" begin
            mktempdir() do td
                cfg = _make_minimal_cfg(td)
                steps = Symbol[]
                flags = RunFlags(false, false, false, false, true)
                result = _with_rm_progress_disabled() do
                    run_pipeline(cfg; flags=flags, on_step=step -> push!(steps, step))
                end
                @test result isa InstrumentalEffect.FilterPassResult
                @test steps == [:filter_pass, :run_lic]
            end
        end

        @testset "pmax section hook reached only when enabled" begin
            mktempdir() do td
                cfg = _make_minimal_cfg(td)
                steps = Symbol[]

                err = try
                    _with_rm_progress_disabled() do
                        run_pipeline(cfg;
                            flags=RunFlags(true, false, false, false, false),
                            on_step=step -> begin
                                push!(steps, step)
                                step == :run_pmax_maps && throw(_StopPipeline(step))
                            end,
                        )
                    end
                    nothing
                catch ex
                    ex
                end

                @test err isa _StopPipeline
                @test err.step == :run_pmax_maps
                @test steps == [:filter_pass, :run_pmax_maps]
            end
        end

        @testset "channel-angle section enabled writes outputs" begin
            mktempdir() do td
                cfg = _make_minimal_cfg(td)
                steps = Symbol[]
                flags = RunFlags(
                    run_pmax_maps=false,
                    run_psd=false,
                    run_q_u_p_q2=false,
                    run_phi_q_u_p=false,
                    run_lic=false,
                    run_channel_b_alignment=true,
                )

                result = _with_rm_progress_disabled() do
                    run_pipeline(cfg; flags=flags, on_step=step -> push!(steps, step))
                end

                @test result isa InstrumentalEffect.FilterPassResult
                @test steps == [:filter_pass, :run_channel_b_alignment]
                @test isfile(joinpath(cfg.base_out, "channel_alignment_summary.csv"))
                @test isfile(joinpath(cfg.base_out, "figures", "channel_alignment_vs_filter.pdf"))

                csv = read(joinpath(cfg.base_out, "channel_alignment_summary.csv"), String)
                @test occursin("reference", csv)
                @test occursin("B_perp", csv)
            end
        end
    end
end
