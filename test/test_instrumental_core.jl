using Test

@testset "Instrumental core contracts" begin
    @testset "channel extraction" begin
        cube_last = reshape(Float64.(1:48), 4, 3, 4)
        slice_last = InstrumentalEffect.get_chan_xy(cube_last, 2, 4, 3)
        @test size(slice_last) == (4, 3)

        cube_first = permutedims(cube_last, (3, 1, 2))
        slice_first = InstrumentalEffect.get_chan_xy(cube_first, 2, 4, 3)
        @test slice_first == slice_last

        @test_throws ErrorException InstrumentalEffect.get_chan_xy(cube_last, 99, 4, 3)
    end

    @testset "RMSynthesis shapes" begin
        nu = [1.0e8, 1.1e8, 1.2e8, 1.3e8]
        phi = [-1.0, 0.0, 1.0]

        q1 = randn(4)
        u1 = randn(4)
        absf1, _, _ = InstrumentalEffect.RMSynthesis(q1, u1, nu, phi; log_progress=false)
        @test size(absf1) == (3,)

        q2 = randn(5, 4)
        u2 = randn(5, 4)
        absf2, _, _ = InstrumentalEffect.RMSynthesis(q2, u2, nu, phi; log_progress=false)
        @test size(absf2) == (5, 3)

        q3 = randn(2, 3, 4)
        u3 = randn(2, 3, 4)
        absf3, _, _ = InstrumentalEffect.RMSynthesis(q3, u3, nu, phi; log_progress=false)
        @test size(absf3) == (2, 3, 3)
    end

    @testset "RMSynthesis pmax-only path" begin
        nu = [1.0e8, 1.1e8, 1.2e8, 1.3e8]
        phi = [-1.0, 0.0, 1.0, 2.0]

        q3 = randn(Float32, 2, 3, 4)
        u3 = randn(Float32, 2, 3, 4)
        absf3, _, _ = InstrumentalEffect.RMSynthesis(q3, u3, nu, phi; log_progress=false)
        pmax_ref = Float32.(InstrumentalEffect.Pphi_max_map(absf3))
        pmax_fast = InstrumentalEffect.RMSynthesis_pmax_map(q3, u3, nu, phi; log_progress=false)
        @test isapprox(pmax_fast, pmax_ref; rtol=1e-5, atol=1e-6)

        q3f = permutedims(q3, (3, 1, 2))
        u3f = permutedims(u3, (3, 1, 2))
        absf3f, _, _ = InstrumentalEffect.RMSynthesis(q3f, u3f, nu, phi; log_progress=false)
        pmax_ref_f = Float32.(InstrumentalEffect.Pphi_max_map(absf3f))
        pmax_fast_f = InstrumentalEffect.RMSynthesis_pmax_map(q3f, u3f, nu, phi; log_progress=false)
        @test isapprox(pmax_fast_f, pmax_ref_f; rtol=1e-5, atol=1e-6)
    end

    @testset "filter application guards" begin
        img = randn(4, 4)
        H = ones(Float32, 4, 4)
        out = InstrumentalEffect.apply_to_array_xy(img, H; n=4, m=4)
        @test size(out) == (4, 4)

        @test_throws ErrorException InstrumentalEffect.apply_to_array_xy(img, ones(Float32, 3, 3); n=4, m=4)
    end
end
