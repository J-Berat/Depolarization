using Test
using Statistics

function _instrumental_cfg_stub(; los::String="y", n::Int=4, m::Int=4, lbox_pc::Float64=4.0)
    return InstrumentalEffect.InstrumentalConfig(
        Q_in="",
        U_in="",
        Pmax_nofilter_path="",
        Q_in_phi="",
        U_in_phi="",
        base_out="",
        Bx_in="",
        By_in="",
        Bz_in="",
        dens_in="",
        los=los,
        n=n,
        m=m,
        Lbox_pc=lbox_pc,
        Llarge_list=[20.0],
        PhiArray=[-1.0, 0.0, 1.0],
    )
end

function _main_canal_mask_without_fallback(Pmap::AbstractMatrix, Δx::Real, Δy::Real)
    score_raw = InstrumentalEffect._canal_score_map(Pmap)
    score = InstrumentalEffect._box_mean(score_raw, InstrumentalEffect._CANAL_RIDGE_SMOOTH_RADIUS_PIX)
    thr = InstrumentalEffect._finite_quantile(score, InstrumentalEffect._CANAL_SCORE_QUANTILE)
    gx, gy = InstrumentalEffect._gradients_central(score, Δx, Δy)

    Jxx = InstrumentalEffect._box_mean(gx .^ 2, InstrumentalEffect._CANAL_TENSOR_SMOOTH_RADIUS_PIX)
    Jxy = InstrumentalEffect._box_mean(gx .* gy, InstrumentalEffect._CANAL_TENSOR_SMOOTH_RADIUS_PIX)
    Jyy = InstrumentalEffect._box_mean(gy .^ 2, InstrumentalEffect._CANAL_TENSOR_SMOOTH_RADIUS_PIX)

    mask = falses(size(score))
    @inbounds for j in axes(score, 2), i in axes(score, 1)
        Dij = score[i, j]
        isfinite(Dij) || continue
        Dij >= thr || continue

        λtan, λnorm, tvecx, tvecy, nvecx, nvecy = InstrumentalEffect._symmetric_2x2_eigensystem(
            Jxx[i, j], Jxy[i, j], Jyy[i, j]
        )
        λnorm > 0 || continue
        anis = λnorm / max(λtan, 1e-6)
        anis >= InstrumentalEffect._CANAL_RIDGE_ANIS_MIN || continue

        abs(tvecx) + abs(tvecy) > 1e-12 || continue
        abs(nvecx) + abs(nvecy) > 1e-12 || continue

        sx = round(Int, clamp(sign(nvecx), -1, 1))
        sy = round(Int, clamp(sign(nvecy), -1, 1))
        (sx != 0 || sy != 0) || continue

        i1 = clamp(i - sx, first(axes(score, 1)), last(axes(score, 1)))
        i2 = clamp(i + sx, first(axes(score, 1)), last(axes(score, 1)))
        j1 = clamp(j - sy, first(axes(score, 2)), last(axes(score, 2)))
        j2 = clamp(j + sy, first(axes(score, 2)), last(axes(score, 2)))
        (Dij >= score[i1, j1] && Dij >= score[i2, j2]) || continue

        mask[i, j] = true
    end

    return mask
end

function _alignment_delta_median_deg(θ::AbstractMatrix, mask::BitMatrix, expected::Real)
    θref = fill(Float64(expected), size(θ))
    deltas = abs.(InstrumentalEffect._alignment_deltas_deg(θ, θref, mask))
    return isempty(deltas) ? NaN : median(deltas)
end

@testset "Instrumental spectra and helper utilities" begin
    @testset "PSD normalization and mean removal" begin
        img = fill(3.0, 4, 3)
        kx = [-2.0, -1.0, 0.0, 1.0]

        k_keep, p_keep = InstrumentalEffect.psd1d_x_mean_over_y(img, kx; remove_mean=false)
        _, p_rm = InstrumentalEffect.psd1d_x_mean_over_y(img, kx; remove_mean=true)

        @test k_keep == kx
        @test isapprox(p_keep[3], 36.0; atol=1e-12)
        @test all(isapprox.(p_keep[[1, 2, 4]], 0.0; atol=1e-12))
        @test all(isapprox.(p_rm, 0.0; atol=1e-12))

        img2 = reshape(Float64.(1:12), 3, 4)
        p2 = InstrumentalEffect.psd2d(img2)
        @test isapprox(sum(p2), sum(abs2, img2); rtol=1e-12, atol=1e-12)
    end

    @testset "label helpers" begin
        @test InstrumentalEffect.Llabel_pc(20.0) == 4
        @test InstrumentalEffect.Llabel_pc(30.0) == 6
        @test InstrumentalEffect._filter_display_label("nofilter", 0.0) == "No filter"
        @test InstrumentalEffect._filter_display_label("filtered", 3.6) == "L_eff = 4 pc"
        @test string(InstrumentalEffect._filter_display_label_latex("nofilter", 0.0)) == "\\mathrm{No\\ filter}"
        @test string(InstrumentalEffect._filter_display_label_latex("filtered", 3.6)) == "L_{\\mathrm{eff}}=4\\,\\mathrm{pc}"
    end

    @testset "plot axis helpers" begin
        labels_x = InstrumentalEffect._sky_plane_labels(_instrumental_cfg_stub(los="x"))
        labels_y = InstrumentalEffect._sky_plane_labels(_instrumental_cfg_stub(los="y"))
        labels_z = InstrumentalEffect._sky_plane_labels(_instrumental_cfg_stub(los="z"))
        @test string(labels_x[1]) == "y\\,[\\mathrm{pc}]"
        @test string(labels_x[2]) == "z\\,[\\mathrm{pc}]"
        @test string(labels_y[1]) == "z\\,[\\mathrm{pc}]"
        @test string(labels_y[2]) == "x\\,[\\mathrm{pc}]"
        @test string(labels_z[1]) == "x\\,[\\mathrm{pc}]"
        @test string(labels_z[2]) == "y\\,[\\mathrm{pc}]"

        @test InstrumentalEffect._cell_edges_for_plot([0.0, 2.0, 4.0]) == [-1.0, 1.0, 3.0, 5.0]
        @test_throws ErrorException InstrumentalEffect._cell_edges_for_plot([1.0])

        ux, uy = InstrumentalEffect._bperp_arrow_unit_vector(ones(3, 3), zeros(3, 3))
        @test isapprox(ux, 0.0; atol=1e-12)
        @test isapprox(uy, 1.0; atol=1e-12)

        ux_nan, uy_nan = InstrumentalEffect._bperp_arrow_unit_vector(fill(NaN, 2, 2), fill(NaN, 2, 2))
        @test isapprox(ux_nan, 1.0; atol=1e-12)
        @test isapprox(uy_nan, 0.0; atol=1e-12)
    end

    @testset "orientation and projection helpers" begin
        θ_h = InstrumentalEffect._orientation_map_bperp_from_b1_b2(zeros(3, 3), ones(3, 3))
        @test maximum(abs.(θ_h)) < 1e-12

        θ_v = InstrumentalEffect._orientation_map_bperp_from_b1_b2(ones(3, 3), zeros(3, 3))
        @test maximum(abs.(θ_v .- (π / 2))) < 1e-12

        θ_d = InstrumentalEffect._orientation_map_bperp_from_b1_b2(ones(3, 3), ones(3, 3))
        @test maximum(abs.(θ_d .- (π / 4))) < 1e-12

        Bx = fill(1.0f0, 2, 2, 2)
        By = fill(2.0f0, 2, 2, 2)
        Bz = fill(3.0f0, 2, 2, 2)
        b1x, b2x = InstrumentalEffect._select_bsky_component_cubes(Bx, By, Bz, "x")
        b1y, b2y = InstrumentalEffect._select_bsky_component_cubes(Bx, By, Bz, "y")
        b1z, b2z = InstrumentalEffect._select_bsky_component_cubes(Bx, By, Bz, "z")
        @test b1x === By && b2x === Bz
        @test b1y === Bx && b2y === Bz
        @test b1z === Bx && b2z === By
        @test_throws ErrorException InstrumentalEffect._select_bsky_component_cubes(Bx, By, Bz, "bad")

        B1_cube = zeros(Float64, 2, 2, 1)
        B2_cube = zeros(Float64, 2, 2, 1)
        B1_cube[1, :, 1] .= [3.0, 0.0]
        B2_cube[1, :, 1] .= [0.0, 4.0]
        B1_cube[2, :, 1] .= [5.0, 0.0]
        B2_cube[2, :, 1] .= [0.0, 12.0]

        bperp = InstrumentalEffect._project_hypot_mean(B1_cube, B2_cube, "y", (2, 1))
        @test size(bperp) == (2, 1)
        @test isapprox(bperp[1, 1], 3.5; atol=1e-12)
        @test isapprox(bperp[2, 1], 8.5; atol=1e-12)

        b1mean = InstrumentalEffect._project_cube_mean(B1_cube, "y", (2, 1))
        b2mean = InstrumentalEffect._project_cube_mean(B2_cube, "y", (2, 1))
        @test bperp[1, 1] > hypot(b1mean[1, 1], b2mean[1, 1])
        @test bperp[2, 1] > hypot(b1mean[2, 1], b2mean[2, 1])
    end

    @testset "canal score and orientation" begin
        Pscore = ones(9, 9)
        Pscore[:, 5] .= 0.05
        Pscore[2, 2] = NaN
        score = InstrumentalEffect._canal_score_map(Pscore)
        @test isnan(score[2, 2])
        @test score[5, 5] > score[1, 1]

        Ph = ones(21, 21)
        Ph[:, 11] .= 0.0
        θh, maskh, _, _ = InstrumentalEffect._canal_orientation_map(Ph, 1.0, 1.0)
        @test count(maskh) > 0
        @test _alignment_delta_median_deg(θh, maskh, 0.0) <= 10.0

        Pv = ones(21, 21)
        Pv[11, :] .= 0.0
        θv, maskv, _, _ = InstrumentalEffect._canal_orientation_map(Pv, 1.0, 1.0)
        @test count(maskv) > 0
        @test _alignment_delta_median_deg(θv, maskv, π / 2) <= 10.0
    end

    @testset "canal fallback path on isotropic blob" begin
        Pblob = Matrix{Float64}(undef, 31, 31)
        for j in axes(Pblob, 2), i in axes(Pblob, 1)
            r2 = (i - 16)^2 + (j - 16)^2
            Pblob[i, j] = 1.0 - exp(-r2 / 18.0)
        end

        main_mask = _main_canal_mask_without_fallback(Pblob, 1.0, 1.0)
        @test count(main_mask) == 0

        θblob, maskblob, _, _ = InstrumentalEffect._canal_orientation_map(Pblob, 1.0, 1.0)
        @test count(maskblob) > 0
        @test all(isfinite, θblob[maskblob])
    end

    @testset "channel alignment summary on controlled geometry" begin
        mktempdir() do td
            n = 21
            m = 21
            cfg = InstrumentalEffect.InstrumentalConfig(
                Q_in="",
                U_in="",
                Pmax_nofilter_path="",
                Q_in_phi="",
                U_in_phi="",
                base_out=td,
                Bx_in=joinpath(td, "Bx.fits"),
                By_in=joinpath(td, "By.fits"),
                Bz_in=joinpath(td, "Bz.fits"),
                dens_in=joinpath(td, "density.fits"),
                los="y",
                n=n,
                m=m,
                Lbox_pc=21.0,
                Llarge_list=[20.0],
            )

            Bx = zeros(Float32, n, 2, m)
            By = zeros(Float32, n, 2, m)
            Bz = ones(Float32, n, 2, m)
            write_FITS(cfg.Bx_in, Bx)
            write_FITS(cfg.By_in, By)
            write_FITS(cfg.Bz_in, Bz)
            write_FITS(cfg.dens_in, ones(Float32, n, 2, m))

            Pmap = ones(Float32, n, m)
            Pmap[:, 11] .= 0.0f0

            scales = InstrumentalEffect.grid_scales(cfg)
            axes = InstrumentalEffect.spectral_axes(cfg, scales.Δx, scales.Δy)
            data = InstrumentalEffect.FilterPassResult(
                zeros(Float32, n, m, 1),
                zeros(Float32, n, m, 1),
                Dict(20.0 => copy(Pmap)),
                Dict(20.0 => copy(Pmap)),
                copy(Pmap),
                Dict(20.0 => copy(Pmap)),
                [20.0],
                scales,
                axes,
                [1.0],
                [-1.0, 0.0, 1.0],
                Dict{String,Any}("status" => "pass"),
            )

            result = InstrumentalEffect._run_channel_b_alignment(cfg, data)
            @test isfile(result.summary_path)
            @test isfile(result.figure_path)
            @test length(result.rows) == 2

            row_nofilter = only(filter(r -> r.map_tag == "nofilter", result.rows))
            row_filtered = only(filter(r -> r.map_tag != "nofilter", result.rows))
            @test row_nofilter.reference == "B_perp"
            @test row_nofilter.npix > 0
            @test row_nofilter.frac_parallel_15deg >= 0.95
            @test row_nofilter.median_abs_delta_deg <= 5.0
            @test row_filtered.llarge_eff_pc == 20.0
            @test row_filtered.frac_parallel_15deg >= 0.95

            csv = read(result.summary_path, String)
            @test occursin("nofilter,0.000000,B_perp", csv)
            @test occursin("HardBandPass_remove_L0_to_1pc_and_20to50pc,20.000000,B_perp", csv)
        end
    end

    @testset "canal morphology on controlled straight filaments" begin
        Pmap = Matrix{Float64}(undef, 61, 61)
        σ1 = 1.6
        @inbounds for j in axes(Pmap, 2), i in axes(Pmap, 1)
            v1 = 0.90 * exp(-((j - 20)^2) / (2 * σ1^2))
            Pmap[i, j] = 1.0 - v1
        end

        result = InstrumentalEffect._analyze_canal_morphology_map(Pmap, 1.0, 1.0)
        @test length(result.rows) >= 1

        rows = sort(result.rows, by=row -> row.length_pc, rev=true)
        @test rows[1].length_pc >= 40.0
        @test rows[1].width_fwhm_median_pc > 0
        @test rows[1].curvature_median_pc_inv <= 0.2

        summary = InstrumentalEffect._build_canal_morphology_summary(result.rows, 1.0, 1.0)
        @test length(summary) == 1
        @test summary[1].n_branches >= 1
        @test isfinite(summary[1].width_xmin_pc)
    end
end
