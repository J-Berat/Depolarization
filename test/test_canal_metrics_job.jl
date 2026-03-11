using Test
using Statistics

module CanalMetricsJobUnderTest
include(joinpath(@__DIR__, "..", "src", "code", "jobs", "run_canal_metrics_job.jl"))
end

const CMJ = CanalMetricsJobUnderTest

@testset "Canal metrics helpers" begin
    function _map_shape_for_los(cube::AbstractArray{<:Real,3}, los::AbstractString)
        if los == "x"
            return size(cube, 2), size(cube, 3)
        elseif los == "y"
            return size(cube, 1), size(cube, 3)
        else
            return size(cube, 1), size(cube, 2)
        end
    end

    function _manual_cb_map(Bpar::AbstractArray{<:Real,3}, los::AbstractString)
        _, profile_fun = DepolLib.los_config(los)
        nx, ny = _map_shape_for_los(Bpar, los)
        out = zeros(Float64, nx, ny)
        for j in 1:ny, i in 1:nx
            prof = Float64.(profile_fun(i, j, Bpar))
            out[i, j] = abs(sum(prof)) / (sum(abs.(prof)) + 1e-30)
        end
        return out
    end

    function _manual_cphi_map(ne::AbstractArray{<:Real,3}, Bpar::AbstractArray{<:Real,3}, los::AbstractString)
        _, profile_fun = DepolLib.los_config(los)
        nx, ny = _map_shape_for_los(Bpar, los)
        out = zeros(Float64, nx, ny)
        for j in 1:ny, i in 1:nx
            ne_prof = Float64.(profile_fun(i, j, ne))
            b_prof = Float64.(profile_fun(i, j, Bpar))
            x = ne_prof .* b_prof
            out[i, j] = abs(sum(x)) / (sum(abs.(x)) + 1e-30)
        end
        return out
    end

    function _manual_nrev_map(Bpar::AbstractArray{<:Real,3}, los::AbstractString; thresh::Real=0.0)
        _, profile_fun = DepolLib.los_config(los)
        nx, ny = _map_shape_for_los(Bpar, los)
        out = zeros(Int16, nx, ny)
        for j in 1:ny, i in 1:nx
            prof = Float64.(profile_fun(i, j, Bpar))
            prev = 0
            n = 0
            for v in prof
                s = abs(v) <= thresh ? 0 : (v > 0 ? 1 : -1)
                if prev != 0 && s != 0 && s != prev
                    n += 1
                end
                if s != 0
                    prev = s
                end
            end
            out[i, j] = Int16(n)
        end
        return out
    end

    Bx = reshape(Float32.(1:24), 2, 3, 4)
    By = Bx .- 8.0f0
    Bz = Bx .* (-1.0f0)
    ne = reshape(Float32.(range(0.5, 3.0; length=24)), 2, 3, 4)

    for los in ("x", "y", "z")
        Bpar = CMJ._bparallel_from_components(Bx, By, Bz, los)
        cb_expected = _manual_cb_map(Bpar, los)
        cphi_expected = _manual_cphi_map(ne, Bpar, los)
        nrev_expected = _manual_nrev_map(Bpar, los)

        @test CMJ._cb_map(Bpar, los) ≈ cb_expected atol=1e-7 rtol=1e-7
        @test CMJ._cphi_map(ne, Bpar, los) ≈ cphi_expected atol=1e-7 rtol=1e-7
        @test CMJ._nrev_map(Bpar, los) == nrev_expected
    end

    Pmax = Float32[1 2 3; 4 5 6; 7 8 9]
    mask, thr = CMJ._canal_mask(Pmax; decile=0.25)
    @test isapprox(thr, quantile(vec(Float64.(Pmax)), 0.25); atol=1e-12)
    @test count(mask) >= 2
end
