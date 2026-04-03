"""
    Shared magnetic-field projection helpers used by plotting and alignment.
"""

"""
    _project_cube_mean(...)

LOS-average projection with map-shape reconciliation.
"""
function _project_cube_mean(cube, los::String, target_size::Tuple{Int,Int})
    ndims(cube) == 3 || error("Expected 3D cube for LOS projection, got ndims=$(ndims(cube)) size=$(size(cube))")
    los in ("x", "y", "z") || error("LOS must be x/y/z, got $los")

    dim = los == "x" ? 1 : (los == "y" ? 2 : 3)
    proj = dropdims(mean(Float64.(cube); dims=dim), dims=dim)

    if size(proj) == target_size
        return Matrix{Float64}(proj)
    end

    proj_t = permutedims(proj, (2, 1))
    if size(proj_t) == target_size
        return Matrix{Float64}(proj_t)
    end

    error("Projected map has size $(size(proj)) (or transposed $(size(proj_t))), expected $target_size")
end

"""
    _select_bsky_component_cubes(...)

Returns the two magnetic-field cubes that lie in the sky plane for the chosen LOS.
"""
function _select_bsky_component_cubes(Bx, By, Bz, los::String)
    los in ("x", "y", "z") || error("LOS must be x/y/z, got $los")
    if los == "x"
        return By, Bz
    elseif los == "y"
        return Bx, Bz
    end
    return Bx, By
end

"""
    _project_hypot_mean(...)

LOS-average projection of
`B_perp(B1,B2) = sqrt(B1^2 + B2^2)`
without allocating a full intermediate cube.
"""
function _project_hypot_mean(B1_cube, B2_cube, los::String, target_size::Tuple{Int,Int})
    size(B1_cube) == size(B2_cube) || error("B1/B2 size mismatch: B1=$(size(B1_cube)) B2=$(size(B2_cube))")
    ndims(B1_cube) == 3 || error("B1_cube must be 3D, got ndims=$(ndims(B1_cube))")

    nx, ny, nz = size(B1_cube)
    if los == "x"
        out = Matrix{Float64}(undef, ny, nz)
        @inbounds for j in 1:nz, i in 1:ny
            acc = 0.0
            for k in 1:nx
                acc += hypot(B1_cube[k, i, j], B2_cube[k, i, j])
            end
            out[i, j] = acc / nx
        end
    elseif los == "y"
        out = Matrix{Float64}(undef, nx, nz)
        @inbounds for j in 1:nz, i in 1:nx
            acc = 0.0
            for k in 1:ny
                acc += hypot(B1_cube[i, k, j], B2_cube[i, k, j])
            end
            out[i, j] = acc / ny
        end
    elseif los == "z"
        out = Matrix{Float64}(undef, nx, ny)
        @inbounds for j in 1:ny, i in 1:nx
            acc = 0.0
            for k in 1:nz
                acc += hypot(B1_cube[i, j, k], B2_cube[i, j, k])
            end
            out[i, j] = acc / nz
        end
    else
        error("LOS must be x/y/z, got $los")
    end

    if size(out) == target_size
        return out
    end

    out_t = permutedims(out, (2, 1))
    if size(out_t) == target_size
        return out_t
    end
    error("Projected hypot map has size $(size(out)) (or transposed $(size(out_t))), expected $target_size")
end

"""
    _project_bperp_maps(...)

Projects magnetic field onto the sky plane according to LOS.
"""
function _project_bperp_maps(cfg::InstrumentalConfig)
    Bx = read_FITS(cfg.Bx_in)
    By = read_FITS(cfg.By_in)
    Bz = read_FITS(cfg.Bz_in)

    ndims(Bx) == 3 || error("Bx must be 3D, got size=$(size(Bx))")
    ndims(By) == 3 || error("By must be 3D, got size=$(size(By))")
    ndims(Bz) == 3 || error("Bz must be 3D, got size=$(size(Bz))")
    size(Bx) == size(By) || error("Bx/By shape mismatch: Bx=$(size(Bx)) By=$(size(By))")
    size(By) == size(Bz) || error("By/Bz shape mismatch: By=$(size(By)) Bz=$(size(Bz))")

    target = (cfg.n, cfg.m)
    B1_cube, B2_cube = _select_bsky_component_cubes(Bx, By, Bz, cfg.los)

    b1 = _project_cube_mean(B1_cube, cfg.los, target)
    b2 = _project_cube_mean(B2_cube, cfg.los, target)
    bperp = _project_hypot_mean(B1_cube, B2_cube, cfg.los, target)

    return b1, b2, bperp
end
