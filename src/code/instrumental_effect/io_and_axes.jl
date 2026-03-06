function grid_scales(cfg::InstrumentalConfig)
    Δx = cfg.Lbox_pc / cfg.n
    Δy = cfg.Lbox_pc / cfg.m
    fNy_x = 1 / (2Δx)
    fNy_y = 1 / (2Δy)
    fNy = min(fNy_x, fNy_y)
    return (; Δx, Δy, fNy)
end

function spectral_axes(cfg::InstrumentalConfig, Δx::Real, Δy::Real)
    kx = (2π) .* fftshift(fftfreq(cfg.n, Δx))
    ky = (2π) .* fftshift(fftfreq(cfg.m, Δy))
    x_pc = (0:cfg.n-1) .* Δx
    y_pc = (0:cfg.m-1) .* Δy
    return (; kx, ky, x_pc, y_pc)
end

Pphi_max_map(absF::AbstractArray) =
    dropdims(maximum(absF; dims=ndims(absF)); dims=ndims(absF))

function get_chan_xy(cube, ichan::Int, n::Int, m::Int)
    @assert ndims(cube) == 3 "cube must be 3D. Got ndims=$(ndims(cube)) size=$(size(cube))"
    sz = size(cube)

    if sz[1] == n && sz[2] == m
        @assert 1 <= ichan <= sz[3] "ichan=$(ichan) out of bounds for third dim=$(sz[3])"
        return float.(@view cube[:, :, ichan])
    elseif sz[2] == n && sz[3] == m
        @assert 1 <= ichan <= sz[1] "ichan=$(ichan) out of bounds for first dim=$(sz[1])"
        return float.(reshape(@view(cube[ichan, :, :]), n, m))
    else
        error("Unsupported cube shape $(sz). Expected (n,m,nν) or (nν,n,m) with n=$n m=$m.")
    end
end

function interpol2d(map, xmap, ymap)
    nx, ny = size(map)
    x = range(1, nx, length=nx)
    y = range(1, ny, length=ny)
    itp = linear_interpolation((x, y), map, extrapolation_bc=Line())
    mapout = zeros(nx, ny)
    for p in 1:length(map)
        mapout[p] = itp(xmap[p], ymap[p])
    end
    return mapout
end

function xymap(nx, ny)
    xmap = collect(1:nx) .* ones(Int64, ny)'
    ymap = ones(Int64, nx) .* collect(1:ny)'
    return xmap, ymap
end

function lic(vx, vy; niter=1, len=8, normalize=false, amplitude=nothing, level=0.1, scalar=nothing)
    nx, ny = size(vx)

    uu = @. sqrt(vx^2 + vy^2)
    ii = findall(uu .== 0.0)
    if !isempty(ii)
        uu[ii] .= 1.0
    end

    if normalize
        ux = vx ./ uu
        uy = vy ./ uu
    else
        ux = vx ./ maximum(uu)
        uy = vy ./ maximum(uu)
    end

    vl = rand(nx, ny)

    for _ in 1:niter
        texture = copy(vl)

        vv = zeros(nx, ny)
        pi, pj = xymap(nx, ny)
        mi = copy(pi)
        mj = copy(pj)

        ppi = 1.0 * pi
        ppj = 1.0 * pj
        mmi = 1.0 * mi
        mmj = 1.0 * mj

        for _ in 0:len
            dpi = interpol2d(ux, ppi, ppj)
            dpj = interpol2d(uy, ppi, ppj)
            dmi = interpol2d(ux, mmi, mmj)
            dmj = interpol2d(uy, mmi, mmj)

            ppi = @. ppi + 0.25 * dpi
            ppj = @. ppj + 0.25 * dpj
            mmi = @. mmi - 0.25 * dmi
            mmj = @. mmj - 0.25 * dmj

            pi = @. mod(round(ppi) + nx - 1, nx) + 1
            pj = @. mod(round(ppj) + ny - 1, ny) + 1
            mi = @. mod(round(mmi) + nx - 1, nx) + 1
            mj = @. mod(round(mmj) + ny - 1, ny) + 1

            ppi = @. pi + (ppi - round(ppi))
            ppj = @. pj + (ppj - round(ppj))
            mmi = @. mi + (mmi - round(mmi))
            mmj = @. mj + (mmj - round(mmj))

            tp = interpol2d(texture, ppi, ppj)
            tm = interpol2d(texture, mmi, mmj)
            vv = @. vv + tp + tm
        end

        vl = 0.25 * vv / len
    end

    if amplitude !== nothing
        if scalar !== nothing && length(scalar) > 1
            uu .= scalar
        end

        if length(level) != 1
            level = 0.1
        end

        if !isempty(ii)
            uu[ii] .= 0.0
        end

        level = max([0.0, min([level, 1.0])])
        uu .= (1.0 - level) * uu ./ maximum(uu) .+ level
        vl .= vl .* uu
    end

    if !isempty(ii)
        vl[ii] .= 0.0
    end

    return vl
end

Llabel_pc(Llarge) = round(50 / 256 * Llarge)
