"""
    grid_scales(...)

    Computes grid spacing and Nyquist limits.
"""
function grid_scales(cfg::InstrumentalConfig)
    cfg.n > 0 || error("cfg.n must be > 0, got $(cfg.n)")
    cfg.m > 0 || error("cfg.m must be > 0, got $(cfg.m)")
    cfg.Lbox_pc > 0 || error("cfg.Lbox_pc must be > 0, got $(cfg.Lbox_pc)")
    Δx = cfg.Lbox_pc / cfg.n
    Δy = cfg.Lbox_pc / cfg.m
    fNy_x = 1 / (2Δx)
    fNy_y = 1 / (2Δy)
    fNy = min(fNy_x, fNy_y)
    return (; Δx, Δy, fNy)
end

"""
    spectral_axes(...)

    Builds `kx`, `ky`, `x_pc`, `y_pc` axes.
"""
function spectral_axes(cfg::InstrumentalConfig, Δx::Real, Δy::Real)
    isfinite(Δx) && Δx > 0 || error("Δx must be positive finite, got $Δx")
    isfinite(Δy) && Δy > 0 || error("Δy must be positive finite, got $Δy")
    kx = (2π) .* fftshift(fftfreq(cfg.n, Δx))
    ky = (2π) .* fftshift(fftfreq(cfg.m, Δy))
    x_pc = (0:cfg.n-1) .* Δx
    y_pc = (0:cfg.m-1) .* Δy
    return (; kx, ky, x_pc, y_pc)
end

"""
    Pphi_max_map(...)

    Collapses over `phi` by maximum to get `Pphi_max` map.
"""
Pphi_max_map(absF::AbstractArray) =
    dropdims(maximum(absF; dims=ndims(absF)); dims=ndims(absF))

"""
    get_chan_xy(...)

    Extracts a 2D channel from a 3D cube with layout detection.
"""
function get_chan_xy(cube, ichan::Int, n::Int, m::Int)
    ndims(cube) == 3 || error("cube must be 3D. Got ndims=$(ndims(cube)) size=$(size(cube))")
    sz = size(cube)

    if sz[1] == n && sz[2] == m
        1 <= ichan <= sz[3] || error("ichan=$ichan out of bounds for third dimension=$(sz[3])")
        return float.(@view cube[:, :, ichan])
    elseif sz[2] == n && sz[3] == m
        1 <= ichan <= sz[1] || error("ichan=$ichan out of bounds for first dimension=$(sz[1])")
        return float.(reshape(@view(cube[ichan, :, :]), n, m))
    else
        error("Unsupported cube shape $(sz). Expected (n,m,nν) or (nν,n,m) with n=$n m=$m.")
    end
end

"""
    interpol2d(...)

    Bilinear interpolation on target index maps.
"""
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

"""
    xymap(...)

    Generates `x`/`y` index grids.
"""
function xymap(nx, ny)
    xmap = collect(1:nx) .* ones(Int64, ny)'
    ymap = ones(Int64, nx) .* collect(1:ny)'
    return xmap, ymap
end

"""
    lic(...)

    Computes line-integral-convolution texture from a vector field.
"""
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

"""
    Llabel_pc(...)

    Converts filter scale label to parsec units.
"""
Llabel_pc(Llarge) = round(50 / 256 * Llarge)

const C_m = 299_792_458.0

"""
    print_progress(progress::Int, total::Int)

    Prints a dynamic ANSI progress bar.
"""
function print_progress(progress::Int, total::Int)
    @assert total > 0 "Total must be greater than 0."
    @assert 0 <= progress <= total "Progress must be between 0 and total."

    bar_width = 50
    progress_ratio = progress / total
    filled_length = Int(round(bar_width * progress_ratio))
    empty_length = bar_width - filled_length

    green = "\u001b[42m"
    reset = "\u001b[0m"
    gray = "\u001b[47m"

    filled_bar = repeat(" ", filled_length) |> x -> green * x * reset
    empty_bar = repeat(" ", empty_length) |> x -> gray * x * reset

    percentage = Int(round(100 * progress_ratio))

    print("\rProgress: |$filled_bar$empty_bar| $progress/$total ($percentage%)")
    flush(stdout)

    if progress == total
        println()
    end
end

"""
    RMSynthesis(Q::AbstractArray, U::AbstractArray, nuArray::AbstractArray, PhiArray::AbstractArray) -> Tuple

Perform RM synthesis on Stokes `Q` and `U`.

Supported inputs:
- 1D: `(nν,)`
- 2D: `(nx, nν)` (interpreted as a single-row map)
- 3D: `(nx, ny, nν)` or `(nν, nx, ny)` (auto-detected by `length(nuArray)`).
"""
function RMSynthesis(Q::AbstractArray, U::AbstractArray, nuArray::AbstractArray, PhiArray::AbstractArray; log_progress::Bool=false)
    length(nuArray) > 0 || error("nuArray must be non-empty")
    length(PhiArray) > 0 || error("PhiArray must be non-empty")

    log_progress && @info "Starting RM synthesis" n_phi=length(PhiArray) n_lambda=length(nuArray)

    LambdaSqArray = @. (C_m / Float64.(nuArray))^2
    nPhi = length(PhiArray)
    nLambda = length(LambdaSqArray)
    nDims = ndims(Q)

    size(Q) == size(U) || error("Q/U size mismatch: size(Q)=$(size(Q)) size(U)=$(size(U))")

    WeightArray = ones(Float64, nLambda)
    K = 1.0 / sum(WeightArray)

    # Normalize Q/U to (nx, ny, nν)
    if nDims == 1
        length(Q) == nLambda || error("1D Q/U must have length(nuArray) elements")
        Q = reshape(Float64.(Q), (1, 1, length(Q)))
        U = reshape(Float64.(U), (1, 1, length(U)))
    elseif nDims == 2
        size(Q, 2) == nLambda || error("2D Q/U expected shape (nx, nν) with nν=$(nLambda), got size(Q)=$(size(Q))")
        Q = reshape(Float64.(Q), (1, size(Q, 1), size(Q, 2)))
        U = reshape(Float64.(U), (1, size(U, 1), size(U, 2)))
    elseif nDims == 3
        if size(Q, 3) == nLambda
            Q = Float64.(Q)
            U = Float64.(U)
        elseif size(Q, 1) == nLambda
            Q = permutedims(Float64.(Q), (2, 3, 1))
            U = permutedims(Float64.(U), (2, 3, 1))
        else
            error("3D Q/U must have frequency axis first or last matching nν=$(nLambda), got size(Q)=$(size(Q))")
        end
    else
        error("Unsupported Q/U dimensions: ndims=$(nDims)")
    end

    P = @. (Q + 1im * U) * WeightArray[[CartesianIndex()], [CartesianIndex()], :]
    nx, ny = size(Q, 1), size(Q, 2)

    Lambda0Sq = K * sum(WeightArray .* LambdaSqArray)
    a = LambdaSqArray .- Lambda0Sq

    F = complex(zeros(Float64, nx, ny, nPhi))
    for i in 1:nPhi
        arg = exp.((-2.0im * Float64(PhiArray[i])) .* a)[[CartesianIndex()], [CartesianIndex()], :]
        F[:, :, i] = K .* sum(P .* arg, dims=3)
        if log_progress
            print_progress(i, nPhi)
            @debug "RM synthesis accumulation" idx=i total=nPhi
        end
    end

    if nDims == 1
        F = dropdims(dropdims(F, dims=1), dims=1)
    elseif nDims == 2
        F = dropdims(F, dims=1)
    end

    log_progress && @info "RM synthesis complete" output_size=size(F)
    return abs.(F), real.(F), imag.(F)
end

"""
    getRMSF(nuArray::AbstractArray, PhiArray::AbstractArray) -> Tuple{AbstractArray, Float64}

Compute RMSF amplitude and FWHM.
"""
function getRMSF(nuArray::AbstractArray, PhiArray::AbstractArray; log_progress::Bool=false)
    length(nuArray) > 0 || error("nuArray must be non-empty")
    length(PhiArray) > 0 || error("PhiArray must be non-empty")

    log_progress && @info "Starting RMSF computation" n_phi=length(PhiArray)

    LambdaSqArray = @. (C_m / Float64.(nuArray))^2
    nPhi = length(PhiArray)

    WeightArray = ones(Float64, length(LambdaSqArray))
    K = 1.0 / sum(WeightArray)

    Lambda0Sq = K * sum(WeightArray .* LambdaSqArray)
    a = LambdaSqArray .- Lambda0Sq

    fwhmRMSF = 3.8 / (maximum(LambdaSqArray) - minimum(LambdaSqArray))

    RMSF = complex(zeros(Float64, nPhi))
    for i in 1:nPhi
        arg = exp.((-2.0im * Float64(PhiArray[i])) .* a)
        RMSF[i] = K * sum(WeightArray .* arg)
        if log_progress
            print_progress(i, nPhi)
            @debug "RMSF accumulation" idx=i total=nPhi
        end
    end

    log_progress && @info "RMSF computation complete" output_size=length(RMSF) fwhm=fwhmRMSF
    return abs.(RMSF), fwhmRMSF
end
