"""
    Shared support for canal detection, smoothing, and connected components.

These helpers are used by the plotting, alignment, and morphology passes so
the canal logic lives in one place.
"""

const _CANAL_LOCAL_RADIUS_PIX = 2
const _CANAL_SCORE_QUANTILE = 0.90
const _CANAL_RIDGE_SMOOTH_RADIUS_PIX = 9
const _CANAL_DERIVATIVE_RADIUS_PIX = 4
const _CANAL_TENSOR_SMOOTH_RADIUS_PIX = 4
const _CANAL_RIDGE_ANIS_MIN = 2.0

function _finite_quantile_or_nan(A, q::Real)
    0.0 <= q <= 1.0 || error("q must be in [0,1], got $q")
    vals = Float64[]
    sizehint!(vals, _count_finite_values(A))
    _append_finite_float_values!(vals, A)
    isempty(vals) && return NaN
    return quantile!(vals, q)
end

@inline _finite_quantile(A, q::Real) = _finite_quantile_or_nan(A, q)

const _ISOTROPIC_OFFSET_CACHE = Dict{Int, Vector{Tuple{Int, Int}}}()

function _disk_offsets(radius::Int)
    radius >= 0 || error("radius must be >= 0, got $radius")
    get!(_ISOTROPIC_OFFSET_CACHE, radius) do
        offs = Tuple{Int, Int}[]
        for dj in -radius:radius, di in -radius:radius
            di^2 + dj^2 <= radius^2 || continue
            push!(offs, (di, dj))
        end
        offs
    end
end

"""
    _box_mean(... )

Backward-compatible name for the isotropic local average used by the canal
method. The kernel support is circular, not square, so the smoothing does not
favor the grid axes.
"""
function _box_mean(A::AbstractMatrix, radius::Int)
    radius >= 0 || error("radius must be >= 0, got $radius")
    radius == 0 && return Float64.(A)

    n, m = size(A)
    offsets = _disk_offsets(radius)
    out = Matrix{Float64}(undef, n, m)

    @inbounds for j in 1:m, i in 1:n
        acc = 0.0
        wsum = 0
        for (di, dj) in offsets
            ii = i + di
            jj = j + dj
            if 1 <= ii <= n && 1 <= jj <= m
                acc += Float64(A[ii, jj])
                wsum += 1
            end
        end
        out[i, j] = acc / max(wsum, 1)
    end
    return out
end

"""
    _canal_score_map(... )

Builds a local-deficit map that highlights dark narrow valleys relative to
their neighbourhood, which follows the visual canal morphology better than
a global intensity threshold.
"""
function _canal_score_map(Pmap::AbstractMatrix; radius_pix::Int=_CANAL_LOCAL_RADIUS_PIX)
    P = Float64.(Pmap)
    valid = isfinite.(P)
    fillval = _finite_quantile(P, 0.50)
    Pfill = copy(P)
    Pfill[.!valid] .= fillval
    Plocal = _box_mean(Pfill, radius_pix)
    score = Plocal .- Pfill
    score[.!valid] .= NaN
    return score
end

"""
    _symmetric_2x2_eigensystem(... )

Returns ordered eigenvalues and eigenvectors for a symmetric `2x2` matrix
`[a b; b c]`. The first eigenpair corresponds to the smallest eigenvalue.
"""
function _symmetric_2x2_eigensystem(a::Real, b::Real, c::Real)
    tr = float(a) + float(c)
    disc = sqrt(max((float(a) - float(c))^2 + 4.0 * float(b)^2, 0.0))
    λ1 = 0.5 * (tr - disc)
    λ2 = 0.5 * (tr + disc)

    # Angle of the eigenvector associated with the largest eigenvalue.
    θ2 = 0.5 * atan(2.0 * float(b), float(a) - float(c))
    v2x = cos(θ2)
    v2y = sin(θ2)
    v1x = -sin(θ2)
    v1y = cos(θ2)

    return λ1, λ2, v1x, v1y, v2x, v2y
end

@inline _wrap_orientation_pi(θ::Real) = mod(float(θ), π)

@inline function _sample_bilinear(A::AbstractMatrix, x::Real, y::Real)
    n, m = size(A)
    if !(1.0 <= x <= n && 1.0 <= y <= m)
        return NaN
    end

    x1 = clamp(floor(Int, x), 1, n)
    y1 = clamp(floor(Int, y), 1, m)
    x2 = clamp(x1 + 1, 1, n)
    y2 = clamp(y1 + 1, 1, m)

    tx = clamp(float(x) - x1, 0.0, 1.0)
    ty = clamp(float(y) - y1, 0.0, 1.0)

    v11 = Float64(A[x1, y1])
    v21 = Float64(A[x2, y1])
    v12 = Float64(A[x1, y2])
    v22 = Float64(A[x2, y2])
    if !(isfinite(v11) && isfinite(v21) && isfinite(v12) && isfinite(v22))
        return NaN
    end

    return (1 - tx) * (1 - ty) * v11 +
           tx * (1 - ty) * v21 +
           (1 - tx) * ty * v12 +
           tx * ty * v22
end

function _gaussian_derivative_kernel(Δx::Real, Δy::Real; radius::Int=_CANAL_DERIVATIVE_RADIUS_PIX)
    radius >= 1 || error("radius must be >= 1, got $radius")
    isfinite(Δx) && Δx > 0 || error("Δx must be positive finite, got $Δx")
    isfinite(Δy) && Δy > 0 || error("Δy must be positive finite, got $Δy")

    σ = max(radius / 2, 1.0) * sqrt(float(Δx) * float(Δy))
    norm = 1.0 / (2π * σ^2)
    Kx = Matrix{Float64}(undef, 2radius + 1, 2radius + 1)
    Ky = Matrix{Float64}(undef, 2radius + 1, 2radius + 1)

    @inbounds for (jj, dj) in enumerate(-radius:radius), (ii, di) in enumerate(-radius:radius)
        x = di * float(Δx)
        y = dj * float(Δy)
        g = norm * exp(-0.5 * (x^2 + y^2) / σ^2)
        Kx[ii, jj] = -(x / σ^2) * g
        Ky[ii, jj] = -(y / σ^2) * g
    end
    return Kx, Ky
end

function _convolve_clamped(A::AbstractMatrix, K::AbstractMatrix)
    n, m = size(A)
    rkx, rky = size(K)
    radius_i = (rkx - 1) ÷ 2
    radius_j = (rky - 1) ÷ 2
    out = Matrix{Float64}(undef, n, m)

    @inbounds for j in 1:m, i in 1:n
        acc = 0.0
        for kj in 1:rky, ki in 1:rkx
            ii = clamp(i + (ki - radius_i - 1), 1, n)
            jj = clamp(j + (kj - radius_j - 1), 1, m)
            acc += Float64(A[ii, jj]) * K[ki, kj]
        end
        out[i, j] = acc
    end
    return out
end

"""
    _gradients_gaussian(...)

First derivatives obtained by convolving the map with the derivatives of an
isotropic Gaussian kernel.
"""
function _gradients_gaussian(A::AbstractMatrix, Δx::Real, Δy::Real; radius::Int=_CANAL_DERIVATIVE_RADIUS_PIX)
    Kx, Ky = _gaussian_derivative_kernel(Δx, Δy; radius=radius)
    gx = _convolve_clamped(A, Kx)
    gy = _convolve_clamped(A, Ky)
    return gx, gy
end

"""
    _gradients_central(...)

Compatibility helper for tests and quick diagnostics based on centered finite
differences with edge clamping.
"""
function _gradients_central(A::AbstractMatrix, Δx::Real, Δy::Real)
    isfinite(Δx) && Δx > 0 || error("Δx must be positive finite, got $Δx")
    isfinite(Δy) && Δy > 0 || error("Δy must be positive finite, got $Δy")

    n, m = size(A)
    gx = fill(NaN, n, m)
    gy = fill(NaN, n, m)

    @inbounds for j in 1:m, i in 1:n
        im = max(i - 1, 1)
        ip = min(i + 1, n)
        jm = max(j - 1, 1)
        jp = min(j + 1, m)

        axm = Float64(A[im, j])
        axp = Float64(A[ip, j])
        aym = Float64(A[i, jm])
        ayp = Float64(A[i, jp])

        if isfinite(axm) && isfinite(axp)
            gx[i, j] = (axp - axm) / ((ip - im) * float(Δx))
        end
        if isfinite(aym) && isfinite(ayp)
            gy[i, j] = (ayp - aym) / ((jp - jm) * float(Δy))
        end
    end

    return gx, gy
end

function _zone_connected_components(mask::AbstractMatrix{Bool})
    nx, ny = size(mask)
    labels = zeros(Int, nx, ny)
    comps = Vector{Vector{CartesianIndex{2}}}()
    neighbors = (
        CartesianIndex(-1, -1), CartesianIndex(-1, 0), CartesianIndex(-1, 1),
        CartesianIndex(0, -1),                         CartesianIndex(0, 1),
        CartesianIndex(1, -1),  CartesianIndex(1, 0), CartesianIndex(1, 1),
    )

    next_label = 0
    for j in 1:ny, i in 1:nx
        if !mask[i, j] || labels[i, j] != 0
            continue
        end
        next_label += 1
        queue = CartesianIndex{2}[CartesianIndex(i, j)]
        labels[i, j] = next_label
        comp = CartesianIndex{2}[]
        qhead = 1
        while qhead <= length(queue)
            idx = queue[qhead]
            qhead += 1
            push!(comp, idx)
            for δ in neighbors
                nb = idx + δ
                ii, jj = Tuple(nb)
                if 1 <= ii <= nx && 1 <= jj <= ny && mask[ii, jj] && labels[ii, jj] == 0
                    labels[ii, jj] = next_label
                    push!(queue, nb)
                end
            end
        end
        push!(comps, comp)
    end
    return comps
end
