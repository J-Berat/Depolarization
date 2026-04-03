"""
    psd2d(...)

    2D power spectral density.
"""
function psd2d(img::AbstractMatrix)
    A = float.(img)
    F = fft(A) ./ sqrt(float(length(A)))
    P = abs2.(F)
    return fftshift(P)
end

"""
    psd1d_isotropic(...)

    Radially binned isotropic 1D PSD.
"""
function psd1d_isotropic(img::AbstractMatrix, kx::AbstractVector, ky::AbstractVector;
                         nbins::Int=80, kmin::Real=0.0, kmax::Real=Inf)
    n, m = size(img)
    @assert length(kx) == n
    @assert length(ky) == m

    P2 = psd2d(img)
    start = max(float(kmin), 0.0)

    kx2max = 0.0
    for v in kx
        if isfinite(v)
            vv = abs2(float(v))
            vv > kx2max && (kx2max = vv)
        end
    end

    ky2max = 0.0
    for v in ky
        if isfinite(v)
            vv = abs2(float(v))
            vv > ky2max && (ky2max = vv)
        end
    end

    stop = min(float(kmax), sqrt(kx2max + ky2max))
    edges = range(start, stop; length=nbins + 1)
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    Pk = fill(NaN, nbins)

    if !(isfinite(start) && isfinite(stop) && stop > start)
        return centers, Pk
    end

    sums = zeros(Float64, nbins)
    counts = zeros(Int, nbins)
    invw = nbins / (stop - start)

    @inbounds for j in 1:m
        kyj = float(ky[j])
        for i in 1:n
            kxi = float(kx[i])
            p = P2[i, j]
            if !(isfinite(kxi) && isfinite(kyj) && isfinite(p))
                continue
            end
            kval = hypot(kxi, kyj)
            if !(kval >= start && kval <= stop)
                continue
            end
            bin = Int(floor((kval - start) * invw)) + 1
            if 1 <= bin <= nbins
                sums[bin] += p
                counts[bin] += 1
            end
        end
    end

    @inbounds for b in 1:nbins
        if counts[b] > 0
            Pk[b] = sums[b] / counts[b]
        end
    end

    return centers, Pk
end

"""
    psd1d_x_mean_over_y(...)

    1D PSD along `kx`, averaged over `y`.
"""
function psd1d_x_mean_over_y(img::AbstractMatrix, kx::AbstractVector; remove_mean::Bool=false)
    n, _ = size(img)
    @assert length(kx) == n

    A = float.(img)
    if remove_mean
        A .-= mean(A)
    end

    Fx = fft(A, 1) ./ sqrt(float(n))
    _, m = size(Fx)
    Pkx = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        acc = 0.0
        for j in 1:m
            acc += abs2(Fx[i, j])
        end
        Pkx[i] = acc / m
    end
    return kx, fftshift(Pkx)
end

"""
    psd1d_x_stats_over_y(...)

1D PSD along `kx`, with mean and standard deviation across the transverse cuts.
"""
function psd1d_x_stats_over_y(img::AbstractMatrix, kx::AbstractVector; remove_mean::Bool=false)
    n, m = size(img)
    @assert length(kx) == n

    A = float.(img)
    if remove_mean
        A .-= mean(A)
    end

    Fx = fft(A, 1) ./ sqrt(float(n))
    Prow = abs2.(Fx)
    Prow = fftshift(Prow, 1)

    Pmean = Vector{Float64}(undef, n)
    Pstd = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        vals = view(Prow, i, :)
        μ = mean(vals)
        Pmean[i] = μ
        Pstd[i] = std(vals; corrected=false)
    end
    return kx, Pmean, Pstd
end

"""
    _smooth_for_peak(...)

Light smoothing used only for peak detection.
"""
function _smooth_for_peak(y::AbstractVector; radius::Int=2)
    n = length(y)
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        acc = 0.0
        w = 0
        i1 = max(1, i - radius)
        i2 = min(n, i + radius)
        for j in i1:i2
            yj = y[j]
            if isfinite(yj)
                acc += Float64(yj)
                w += 1
            end
        end
        out[i] = w > 0 ? acc / w : NaN
    end
    return out
end

"""
    kpeak_in_window(...)

    Finds dominant spectral peak in a `k` window.
"""
function kpeak_in_window(k::AbstractVector, y::AbstractVector; kmin::Real=0.0, kmax::Real=Inf, smooth_radius::Int=2)
    ys = smooth_radius > 0 ? _smooth_for_peak(y; radius=smooth_radius) : Float64.(y)
    best_k = NaN
    best_y = -Inf
    @inbounds for i in eachindex(k, ys)
        ki = k[i]
        yi = ys[i]
        if isfinite(ki) && isfinite(yi) && ki > 0 && yi > 0 && ki >= kmin && ki <= kmax
            if yi > best_y
                best_y = yi
                best_k = ki
            end
        end
    end
    return best_k
end
