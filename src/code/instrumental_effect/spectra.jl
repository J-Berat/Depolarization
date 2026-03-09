"""
    psd2d(...)

    2D power spectral density.
"""
function psd2d(img::AbstractMatrix)
    A = float.(img)
    A .-= mean(A)
    F = fft(A)
    P = abs2.(F) ./ (length(A)^2)
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
function psd1d_x_mean_over_y(img::AbstractMatrix, kx::AbstractVector; remove_mean::Bool=true)
    n, _ = size(img)
    @assert length(kx) == n

    A = float.(img)
    if remove_mean
        A .-= mean(A)
    end

    Fx = fft(A, 1)
    invn2 = inv(float(n)^2)
    _, m = size(Fx)
    Pkx = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        acc = 0.0
        for j in 1:m
            acc += abs2(Fx[i, j])
        end
        Pkx[i] = (acc / m) * invn2
    end
    return kx, fftshift(Pkx)
end

"""
    kpeak_in_window(...)

    Finds dominant spectral peak in a `k` window.
"""
function kpeak_in_window(k::AbstractVector, y::AbstractVector; kmin::Real=0.0, kmax::Real=Inf)
    best_k = NaN
    best_y = -Inf
    @inbounds for i in eachindex(k, y)
        ki = k[i]
        yi = y[i]
        if isfinite(ki) && isfinite(yi) && ki > 0 && yi > 0 && ki >= kmin && ki <= kmax
            if yi > best_y
                best_y = yi
                best_k = ki
            end
        end
    end
    return best_k
end
