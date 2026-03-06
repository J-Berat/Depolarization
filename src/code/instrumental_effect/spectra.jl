function psd2d(img::AbstractMatrix)
    A = float.(img)
    A .-= mean(A)
    F = fft(A)
    P = abs2.(F) ./ (length(A)^2)
    return fftshift(P)
end

function psd1d_isotropic(img::AbstractMatrix, kx::AbstractVector, ky::AbstractVector;
                         nbins::Int=80, kmin::Real=0.0, kmax::Real=Inf)
    n, m = size(img)
    @assert length(kx) == n
    @assert length(ky) == m

    P2 = psd2d(img)

    KX = repeat(kx, 1, m)
    KY = repeat(ky', n, 1)
    K = sqrt.(KX.^2 .+ KY.^2)

    kvec = vec(K)
    pvec = vec(P2)

    good = isfinite.(kvec) .& isfinite.(pvec) .& (kvec .>= kmin) .& (kvec .<= kmax)
    kvec = kvec[good]
    pvec = pvec[good]

    kmax_eff = isempty(kvec) ? 0.0 : maximum(kvec)
    edges = range(max(kmin, 0.0), min(kmax, kmax_eff), length=nbins + 1)
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])

    Pk = fill(NaN, nbins)
    @inbounds for i in 1:nbins
        lo, hi = edges[i], edges[i + 1]
        sel = (kvec .>= lo) .& (kvec .< hi)
        if any(sel)
            Pk[i] = mean(pvec[sel])
        end
    end

    return centers, Pk
end

function psd1d_x_mean_over_y(img::AbstractMatrix, kx::AbstractVector; remove_mean::Bool=true)
    n, _ = size(img)
    @assert length(kx) == n

    A = float.(img)
    if remove_mean
        A .-= mean(A)
    end

    Fx = fft(A, 1)
    Px = abs2.(Fx) ./ (n^2)
    Pkx = vec(mean(Px; dims=2))
    return kx, fftshift(Pkx)
end

function kpeak_in_window(k::AbstractVector, y::AbstractVector; kmin::Real=0.0, kmax::Real=Inf)
    ok = isfinite.(k) .& isfinite.(y) .& (k .> 0) .& (y .> 0) .& (k .>= kmin) .& (k .<= kmax)
    if !any(ok)
        return NaN
    end
    kk = k[ok]
    yy = y[ok]
    return kk[argmax(yy)]
end
