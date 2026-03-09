"""
    apply_instrument_2d(...)

    Applies 2D Fourier-domain instrument filter.
"""
function apply_instrument_2d(img::AbstractMatrix, H::AbstractMatrix)
    size(img) == size(H) || error("Filter shape mismatch: image size=$(size(img)) filter size=$(size(H))")
    real.(ifft(fft(img) .* H))
end

"""
    apply_to_array_xy(...)

    Applies filter to 2D/3D arrays with supported axis layouts.
"""
function apply_to_array_xy(data, H; n::Int=256, m::Int=256)
    size(H) == (n, m) || error("Filter mask H must have size ($n,$m), got $(size(H))")
    nd = ndims(data)

    if nd == 2
        size(data) == (n, m) || error("2D input must have size ($n,$m), got $(size(data))")
        return apply_instrument_2d(data, H)
    elseif nd == 3
        sz = size(data)
        Tout = float(eltype(data))
        out = similar(data, Tout, sz)

        if sz[1] == n && sz[2] == m
            @views for k in 1:sz[3]
                out[:, :, k] = apply_instrument_2d(data[:, :, k], H)
            end
            return out
        elseif sz[2] == n && sz[3] == m
            @views for k in 1:sz[1]
                out[k, :, :] = apply_instrument_2d(data[k, :, :], H)
            end
            return out
        else
            error("Unsupported 3D shape $(sz). Expected (n,m,nν) or (nν,n,m) with n=$n m=$m.")
        end
    else
        error("Unsupported ndims(data)=$nd")
    end
end

"""
    instrument_bandpass_L(...)

    Builds spatial band-pass mask.
"""
function instrument_bandpass_L(n::Int, m::Int;
                               Δx::Real, Δy::Real=Δx,
                               Lcut_small::Real,
                               Llarge::Real,
                               fNy::Real)
    fx = fftfreq(n, Δx)
    fy = fftfreq(m, Δy)

    flo = 1 / Llarge
    fhi_raw = 1 / Lcut_small
    fhi = min(fhi_raw, fNy)

    @debug "Filter" Lcut_small Llarge flo fhi fhi_raw

    flo2 = flo^2
    fhi2 = fhi^2
    H = Matrix{Float32}(undef, n, m)
    @inbounds for j in 1:m
        fy2 = float(fy[j])^2
        for i in 1:n
            f2 = float(fx[i])^2 + fy2
            H[i, j] = (f2 >= flo2 && f2 <= fhi2) ? 1f0 : 0f0
        end
    end
    return H, fftshift(H)
end
