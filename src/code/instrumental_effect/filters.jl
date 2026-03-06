function apply_instrument_2d(img::AbstractMatrix, H::AbstractMatrix)
    real.(ifft(fft(img) .* H))
end

function apply_to_array_xy(data, H; n::Int=256, m::Int=256)
    nd = ndims(data)

    if nd == 2
        return apply_instrument_2d(data, H)
    elseif nd == 3
        sz = size(data)
        out = similar(float.(data))

        if sz[1] == n && sz[2] == m
            @views for k in 1:sz[3]
                out[:, :, k] = apply_instrument_2d(data[:, :, k], H)
            end
            return out
        elseif sz[2] == n && sz[3] == m
            @views for k in 1:sz[1]
                tmp = apply_instrument_2d(reshape(data[k, :, :], n, m), H)
                out[k, :, :] .= reshape(tmp, 1, n, m)
            end
            return out
        else
            error("Unsupported 3D shape $(sz)")
        end
    else
        error("Unsupported ndims(data)=$nd")
    end
end

function instrument_bandpass_L(n::Int, m::Int;
                               Δx::Real, Δy::Real=Δx,
                               Lcut_small::Real,
                               Llarge::Real,
                               fNy::Real)
    fx = fftfreq(n, Δx)
    fy = fftfreq(m, Δy)
    FX = repeat(fx, 1, m)
    FY = repeat(fy', n, 1)
    F = sqrt.(FX.^2 .+ FY.^2)

    flo = 1 / Llarge
    fhi_raw = 1 / Lcut_small
    fhi = min(fhi_raw, fNy)

    @info "Filter" Lcut_small Llarge flo fhi fhi_raw

    H = Float32.((F .>= flo) .& (F .<= fhi))
    return H, fftshift(H)
end
