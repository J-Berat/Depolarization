using Statistics

function los_config(los::String)
    if los == "x"
        return ("Bx", (i, j, B) -> Array(@view B[:, i, j]))
    elseif los == "y"
        return ("By", (i, j, B) -> Array(@view B[i, :, j]))
    elseif los == "z"
        return ("Bz", (i, j, B) -> Array(@view B[i, j, :]))
    else
        error("LOS must be x/y/z")
    end
end

function smooth_moving_average(x::AbstractVector, w::Int)
    w = max(w, 1)
    w = isodd(w) ? w : (w + 1)
    n = length(x)
    if w == 1
        return collect(x)
    end
    y = similar(x, n)
    h = (w - 1) ÷ 2
    @inbounds for i in 1:n
        i1 = max(1, i - h)
        i2 = min(n, i + h)
        y[i] = mean(@view x[i1:i2])
    end
    return y
end

@inline function sign_eps(x::Real; eps::Real=0.0)
    if x > eps
        return 1
    elseif x < -eps
        return -1
    else
        return 0
    end
end

function reversal_indices(B::AbstractVector; eps::Real=0.0)
    idx = Int[]
    @inbounds for k in 1:(length(B)-1)
        s1 = sign_eps(B[k]; eps=eps)
        s2 = sign_eps(B[k+1]; eps=eps)
        if s1 != 0 && s2 != 0 && s1 != s2
            push!(idx, k)
        end
    end
    return idx
end
