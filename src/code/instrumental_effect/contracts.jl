"""
    FilterPassResult

Typed payload exchanged between compute and plotting stages.
"""
struct FilterPassResult{
    TQ<:AbstractArray,
    TU<:AbstractArray,
    TP<:AbstractMatrix,
    TS<:NamedTuple,
    TA<:NamedTuple,
}
    Qdata::TQ
    Udata::TU
    Qslice_filt::Dict{Float64, Matrix{Float32}}
    Uslice_filt::Dict{Float64, Matrix{Float32}}
    Pmax0::TP
    Pmax_filt::Dict{Float64, Matrix{Float32}}
    L_ok::Vector{Float64}
    scales::TS
    axes::TA
    nuArray::Vector{Float64}
    PhiArray::Vector{Float64}
end

"""
    validate_instrumental_config!(...)

Validates basic instrumental configuration constraints.
"""
function validate_instrumental_config!(cfg::InstrumentalConfig)
    cfg.n > 0 || error("cfg.n must be > 0, got $(cfg.n)")
    cfg.m > 0 || error("cfg.m must be > 0, got $(cfg.m)")
    cfg.Lbox_pc > 0 || error("cfg.Lbox_pc must be > 0, got $(cfg.Lbox_pc)")
    cfg.Lcut_small > 0 || error("cfg.Lcut_small must be > 0, got $(cfg.Lcut_small)")
    isempty(cfg.Llarge_list) && error("cfg.Llarge_list must not be empty")
    all(isfinite, cfg.Llarge_list) || error("cfg.Llarge_list must contain finite values")
    all(>(0), cfg.Llarge_list) || error("cfg.Llarge_list must contain positive values")
    cfg.Δν_MHz > 0 || error("cfg.Δν_MHz must be > 0, got $(cfg.Δν_MHz)")
    cfg.νmax_MHz >= cfg.νmin_MHz || error("cfg.νmax_MHz must be >= cfg.νmin_MHz")
    cfg.ichan >= 1 || error("cfg.ichan must be >= 1, got $(cfg.ichan)")
    cfg.iphi >= 1 || error("cfg.iphi must be >= 1, got $(cfg.iphi)")
    return true
end

"""
    require_channel_index(...)

Checks channel index bounds.
"""
function require_channel_index(idx::Int, nmax::Int, label::AbstractString)
    1 <= idx <= nmax || error("$label index=$idx out of bounds 1:$nmax")
    return true
end

"""
    require_stokes_cube_layout(...)

Validates a Stokes cube against configured `(n,m,nν)` or `(nν,n,m)` layout.
Returns frequency-axis index (`1` or `3`).
"""
function require_stokes_cube_layout(cube, label::AbstractString, cfg::InstrumentalConfig, nν::Int)
    ndims(cube) == 3 || error("$label must be 3-D, got ndims=$(ndims(cube)) size=$(size(cube))")
    sz = size(cube)
    if sz[1] == cfg.n && sz[2] == cfg.m
        sz[3] == nν || error("$label frequency axis mismatch: expected nν=$nν on dim=3, got size=$sz")
        return 3
    elseif sz[2] == cfg.n && sz[3] == cfg.m
        sz[1] == nν || error("$label frequency axis mismatch: expected nν=$nν on dim=1, got size=$sz")
        return 1
    else
        error("$label must match (n,m,nν) or (nν,n,m) with n=$(cfg.n) m=$(cfg.m), got size=$sz")
    end
end

"""
    require_map_layout(...)

Validates a `(n,m)` map.
"""
function require_map_layout(map, label::AbstractString, cfg::InstrumentalConfig)
    ndims(map) == 2 || error("$label must be 2-D, got ndims=$(ndims(map)) size=$(size(map))")
    size(map) == (cfg.n, cfg.m) || error("$label must have size ($(cfg.n),$(cfg.m)), got $(size(map))")
    return true
end
