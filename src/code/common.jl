include(joinpath(@__DIR__, "lib", "DepolLib.jl"))
using .DepolLib

# Backward-compatible aliases for existing scripts.
const DEFAULT_CONFIG_PATH = DepolLib.DEFAULT_CONFIG_PATH

"""
    load_runtime_config(...)

    Alias to `DepolLib.load_runtime_config`.
"""
load_runtime_config(task::String; args::Vector{String}=ARGS) = DepolLib.load_runtime_config(task; args=args)
"""
    build_cfg(...)

    Alias to `DepolLib.build_cfg`.
"""
build_cfg(task::String; config_path::AbstractString=DEFAULT_CONFIG_PATH, overrides::AbstractDict=Dict{String,Any}()) =
    DepolLib.build_cfg(task; config_path=config_path, overrides=overrides)
"""
    cfg_get(...)

    Alias to `DepolLib.cfg_get`.
"""
cfg_get(cfg::AbstractDict, keys::Vector{String}; default=nothing) = DepolLib.cfg_get(cfg, keys; default=default)
"""
    cfg_require(...)

    Alias to `DepolLib.cfg_require`.
"""
cfg_require(cfg::AbstractDict, keys::Vector{String}) = DepolLib.cfg_require(cfg, keys)

"""
    require_los(...)

    Alias to `DepolLib.require_los`.
"""
require_los(los::AbstractString) = DepolLib.require_los(los)
"""
    los_axis(...)

    Alias to `DepolLib.los_axis`.
"""
los_axis(los::AbstractString) = DepolLib.los_axis(los)
"""
    axis_pc(...)

    Alias to `DepolLib.axis_pc`.
"""
axis_pc(n::Int; lbox_pc::Real=50.0) = DepolLib.axis_pc(n; lbox_pc=lbox_pc)
"""
    ticks_pc(...)

    Alias to `DepolLib.ticks_pc`.
"""
ticks_pc(n::Int; lbox_pc::Real=50.0, step_pc::Real=10.0) = DepolLib.ticks_pc(n; lbox_pc=lbox_pc, step_pc=step_pc)

"""
    read_intervals_from_csv(...)

    Alias to `DepolLib.read_intervals_from_csv`.
"""
read_intervals_from_csv(path::String) = DepolLib.read_intervals_from_csv(path)
"""
    merge_intervals(...)

    Alias to `DepolLib.merge_intervals`.
"""
merge_intervals(intervals::Vector{Tuple{Int,Int}}; gap_pix::Int=1) = DepolLib.merge_intervals(intervals; gap_pix=gap_pix)

"""
    normalize_task_name(...)

    Alias to `DepolLib.normalize_task_name`.
"""
normalize_task_name(task::AbstractString) = DepolLib.normalize_task_name(task)
"""
    standard_output_dir(...)

    Alias to `DepolLib.standard_output_dir`.
"""
standard_output_dir(cfg::AbstractDict, task::AbstractString; simu=nothing, los=nothing) = DepolLib.standard_output_dir(cfg, task; simu=simu, los=los)
"""
    standard_output_path(...)

    Alias to `DepolLib.standard_output_path`.
"""
standard_output_path(cfg::AbstractDict, task::AbstractString, artifact::AbstractString, ext::AbstractString; simu=nothing, los=nothing) = DepolLib.standard_output_path(cfg, task, artifact, ext; simu=simu, los=los)

"""
    finite_values(...)

    Alias to `DepolLib.finite_values`.
"""
finite_values(A) = DepolLib.finite_values(A)
"""
    finite_minmax(...)

    Alias to `DepolLib.finite_minmax`.
"""
finite_minmax(A) = DepolLib.finite_minmax(A)
"""
    resolve_simulations_root(...)

    Resolves `paths.simulations_root` with fallback support.
"""
resolve_simulations_root(cfg::AbstractDict) = DepolLib.resolve_simulations_root(cfg)
"""
    require_existing_files(...)

    Verifies all listed files exist.
"""
require_existing_files(paths::AbstractVector{<:AbstractString}; context::AbstractString="required inputs") = DepolLib.require_existing_files(paths; context=context)
"""
    read_fits_f32(...)

    Alias to `DepolLib.read_fits_f32`.
"""
read_fits_f32(path::AbstractString) = DepolLib.read_fits_f32(path)
"""
    require_ndims(...)

    Alias to `DepolLib.require_ndims`.
"""
require_ndims(A, n::Int, label::AbstractString) = DepolLib.require_ndims(A, n, label)
"""
    require_same_size(...)

    Alias to `DepolLib.require_same_size`.
"""
require_same_size(arrays::AbstractVector, labels::AbstractVector{<:AbstractString}) = DepolLib.require_same_size(arrays, labels)

"""
    read_FITS(...)

    Reads the primary HDU from a FITS file.
"""
read_FITS(path::AbstractString) = DepolLib.read_FITS(path)
"""
    write_FITS(...)

    Writes array data to FITS.
"""
write_FITS(path::AbstractString, data; overwrite::Bool=true) = DepolLib.write_FITS(path, data; overwrite=overwrite)
"""
    los_config(...)

    Returns LOS field name and a closure for 1D profile extraction.
"""
los_config(los::String) = DepolLib.los_config(los)
"""
    smooth_moving_average(...)

    Moving-average smoothing with odd window enforcement.
"""
smooth_moving_average(x::AbstractVector, w::Int) = DepolLib.smooth_moving_average(x, w)
"""
    sign_eps(...)

    Robust sign function with dead zone `[-eps,+eps]`.
"""
sign_eps(x::Real; eps::Real=0.0) = DepolLib.sign_eps(x; eps=eps)
"""
    reversal_indices(...)

    Detects sign-reversal indices in a profile.
"""
reversal_indices(B::AbstractVector; eps::Real=0.0) = DepolLib.reversal_indices(B; eps=eps)
