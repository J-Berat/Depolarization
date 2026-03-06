module DepolLib

using TOML
using Statistics
using FITSIO

include(joinpath(@__DIR__, "..", "fits_io.jl"))
include(joinpath(@__DIR__, "..", "los_utils.jl"))

export DEFAULT_CONFIG_PATH,
       load_runtime_config, build_cfg, cfg_get, cfg_require,
       require_los, los_axis, axis_pc, ticks_pc,
       read_FITS, write_FITS, los_config, smooth_moving_average, sign_eps, reversal_indices,
       read_fits_f32, require_ndims, require_same_size,
       finite_values, finite_minmax,
       resolve_simulations_root, require_existing_files,
       task_enabled, skipped_job_result, run_job_entrypoint,
       simulation_dir, simulation_field_path, synchrotron_dir, synchrotron_path, withfaraday_dir, withfaraday_path,
       read_intervals_from_csv, merge_intervals,
       normalize_task_name, standard_output_dir, standard_output_path,
       run_script_with_same_config, default_simulation_list

const DEFAULT_CONFIG_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "default.toml"))

"""
    _parse_scalar(...)

    Converts a string to a scalar type (`Bool`, `Nothing`, `Int`, `Float64`, otherwise `String`).
"""
function _parse_scalar(value::AbstractString)
    value = String(value)
    low = lowercase(value)
    low == "true" && return true
    low == "false" && return false
    low == "nothing" && return nothing

    try
        return parse(Int, value)
    catch
    end

    try
        return parse(Float64, value)
    catch
    end

    return value
end

"""
    _set_nested!(...)

    Writes a value into nested config dictionaries using a dotted key like `a.b.c`.
"""
function _set_nested!(cfg::Dict{String,Any}, dotted_key::String, value)
    parts = split(dotted_key, '.')
    isempty(parts) && error("Invalid key in --set: '$dotted_key'")

    node = cfg
    for p in parts[1:end-1]
        if !haskey(node, p) || !(node[p] isa AbstractDict)
            node[p] = Dict{String,Any}()
        end
        node = node[p]
    end
    node[parts[end]] = value
end

"""
    _parse_cli_args(...)

    Parses `--config` and `--set` CLI arguments.
"""
function _parse_cli_args(args::Vector{String})
    config_path = DEFAULT_CONFIG_PATH
    set_args = String[]

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i += 1
            i <= length(args) || error("Missing value after --config")
            config_path = args[i]
        elseif arg == "--set"
            i += 1
            i <= length(args) || error("Missing value after --set")
            push!(set_args, args[i])
        else
            error("Unknown CLI argument '$arg'. Allowed: --config, --set")
        end
        i += 1
    end

    return config_path, set_args
end

"""
    _read_config_dict(...)

    Reads TOML config and returns `(cfg, absolute_config_path)`.
"""
function _read_config_dict(config_path::AbstractString)
    path = String(config_path)
    isfile(path) || error("Config file not found: $path")
    cfg = Dict{String,Any}(TOML.parsefile(path))
    return cfg, abspath(path)
end

"""
    _set_runtime_metadata!(...)

    Stores runtime metadata fields under `cfg[\"runtime\"]`.
"""
function _set_runtime_metadata!(cfg::Dict{String,Any}, task::String, config_path::AbstractString, set_args::Vector{String})
    cfg["runtime"] = get(cfg, "runtime", Dict{String,Any}())
    cfg["runtime"]["task"] = task
    cfg["runtime"]["config_path"] = abspath(String(config_path))
    cfg["runtime"]["set_args"] = copy(set_args)
    return cfg
end

"""
    _apply_cli_set_arg!(...)

    Applies one CLI `--set key=value` override.
"""
function _apply_cli_set_arg!(cfg::Dict{String,Any}, kv::AbstractString)
    occursin("=", kv) || error("Invalid --set value '$kv'. Expected key=value")
    k, v = split(kv, "="; limit=2)
    _set_nested!(cfg, String(k), _parse_scalar(v))
    return true
end

"""
    load_runtime_config(...)

    Loads runtime TOML config, applies CLI overrides, and stores runtime metadata.
"""
function load_runtime_config(task::String; args::Vector{String}=ARGS)
    config_path, set_args = _parse_cli_args(args)
    cfg, abs_config_path = _read_config_dict(config_path)

    for kv in set_args
        _apply_cli_set_arg!(cfg, kv)
    end

    return _set_runtime_metadata!(cfg, task, abs_config_path, set_args)
end

"""
    build_cfg(...)

    Loads runtime TOML config and applies dictionary overrides using dotted keys.
    This provides the same runtime metadata contract as `load_runtime_config`,
    without requiring CLI-style `args`.
"""
function build_cfg(task::String; config_path::AbstractString=DEFAULT_CONFIG_PATH, overrides::AbstractDict=Dict{String,Any}())
    cfg, abs_config_path = _read_config_dict(config_path)
    set_args = String[]

    for (k_raw, v_raw) in pairs(overrides)
        key = if k_raw isa Symbol
            String(k_raw)
        elseif k_raw isa AbstractString
            String(k_raw)
        else
            error("Override keys must be String or Symbol, got $(typeof(k_raw))")
        end
        isempty(strip(key)) && error("Override key cannot be empty")

        value = v_raw isa AbstractString ? _parse_scalar(v_raw) : v_raw
        _set_nested!(cfg, key, value)
        push!(set_args, string(key, "=", v_raw))
    end

    return _set_runtime_metadata!(cfg, task, abs_config_path, set_args)
end

"""
    cfg_get(...)

    Reads a nested config value with a default fallback.
"""
function cfg_get(cfg::AbstractDict, keys::Vector{String}; default=nothing)
    node = cfg
    for k in keys
        if !(node isa AbstractDict) || !haskey(node, k)
            return default
        end
        node = node[k]
    end
    return node
end

"""
    cfg_require(...)

    Reads a nested config value and errors if missing.
"""
function cfg_require(cfg::AbstractDict, keys::Vector{String})
    val = cfg_get(cfg, keys; default=nothing)
    val === nothing && error("Missing required config key: $(join(keys, '.'))")
    return val
end

"""
    require_los(...)

    Validates LOS (`x|y|z`) and returns normalized LOS.
"""
function require_los(los::AbstractString)
    (los == "x" || los == "y" || los == "z") || error("LOS must be one of: x, y, z")
    return String(los)
end

"""
    los_axis(...)

    Maps LOS to numeric axis (`1|2|3`).
"""
function los_axis(los::AbstractString)
    los = require_los(los)
    return los == "x" ? 1 : los == "y" ? 2 : 3
end

"""
    axis_pc(...)

    Builds a linear axis in parsecs.
"""
axis_pc(n::Int; lbox_pc::Real=50.0) = range(0.0, lbox_pc; length=n)

"""
    ticks_pc(...)

    Returns tick pixel positions and parsec labels.
"""
function ticks_pc(n::Int; lbox_pc::Real=50.0, step_pc::Real=10.0)
    pcs = collect(0:step_pc:lbox_pc)
    pos = round.(Int, 1 .+ (n - 1) .* (pcs ./ lbox_pc))
    return pos, string.(Int.(pcs))
end

"""
    read_fits_f32(...)

    Reads FITS data and casts to `Float32`.
"""
read_fits_f32(path::AbstractString) = Float32.(read_FITS(path))

"""
    require_ndims(...)

    Asserts a specific number of dimensions.
"""
function require_ndims(A, n::Int, label::AbstractString)
    ndims(A) == n || error("$label must be $n-D, got ndims=$(ndims(A)) size=$(size(A))")
    return A
end

"""
    require_same_size(...)

    Asserts equal sizes across multiple arrays.
"""
function require_same_size(arrays::AbstractVector, labels::AbstractVector{<:AbstractString})
    isempty(arrays) && return true
    s0 = size(arrays[1])
    for i in 2:length(arrays)
        size(arrays[i]) == s0 || error("Size mismatch: $(labels[1])=$s0, $(labels[i])=$(size(arrays[i]))")
    end
    return true
end

"""
    finite_values(...)

    Extracts finite values from an array.
"""
function finite_values(A)
    v = vec(A)
    return v[isfinite.(v)]
end

"""
    finite_minmax(...)

    Returns finite min/max.
"""
function finite_minmax(A)
    vals = finite_values(A)
    isempty(vals) && error("No finite values found")
    return minimum(vals), maximum(vals)
end

"""
    resolve_simulations_root(...)

    Resolves `paths.simulations_root` with a practical local fallback to
    `~/Desktop/simu_RAMSES` when the configured path is missing.
"""
function resolve_simulations_root(cfg::AbstractDict)
    configured = string(cfg_require(cfg, ["paths", "simulations_root"]))
    cand1 = abspath(configured)
    if isdir(cand1)
        return cand1
    end

    fallback = joinpath(homedir(), "Desktop", "simu_RAMSES")
    if isdir(fallback)
        @warn "Configured simulations_root not found, using fallback" configured=cand1 fallback=fallback
        return fallback
    end

    error("simulations_root not found: '$cand1' (and fallback '$fallback' missing)")
end

"""
    require_existing_files(...)

    Ensures all provided files exist, otherwise raises one explicit error.
"""
function require_existing_files(paths::AbstractVector{<:AbstractString}; context::AbstractString="required inputs")
    missing = [String(p) for p in paths if !isfile(p)]
    isempty(missing) && return true
    msg = "Missing files for $context:\\n  " * join(missing, "\\n  ")
    error(msg)
end

"""
    task_enabled(...)

    Returns a boolean task flag from config with a default.
"""
task_enabled(cfg::AbstractDict, keys::Vector{String}; default::Bool=true) = Bool(cfg_get(cfg, keys; default=default))

"""
    skipped_job_result(...)

    Standard payload for disabled jobs.
"""
function skipped_job_result(task::AbstractString, reason::AbstractString)
    return Dict(
        "task" => String(task),
        "status" => "skipped",
        "reason" => String(reason),
    )
end

"""
    run_job_entrypoint(...)

    Standard CLI entrypoint wrapper for jobs.
"""
function run_job_entrypoint(task::AbstractString, run_job::Function; args::Vector{String}=ARGS)
    cfg = load_runtime_config(String(task); args=args)
    result = run_job(cfg)
    println(result)
    return result
end

"""
    simulation_dir(...)

    Returns `<simulations_root>/<simulation_name>`.
"""
simulation_dir(sim_root::AbstractString, sim::AbstractString) = joinpath(String(sim_root), String(sim))

"""
    simulation_field_path(...)

    Returns `<simulations_root>/<simulation_name>/<field>.fits`.
"""
function simulation_field_path(sim_root::AbstractString, sim::AbstractString, field::AbstractString)
    return joinpath(simulation_dir(sim_root, sim), string(field, ".fits"))
end

"""
    synchrotron_dir(...)

    Returns `<simulations_root>/<simulation_name>/<los>/Synchrotron`.
"""
function synchrotron_dir(sim_root::AbstractString, sim::AbstractString, los::AbstractString)
    return joinpath(simulation_dir(sim_root, sim), require_los(los), "Synchrotron")
end

"""
    synchrotron_path(...)

    Returns `<simulations_root>/<simulation_name>/<los>/Synchrotron/<filename>`.
"""
function synchrotron_path(sim_root::AbstractString, sim::AbstractString, los::AbstractString, filename::AbstractString)
    return joinpath(synchrotron_dir(sim_root, sim, los), String(filename))
end

"""
    withfaraday_dir(...)

    Returns `<simulations_root>/<simulation_name>/<los>/Synchrotron/WithFaraday`.
"""
function withfaraday_dir(sim_root::AbstractString, sim::AbstractString, los::AbstractString)
    return joinpath(synchrotron_dir(sim_root, sim, los), "WithFaraday")
end

"""
    withfaraday_path(...)

    Returns `<simulations_root>/<simulation_name>/<los>/Synchrotron/WithFaraday/<filename>`.
"""
function withfaraday_path(sim_root::AbstractString, sim::AbstractString, los::AbstractString, filename::AbstractString)
    return joinpath(withfaraday_dir(sim_root, sim, los), String(filename))
end

"""
    read_intervals_from_csv(...)

    Reads `kmin,kmax` intervals from CSV.
"""
function read_intervals_from_csv(path::String)
    isfile(path) || error("Missing CSV: $path")
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: $path")

    header = split(strip(lines[1]), ',')
    col = Dict(h => i for (i, h) in enumerate(header))
    (haskey(col, "kmin") && haskey(col, "kmax")) || error("CSV must contain kmin,kmax")

    intervals = Tuple{Int,Int}[]
    for ln in lines[2:end]
        s = strip(ln)
        isempty(s) && continue
        parts = split(s, ',')
        push!(intervals, (parse(Int, parts[col["kmin"]]), parse(Int, parts[col["kmax"]])))
    end

    sort!(intervals, by=x -> x[1])
    return intervals
end

"""
    merge_intervals(...)

    Merges overlapping or near-adjacent intervals.
"""
function merge_intervals(intervals::Vector{Tuple{Int,Int}}; gap_pix::Int=1)
    isempty(intervals) && return Tuple{Int,Int}[]
    out = Tuple{Int,Int}[]
    kL, kR = intervals[1]
    for (L, R) in intervals[2:end]
        if L <= kR + gap_pix
            kR = max(kR, R)
        else
            push!(out, (kL, kR))
            kL, kR = L, R
        end
    end
    push!(out, (kL, kR))
    return out
end

"""
    normalize_task_name(...)

    Normalizes task names by replacing whitespace with underscores.
"""
normalize_task_name(task::AbstractString) = replace(task, r"\s+" => "_")

"""
    standard_output_dir(...)

    Creates/returns a standard output directory.
"""
function standard_output_dir(cfg::AbstractDict, task::AbstractString; simu=nothing, los=nothing)
    out_root = string(cfg_require(cfg, ["paths", "outputs_root"]))
    simu_name = simu === nothing ? string(cfg_require(cfg, ["simulation", "name"])) : string(simu)
    los_name = los === nothing ? string(cfg_require(cfg, ["simulation", "los"])) : string(los)
    require_los(los_name)

    dir = joinpath(out_root, normalize_task_name(task), simu_name, "LOS$(los_name)")
    mkpath(dir)
    return dir
end

"""
    standard_output_path(...)

    Builds a standard output file path.
"""
function standard_output_path(cfg::AbstractDict, task::AbstractString, artifact::AbstractString, ext::AbstractString; simu=nothing, los=nothing)
    simu_name = simu === nothing ? string(cfg_require(cfg, ["simulation", "name"])) : string(simu)
    los_name = los === nothing ? string(cfg_require(cfg, ["simulation", "los"])) : string(los)
    dir = standard_output_dir(cfg, task; simu=simu_name, los=los_name)
    fname = string(normalize_task_name(task), "_", simu_name, "_LOS", los_name, "_", artifact, ".", ext)
    return joinpath(dir, fname)
end

"""
    run_script_with_same_config(...)

    Runs another Julia script with the same runtime config.
"""
function run_script_with_same_config(cfg::AbstractDict, script_path::AbstractString)
    config_path = string(cfg_require(cfg, ["runtime", "config_path"]))
    set_args = cfg_get(cfg, ["runtime", "set_args"]; default=String[])

    julia_exe = get(ENV, "JULIA_EXE", Base.julia_cmd().exec[1])
    cmd_parts = String[julia_exe]
    push!(cmd_parts, "--startup-file=no")
    active_project = Base.active_project()
    if active_project !== nothing
        push!(cmd_parts, "--project=$(dirname(active_project))")
    end
    push!(cmd_parts, script_path, "--config", config_path)
    for kv in set_args
        push!(cmd_parts, "--set")
        push!(cmd_parts, string(kv))
    end

    run(Cmd(cmd_parts))
    return true
end

"""
    default_simulation_list(...)

    Returns the project default simulation list.
"""
function default_simulation_list()
    return [
        "d1cf00bx10rms18000nograv1024",
        "d1cf02bx10rms09000nograv1024",
        "d1cf02bx10rms18000nograv1024",
        "d1cf02bx10rms36000nograv1024",
        "d1cf02bx10rms72000nograv1024",
        "d1cf05bx10rms09000nograv1024",
        "d1cf05bx10rms18000nograv1024",
        "d1cf05bx10rms36000nograv1024",
        "d1cf05bx10rms72000nograv1024",
        "d1cf08bx10rms09000nograv1024",
        "d1cf08bx10rms18000nograv1024",
        "d1cf08bx10rms36000nograv1024",
        "d1cf08bx10rms72000nograv1024",
        "d1cf10bx10rms18000nograv1024",
    ]
end

end
