module DepolLib

using TOML
using Statistics
using FITSIO

include(joinpath(@__DIR__, "..", "fits_io.jl"))
include(joinpath(@__DIR__, "..", "los_utils.jl"))

export DEFAULT_CONFIG_PATH,
       load_runtime_config, build_cfg, cfg_get, cfg_require,
       depol_help, repl_help,
       require_los, los_axis, axis_pc, ticks_pc,
       read_FITS, write_FITS, los_config, smooth_moving_average, sign_eps, reversal_indices,
       read_fits_f32, require_ndims, require_same_size,
       Wolfire_ne, dm_em_maps,
       finite_values, finite_minmax,
       resolve_simulations_root, require_existing_files,
       task_enabled, skipped_job_result, run_job_entrypoint,
       simulation_dir, simulation_field_path, synchrotron_dir, synchrotron_path, withfaraday_dir, withfaraday_path,
       read_intervals_from_csv, merge_intervals,
       normalize_task_name, standard_output_dir, standard_output_path,
       run_script_with_same_config, default_simulation_list

const DEFAULT_CONFIG_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "default.toml"))
const OUTPUTS_ROOT_DEFAULT = joinpath(homedir(), "Desktop", "depolarization_outputs")
const _ANNOUNCED_OUTPUT_DIRS = Set{String}()
const _ANNOUNCED_OUTPUT_PATHS = Set{String}()
const _ANNOUNCED_LEGACY_OUTPUT_ROOT = Ref(false)

const _HELP_JOB_ORDER = [
    "ne_dm_em",
    "mach",
    "reversal",
    "canal_metrics",
    "segmentation",
    "instrumental",
]

const _HELP_TOPIC_ALIASES = Dict(
    "canal" => "canal_metrics",
    "rev" => "reversal",
    "seg" => "segmentation",
    "inst" => "instrumental",
    "overview" => "all",
    "full" => "all",
    "alljobs" => "jobs",
    "pipeline" => "order",
)

const _HELP_JOB_SPECS = Dict(
    "ne_dm_em" => (
        summary="Compute electron density (ne) and integrated DM/EM maps.",
        include_path="src/code/jobs/run_ne_dm_em_job.jl",
        task_name="ne_dm_em",
        run_func="run_ne_dm_em_job",
        config_sections=("tasks.ne_dm_em",),
    ),
    "mach" => (
        summary="Compute <Ms>/<MA> suite metrics and write scatter plots + summary CSV.",
        include_path="src/code/jobs/run_mach_suite.jl",
        task_name="mach",
        run_func="run_mach_suite",
        config_sections=("tasks.mach_suite",),
    ),
    "reversal" => (
        summary="Detect B_LOS reversals, build transition windows, and analyze FDF transitions.",
        include_path="src/code/jobs/run_reversal_transition_job.jl",
        task_name="reversal",
        run_func="run_reversal_transition_job",
        config_sections=("tasks.reversal_transition_job", "tasks.reversals_map", "tasks.b_transition"),
    ),
    "canal_metrics" => (
        summary="Compute Cb/Cphi coherence maps and canal-vs-non-canal diagnostics.",
        include_path="src/code/jobs/run_canal_metrics_job.jl",
        task_name="canal_metrics",
        run_func="run_canal_metrics_job",
        config_sections=("tasks.canal_metrics",),
    ),
    "segmentation" => (
        summary="Chunk cubes, cut transition intervals, and produce Pmax panel outputs.",
        include_path="src/code/jobs/run_segmentation_pipeline_job.jl",
        task_name="segmentation",
        run_func="run_segmentation_pipeline_job",
        config_sections=("tasks.segmentation_pipeline_job", "tasks.cut_transition", "tasks.make_chunk"),
    ),
    "instrumental" => (
        summary="Run instrumental filtering/analysis on Q/U/Pmax/FDF products.",
        include_path="src/code/jobs/run_instrumental_effect_job.jl",
        task_name="instrumental",
        run_func="run_instrumental_effect_job",
        config_sections=("tasks.instrumental_effect",),
    ),
)

function _resolve_help_topic(topic::AbstractString)
    topic_norm = lowercase(strip(topic))
    isempty(topic_norm) && return "all"
    return get(_HELP_TOPIC_ALIASES, topic_norm, topic_norm)
end

function _print_help_setup(io::IO)
    print(io, """
1) Start Julia in project mode
   julia --startup-file=no --project=.

2) Load core library
   include("src/code/lib/DepolLib.jl")
   using .DepolLib

3) Build config (shared pattern for all jobs)
   cfg = build_cfg("instrumental";
       config_path="config/default.toml",
       overrides=Dict(
         "paths.simulations_root" => "/path/to/simu_RAMSES",
         "paths.desktop_output_root" => "./outputs",
         "simulation.name" => "d1cf05bx10rms18000nograv1024",
         "simulation.los" => "y",
       ));

""")
    return nothing
end

function _print_help_job(io::IO, job::String)
    spec = _HELP_JOB_SPECS[job]
    println(io, "   ", job)
    println(io, "     ", spec.summary)
    println(io, "     config: ", join(spec.config_sections, ", "))
    println(io, "     include(\"", spec.include_path, "\")")
    println(io, "     cfg = build_cfg(\"", spec.task_name, "\"; config_path=\"config/default.toml\", overrides=Dict(...))")
    println(io, "     result = ", spec.run_func, "(cfg)")
    println(io)
    return nothing
end

function _print_help_jobs(io::IO; only::Union{Nothing,String}=nothing)
    println(io, "4) Jobs: what they do + how to run in REPL")
    println(io)
    jobs = only === nothing ? _HELP_JOB_ORDER : [only]
    for job in jobs
        _print_help_job(io, job)
    end
    return nothing
end

function _print_help_order(io::IO)
    println(io, "5) Recommended order (fresh simulation)")
    for (i, job) in enumerate(_HELP_JOB_ORDER)
        println(io, "   ", i, ". ", job)
    end
    println(io)
    return nothing
end

function _print_help_cli_overrides(io::IO)
    print(io, """6) CLI-style overrides also work in REPL:
   cfg = load_runtime_config("reversal"; args=[
       "--config", "config/default.toml",
       "--set", "paths.simulations_root=/path/to/simu_RAMSES",
       "--set", "simulation.name=d1cf05bx10rms18000nograv1024",
       "--set", "simulation.los=y",
   ]);
""")
    return nothing
end

function _print_help_cli_commands(io::IO)
    println(io, "7) CLI cheatsheet (one command per job)")
    for job in _HELP_JOB_ORDER
        spec = _HELP_JOB_SPECS[job]
        println(io, "   julia --startup-file=no --project=. ", spec.include_path, " --config config/default.toml")
    end
    println(io)
    return nothing
end

function _print_help_run_all(io::IO)
    print(io, """8) Run all jobs (CLI)
   for job in \\
     src/code/jobs/run_ne_dm_em_job.jl \\
     src/code/jobs/run_mach_suite.jl \\
     src/code/jobs/run_reversal_transition_job.jl \\
     src/code/jobs/run_canal_metrics_job.jl \\
     src/code/jobs/run_segmentation_pipeline_job.jl \\
     src/code/jobs/run_instrumental_effect_job.jl
   do
     julia --startup-file=no --project=. "\$job" --config config/default.toml
   done
""")
    return nothing
end

function _known_help_topics()
    aliases = sort(collect(keys(_HELP_TOPIC_ALIASES)))
    return join(vcat(["all", "jobs", "order", "cli"], _HELP_JOB_ORDER, aliases), ", ")
end

"""
    help([io]; topic="all")
    help(topic)

Prints REPL usage guidance.

`topic` can be:
- `"all"`: setup + jobs + order + CLI examples
- `"jobs"`: setup + all jobs
- `"order"`: recommended execution order
- `"cli"`: CLI examples (single-job and run-all)
- one job name (`"ne_dm_em"`, `"mach"`, `"reversal"`, `"canal_metrics"`, `"segmentation"`, `"instrumental"`)
- aliases (`"canal"`, `"rev"`, `"seg"`, `"inst"`)
"""
function help(io::IO=stdout; topic::AbstractString="all")
    topic_norm = _resolve_help_topic(topic)

    print(io, """Depolarization REPL quick help

Usage:
  depol_help()                     # full guide (preferred, no name conflict)
  depol_help("jobs")               # all jobs
  depol_help("canal")              # alias for canal_metrics
  depol_help("instrumental")       # one job
  depol_help("order")              # recommended order
  depol_help("cli")                # CLI examples
  DepolLib.help("instrumental")    # qualified fallback

""")

    if topic_norm == "all"
        _print_help_setup(io)
        _print_help_jobs(io)
        _print_help_order(io)
        _print_help_cli_overrides(io)
        _print_help_cli_commands(io)
        _print_help_run_all(io)
    elseif topic_norm == "jobs"
        _print_help_setup(io)
        _print_help_jobs(io)
    elseif topic_norm == "order"
        _print_help_order(io)
    elseif topic_norm == "cli"
        _print_help_cli_overrides(io)
        _print_help_cli_commands(io)
        _print_help_run_all(io)
    elseif haskey(_HELP_JOB_SPECS, topic_norm)
        _print_help_setup(io)
        _print_help_jobs(io; only=topic_norm)
    else
        error("Unknown help topic '$topic'. Known topics: $(_known_help_topics())")
    end

    return nothing
end

help(topic::AbstractString) = help(stdout; topic=topic)

"""
    depol_help(...)

Alias for `help(...)` with a non-ambiguous exported name.
"""
depol_help(args...; kwargs...) = help(args...; kwargs...)

"""
    repl_help(...)

Alias for `depol_help(...)`.
"""
repl_help(args...; kwargs...) = depol_help(args...; kwargs...)

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
    Wolfire_ne(...)

    Computes electron density with the Wolfire-style approximation:

    `ne = 2.4e-3*sqrt(zeta/1e-16)*(T/100)^0.25*sqrt(Geff)/omegaPAH + n*XC`.
"""
function Wolfire_ne(zeta::Real, Geff::Real, omegaPAH::Real, XC::Real, T::AbstractArray, n::AbstractArray)
    zeta > 0 || error("zeta must be > 0, got $zeta")
    Geff > 0 || error("Geff must be > 0, got $Geff")
    omegaPAH > 0 || error("omegaPAH must be > 0, got $omegaPAH")
    require_same_size([T, n], ["T", "n"])

    return @. 2.4e-3 * sqrt(zeta / 1e-16) * (T / 100)^0.25 * sqrt(Geff) / omegaPAH + n * XC
end

"""
    dm_em_maps(...)

    Integrates a 3D `ne` cube along LOS to build:
    - `DM = ∫ ne dl` [pc cm^-3]
    - `EM = ∫ ne^2 dl` [pc cm^-6]

    Returns `(dm_map, em_map, dl_pc)`.
"""
function dm_em_maps(ne::AbstractArray{<:Real,3}, los::AbstractString; lbox_pc::Real=50.0)
    ax = los_axis(los)
    nlos = size(ne, ax)
    nlos > 0 || error("LOS size must be positive, got $nlos")
    # Cell-centered discretization: for N pixels over Lbox, use dl = Lbox / N.
    dl_pc = Float64(lbox_pc) / nlos

    dm = dropdims(sum(ne; dims=ax); dims=ax) .* dl_pc
    em = dropdims(sum(abs2, ne; dims=ax); dims=ax) .* dl_pc
    return dm, em, dl_pc
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

    Normalizes/simplifies task names for output folders and filenames.
"""
function normalize_task_name(task::AbstractString)
    name = lowercase(strip(String(task)))
    name = replace(name, r"\s+" => "_")
    name = replace(name, r"^run_" => "")
    name = replace(name, r"_job(?=_|$)" => "")
    name = replace(name, r"_suite(?=_|$)" => "")
    name = replace(name, r"_pipeline(?=_|$)" => "")
    name = replace(name, r"_+" => "_")
    return strip(name, '_')
end

"""
    _resolve_output_root(...)

Resolves output root from config keys with backward-compatible fallback:
- `paths.desktop_output_root` (preferred)
- `paths.outputs_root` (legacy)
- default `~/Desktop/depolarization_outputs`
"""
function _resolve_output_root(cfg::AbstractDict)
    preferred = cfg_get(cfg, ["paths", "desktop_output_root"]; default=nothing)
    if preferred !== nothing
        preferred_str = strip(string(preferred))
        !isempty(preferred_str) && return preferred_str
    end

    legacy = cfg_get(cfg, ["paths", "outputs_root"]; default=nothing)
    if legacy !== nothing
        legacy_str = strip(string(legacy))
        if !isempty(legacy_str)
            if !_ANNOUNCED_LEGACY_OUTPUT_ROOT[]
                @warn "paths.outputs_root is deprecated; prefer paths.desktop_output_root" outputs_root=legacy_str
                _ANNOUNCED_LEGACY_OUTPUT_ROOT[] = true
            end
            return legacy_str
        end
    end

    return OUTPUTS_ROOT_DEFAULT
end

"""
    standard_output_dir(...)

Creates/returns a standard output directory.
"""
function standard_output_dir(cfg::AbstractDict, task::AbstractString; simu=nothing, los=nothing)
    out_root = _resolve_output_root(cfg)
    simu_name = simu === nothing ? string(cfg_require(cfg, ["simulation", "name"])) : string(simu)
    los_name = los === nothing ? string(cfg_require(cfg, ["simulation", "los"])) : string(los)
    require_los(los_name)

    dir = joinpath(out_root, normalize_task_name(task), simu_name, "LOS$(los_name)")
    mkpath(dir)
    if !(dir in _ANNOUNCED_OUTPUT_DIRS)
        push!(_ANNOUNCED_OUTPUT_DIRS, dir)
        @info "Output directory ready" task=normalize_task_name(task) out_dir=dir
    end
    return dir
end

"""
    standard_output_path(...)

    Builds a standard output file path with a short filename (`artifact.ext`).
"""
function standard_output_path(cfg::AbstractDict, task::AbstractString, artifact::AbstractString, ext::AbstractString; simu=nothing, los=nothing)
    simu_name = simu === nothing ? string(cfg_require(cfg, ["simulation", "name"])) : string(simu)
    los_name = los === nothing ? string(cfg_require(cfg, ["simulation", "los"])) : string(los)
    dir = standard_output_dir(cfg, task; simu=simu_name, los=los_name)
    fname = string(replace(String(artifact), r"\s+" => "_"), ".", ext)
    out = joinpath(dir, fname)
    if !(out in _ANNOUNCED_OUTPUT_PATHS)
        push!(_ANNOUNCED_OUTPUT_PATHS, out)
        @info "Output path prepared" task=normalize_task_name(task) file=out
    end
    return out
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
