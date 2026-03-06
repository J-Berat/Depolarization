#!/usr/bin/env julia

using Pkg

project_root = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(project_root)

if Base.JLOptions().startupfile != 0
    @warn "For deterministic runs independent of personal Julia startup config, prefer: julia --startup-file=no --project=. scripts/precompile.jl"
end

required_pkgs = [
    "CairoMakie",
    "FITSIO",
    "LaTeXStrings",
    "StatsBase",
    "Makie",
    "FFTW",
    "Interpolations",
]

proj = Pkg.project()
existing = Set(keys(proj.dependencies))
missing = [p for p in required_pkgs if !(p in existing)]

if !isempty(missing)
    println("[0/3] Adding missing packages to local project: ", join(missing, ", "))
    Pkg.add(missing)
end

println("[1/3] Instantiating project environment...")
Pkg.instantiate()

println("[2/3] Precompiling package dependencies...")
Pkg.precompile()

println("[3/3] Warming up local modules/jobs...")

# Load local code paths (guarded scripts won't execute workloads).
include(joinpath(project_root, "src", "code", "lib", "DepolLib.jl"))
include(joinpath(project_root, "src", "code", "instrumental_effect", "InstrumentalEffect.jl"))

include(joinpath(project_root, "src", "code", "jobs", "run_mach_suite.jl"))
include(joinpath(project_root, "src", "code", "jobs", "run_reversal_transition_job.jl"))
include(joinpath(project_root, "src", "code", "jobs", "run_segmentation_pipeline_job.jl"))
include(joinpath(project_root, "src", "code", "jobs", "run_instrumental_effect_job.jl"))

job_fns = (
    :run_mach_suite,
    :run_reversal_transition_job,
    :run_segmentation_pipeline_job,
    :run_instrumental_effect_job,
)

for fn in job_fns
    if isdefined(Main, fn)
        precompile(getproperty(Main, fn), (Dict{String,Any},))
    end
end

if isdefined(Main, :InstrumentalEffect)
    IE = Main.InstrumentalEffect

    if isdefined(IE, :RMSynthesis)
        precompile(IE.RMSynthesis, (Array{Float64,3}, Array{Float64,3}, Vector{Float64}, Vector{Float64}))
    end
    if isdefined(IE, :getRMSF)
        precompile(IE.getRMSF, (Vector{Float64}, Vector{Float64}))
    end
    if isdefined(IE, :run_filter_pass) && isdefined(IE, :InstrumentalConfig)
        precompile(IE.run_filter_pass, (IE.InstrumentalConfig,))
    end
    if isdefined(IE, :run_pipeline) && isdefined(IE, :InstrumentalConfig)
        precompile(IE.run_pipeline, (IE.InstrumentalConfig,))
    end

    # Tiny synthetic warmup for RM kernels.
    Q = randn(8, 8, 8)
    U = randn(8, 8, 8)
    nu = collect(range(9.5e8, 1.15e9; length=8))
    phi = collect(range(-10.0, 10.0; length=33))
    try
        IE.RMSynthesis(Q, U, nu, phi; log_progress=false)
        IE.getRMSF(nu, phi; log_progress=false)
    catch err
        @warn "RM synthetic warmup skipped" error=string(err)
    end
end

println("Precompile complete.")
