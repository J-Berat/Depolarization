# Getting Started

This guide is for first-time users who want to run an existing workflow with local simulation data.

## Prerequisites

- Julia (tested in CI with Julia 1.12)
- Local clone of this repository
- Access to simulation FITS data under `paths.simulations_root`

## Initial Setup

From repository root:

```bash
julia --startup-file=no --project=. scripts/precompile.jl
```

This command:

- activates the local project
- installs missing dependencies
- precompiles dependencies and project code

If you are new to the repository structure, read `repository_layout.md` before going deeper into `scripts/` or `src/code/`.

## First Job Run (CLI)

Example with explicit runtime overrides:

```bash
julia --startup-file=no --project=. src/code/jobs/run_instrumental_effect_job.jl \
  --config config/default.toml \
  --set paths.simulations_root=/path/to/simu_RAMSES \
  --set paths.desktop_output_root=./outputs \
  --set simulation.name=d1cf05bx10rms18000nograv1024 \
  --set simulation.los=y
```

All jobs use the same CLI contract:

- `--config <path/to/file.toml>`
- repeatable `--set key=value`

Only these two flags are accepted by runtime config parsing.

The most important runtime values to override for a new machine are:

- `paths.simulations_root`
- `paths.desktop_output_root`
- `simulation.name`
- `simulation.los`

## Precompile + Run Wrapper

You can force a precompile before each run:

```bash
./scripts/run_with_precompile.sh reversal_transition --config config/default.toml
```

Accepted job values include:

- full path, e.g. `src/code/jobs/run_reversal_transition_job.jl`
- filename, e.g. `run_reversal_transition_job`
- task alias, e.g. `reversal_transition`

## REPL Workflow

```bash
julia --startup-file=no --project=.
```

Then:

```julia
julia> include("src/code/lib/DepolLib.jl")
julia> using .DepolLib

julia> cfg = build_cfg("instrumental";
           config_path="config/default.toml",
           overrides=Dict(
             "paths.simulations_root" => "/path/to/simu_RAMSES",
             "paths.desktop_output_root" => "./outputs",
             "simulation.name" => "d1cf05bx10rms18000nograv1024",
             "simulation.los" => "y",
           ))

julia> include("src/code/jobs/run_instrumental_effect_job.jl")
julia> result = run_instrumental_effect_job(cfg)
```

Built-in help is available in REPL:

```julia
julia> depol_help()
julia> depol_help("jobs")
julia> depol_help("instrumental")
```

## Running All Jobs

```bash
for job in \
  src/code/jobs/run_ne_dm_em_job.jl \
  src/code/jobs/run_mach_suite.jl \
  src/code/jobs/run_reversal_transition_job.jl \
  src/code/jobs/run_canal_metrics_job.jl \
  src/code/jobs/run_segmentation_pipeline_job.jl \
  src/code/jobs/run_instrumental_effect_job.jl

do
  julia --startup-file=no --project=. "$job" --config config/default.toml
done
```

For a first external read, it is usually better to run one job first and inspect the resulting output tree before launching the full sequence.

## Output Root Behavior

Output root selection order:

1. `paths.desktop_output_root` (preferred)
2. `paths.outputs_root` (legacy fallback)
3. `~/Desktop/depolarization_outputs` (default fallback)

## Troubleshooting

- Missing files: verify `paths.simulations_root`, `simulation.name`, and required FITS files.
- Invalid LOS: only `x`, `y`, `z` are allowed.
- CLI parsing errors: remove unsupported flags and keep only `--config`/`--set`.
- Unsure where a result should land: check `docs/jobs.md` for the standard output layout of each job.
