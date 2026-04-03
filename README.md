# Depolarization

Julia workflows for depolarization, synchrotron, and Faraday-analysis products derived from RAMSES simulation cubes.

This repository is organized around reproducible batch jobs: each job reads FITS inputs from a simulation tree, applies one scientific workflow, and writes structured outputs under a task-specific directory.

## Who This Repo Is For

This repository is easiest to approach if you are:

- running the existing analysis pipeline on local RAMSES simulation data
- inspecting how a given figure or product is generated
- extending one of the existing Julia jobs

If you do not already have the simulation FITS files referenced by `config/default.toml`, you can still read the code and docs, but you will not be able to execute the full workflows end to end.

## Start Here

```bash
# 1) Install dependencies and precompile once
julia --startup-file=no --project=. scripts/precompile.jl

# 2) Run one job
julia --startup-file=no --project=. src/code/jobs/run_instrumental_effect_job.jl \
  --config config/default.toml \
  --set paths.simulations_root=/path/to/simu_RAMSES \
  --set paths.desktop_output_root=./outputs \
  --set simulation.name=d1cf05bx10rms18000nograv1024 \
  --set simulation.los=y

# 3) Inspect outputs
ls -R ./outputs/instrumental/d1cf05bx10rms18000nograv1024/LOSy
```

All executable jobs follow the same CLI contract:

```bash
julia --startup-file=no --project=. <job>.jl --config path/to/file.toml --set key=value
```

Only two runtime flags are accepted:

- `--config <path>`
- repeatable `--set key=value`

## Repository Map

- `src/code/jobs/`: executable job entrypoints.
- `src/code/lib/DepolLib.jl`: shared configuration, path, validation, and REPL helpers.
- `src/code/instrumental_effect/`: internals for the instrumental-effect pipeline.
- `config/`: default runtime configuration.
- `docs/`: user-facing documentation and developer notes.
- `scripts/`: helper scripts for setup, documentation generation, figure generation, and ad hoc analysis.
- `test/`: regression and integration-style tests.
- `outputs/`: generated analysis products for local runs; ignored by Git.

A fuller explanation of what belongs where is available in [docs/repository_layout.md](docs/repository_layout.md).

## Main Jobs

Recommended execution order for a fresh simulation:

1. `run_ne_dm_em_job.jl`
2. `run_mach_suite.jl`
3. `run_reversal_transition_job.jl`
4. `run_canal_metrics_job.jl`
5. `run_segmentation_pipeline_job.jl`
6. `run_instrumental_effect_job.jl`

Job reference:

- [docs/jobs.md](docs/jobs.md): purpose, inputs, outputs, and config section for each job.
- [docs/instrumental_effect.md](docs/instrumental_effect.md): detailed contract for the instrumental-effect pipeline.

## Documentation Guide

- [docs/README.md](docs/README.md): documentation index with reading paths.
- [docs/getting_started.md](docs/getting_started.md): setup, first execution, CLI and REPL usage.
- [docs/configuration.md](docs/configuration.md): configuration keys, overrides, and task sections.
- [docs/jobs.md](docs/jobs.md): job-by-job inputs, outputs, and execution order.
- [docs/repository_layout.md](docs/repository_layout.md): structure of the repository and what is generated locally.
- [docs/development.md](docs/development.md): architecture, tests, and maintenance workflow.
- [docs/functions.md](docs/functions.md): auto-generated function index for `src/code/`.

Regenerate function docs with:

```bash
./scripts/generate_functions_doc.sh
```

## Expected Input Data

Most workflows expect a simulation tree shaped like:

```text
<simulations_root>/
  <simulation_name>/
    Bx.fits
    By.fits
    Bz.fits
    Vx.fits
    Vy.fits
    Vz.fits
    density.fits
    temperature.fits
    <los>/
      Synchrotron/
        ne.fits
        DM.fits
        EM.fits
        WithFaraday/
          Pmax.fits
          Qnu.fits
          Unu.fits
          FDF.fits
          realFDF.fits
          imagFDF.fits
```

The exact requirements differ by job; see [docs/jobs.md](docs/jobs.md) for the per-job breakdown.

## Configuration At A Glance

The main keys in `config/default.toml` are:

- `paths.simulations_root`
- `paths.desktop_output_root`
- `paths.transitions_csv_root`
- `simulation.name`
- `simulation.los`

Task-specific switches live under `tasks.*`.

`paths.desktop_output_root` is the preferred output key. The legacy `paths.outputs_root` key is still accepted as a fallback.

## REPL Usage

Start a project REPL:

```bash
julia --startup-file=no --project=.
```

Then:

```julia
julia> include("src/code/lib/DepolLib.jl")
julia> include("src/code/jobs/run_instrumental_effect_job.jl")
julia> using .DepolLib

julia> cfg = build_cfg("instrumental";
           config_path="config/default.toml",
           overrides=Dict(
             "paths.simulations_root" => "/path/to/simu_RAMSES",
             "paths.desktop_output_root" => "./outputs",
             "simulation.name" => "d1cf05bx10rms18000nograv1024",
             "simulation.los" => "y",
           ));

julia> result = run_instrumental_effect_job(cfg)
```

Built-in help is available from `DepolLib`:

```julia
julia> depol_help()
julia> depol_help("jobs")
julia> depol_help("instrumental")
```

## Local Repository Hygiene

This repository contains both source files and locally generated scientific outputs. To keep the tracked tree readable:

- generated outputs belong under `outputs/`
- scratch files belong under `tmp/`
- LaTeX build artefacts for the paper stay untracked

That way the repository root remains readable for external users while preserving local experimentation.
