# Depolarization

Julia workflows for depolarization and synchrotron/Faraday analysis on RAMSES simulation data.

## Quick Start (CLI)

```bash
# 1) Install + precompile local environment
julia --startup-file=no --project=. scripts/precompile.jl

# 2) Run one job with runtime overrides
julia --startup-file=no --project=. src/code/jobs/run_instrumental_effect_job.jl \
  --config config/default.toml \
  --set paths.simulations_root=/path/to/simu_RAMSES \
  --set paths.desktop_output_root=./outputs \
  --set simulation.name=d1cf05bx10rms18000nograv1024 \
  --set simulation.los=y

# 3) Check generated outputs
ls -R ./outputs/instrumental/d1cf05bx10rms18000nograv1024/LOSy
```

All jobs support:

```bash
julia --startup-file=no --project=. <job>.jl --config config/default.toml --set key=value --set other.key=value
```

Or run with forced precompilation every time:

```bash
./scripts/run_with_precompile.sh reversal_transition --config config/default.toml
```

Accepted CLI args are only `--config` and repeatable `--set key=value`.

## Documentation

- [docs/README.md](docs/README.md): documentation index and reading map.
- [docs/getting_started.md](docs/getting_started.md): setup, first run, CLI/REPL workflows.
- [docs/configuration.md](docs/configuration.md): config keys, overrides, and task sections.
- [docs/jobs.md](docs/jobs.md): complete per-job inputs/outputs and execution order.
- [docs/development.md](docs/development.md): architecture, tests, and maintenance workflow.
- [docs/instrumental_effect.md](docs/instrumental_effect.md): full pipeline contract for `run_instrumental_effect_job.jl` (inputs, flags, outputs, integrity report, logs).
- [docs/functions.md](docs/functions.md): auto-generated function index across `src/code/`.

Regenerate function docs:

```bash
./scripts/generate_functions_doc.sh
```

## Use From REPL

Start REPL in project mode:

```bash
julia --startup-file=no --project=.
```

Example:

```julia
julia> cd("/absolute/path/to/Depolarization")

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
julia> result["output_dir"]
```

CLI-style config in REPL is also supported:

```julia
julia> cfg = load_runtime_config("reversal";
           args=[
             "--config", "config/default.toml",
             "--set", "paths.simulations_root=/path/to/simu_RAMSES",
             "--set", "simulation.name=d1cf05bx10rms18000nograv1024",
             "--set", "simulation.los=y",
           ]);
```

Built-in REPL help (from `DepolLib`):

```julia
julia> include("src/code/lib/DepolLib.jl")
julia> using .DepolLib

julia> depol_help()               # full guide (setup + jobs)
julia> depol_help("jobs")         # all jobs
julia> depol_help("instrumental") # one specific job

julia> repl_help("reversal")       # alias of depol_help(...)
julia> DepolLib.help("reversal")   # qualified fallback if needed
```

## Jobs

- `src/code/jobs/run_ne_dm_em_job.jl`
- `src/code/jobs/run_mach_suite.jl`
- `src/code/jobs/run_reversal_transition_job.jl`
- `src/code/jobs/run_canal_metrics_job.jl`
- `src/code/jobs/run_segmentation_pipeline_job.jl`
- `src/code/jobs/run_instrumental_effect_job.jl`

Recommended order (fresh simulation):

1. `run_ne_dm_em_job.jl` (builds `ne.fits`, `DM.fits`, `EM.fits`)
2. `run_mach_suite.jl`
3. `run_reversal_transition_job.jl`
4. `run_canal_metrics_job.jl`
5. `run_segmentation_pipeline_job.jl`
6. `run_instrumental_effect_job.jl`

CLI cheatsheet (one command per job):

```bash
julia --startup-file=no --project=. src/code/jobs/run_ne_dm_em_job.jl --config config/default.toml
julia --startup-file=no --project=. src/code/jobs/run_mach_suite.jl --config config/default.toml
julia --startup-file=no --project=. src/code/jobs/run_reversal_transition_job.jl --config config/default.toml
julia --startup-file=no --project=. src/code/jobs/run_canal_metrics_job.jl --config config/default.toml
julia --startup-file=no --project=. src/code/jobs/run_segmentation_pipeline_job.jl --config config/default.toml
julia --startup-file=no --project=. src/code/jobs/run_instrumental_effect_job.jl --config config/default.toml
```

REPL task names for `build_cfg(task, ...)` / `load_runtime_config(task, ...)`:

- `ne_dm_em`
- `mach`
- `reversal`
- `canal_metrics`
- `segmentation`
- `instrumental`

Run all jobs:

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

## Minimal Config

Main keys in `config/default.toml`:

- `paths.simulations_root`
- `paths.desktop_output_root` (optional; default is `~/Desktop/depolarization_outputs`)
- `paths.transitions_csv_root`
- `simulation.name`
- `simulation.los` (`x`, `y`, `z`)

Task toggles/params live under `tasks.*` (for example `tasks.canal_metrics`, `tasks.instrumental_effect`, `tasks.ne_dm_em`, `tasks.segmentation_pipeline_job`).

`paths.desktop_output_root` is the preferred key.
Legacy `paths.outputs_root` is still accepted as a fallback.

## Expected Inputs

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
        ne.fits                # generated by run_ne_dm_em_job.jl (or provided)
        DM.fits                # generated by run_ne_dm_em_job.jl
        EM.fits                # generated by run_ne_dm_em_job.jl
        WithFaraday/
          Pmax.fits
          Qnu.fits
          Unu.fits
          FDF.fits
          realFDF.fits
          imagFDF.fits
```

Required by job:

- `run_ne_dm_em_job.jl`: `density.fits`, `temperature.fits`
- `run_mach_suite.jl`: root cubes + LOS `WithFaraday/Pmax.fits` for both `tasks.mach_suite.los_y` and `tasks.mach_suite.los_x`
- `run_reversal_transition_job.jl`: `WithFaraday/Pmax.fits`, `WithFaraday/FDF.fits`, `Synchrotron/ne.fits`, and LOS magnetic cube (`Bx`/`By`/`Bz` by LOS)
- `run_canal_metrics_job.jl`: `WithFaraday/Pmax.fits`, `Synchrotron/ne.fits`, and root `Bx.fits/By.fits/Bz.fits`
- `run_segmentation_pipeline_job.jl`: root cubes + transitions CSV under `paths.transitions_csv_root` with columns `kmin,kmax`
- `run_instrumental_effect_job.jl`: `WithFaraday/Qnu.fits`, `Unu.fits`, `Pmax.fits`, `realFDF.fits`, `imagFDF.fits` + root `Bx/By/Bz/density`

## Outputs

Top-level output folder pattern:

```text
<desktop_output_root>/<task>/<simulation>/LOS<los>/
```

For the instrumental pipeline, each filter also writes:

- `HardBandPass_remove_L0_to_1pc_and_<L>to50pc/Qnu_filtered.fits`
- `HardBandPass_remove_L0_to_1pc_and_<L>to50pc/Unu_filtered.fits`
- `HardBandPass_remove_L0_to_1pc_and_<L>to50pc/RMSynthesis/Pphi_max.fits`
- `figures/*.pdf`
- `channel_alignment_summary.csv` (if `tasks.instrumental_effect.run_channel_b_alignment=true`)
- `figures/channel_alignment_delta_theta_hist.pdf` (if `tasks.instrumental_effect.run_channel_b_alignment=true`)
- `tasks.instrumental_effect.channel_alignment_pdf_plain=true` saves this PDF without slider UI (histograms only).
- `integrity_report.txt`

For `canal_metrics`, outputs include:

- `maps_cb_cphi.pdf`
- `pdfs_2x2.pdf`
- `cb_ab_matrix.pdf`

Task names produced by current jobs:

- `ne_dm_em`
- `mach`
- `reversal`
- `canal_metrics`
- `segmentation`
- `instrumental`

## Troubleshooting

- Missing FITS files: verify `paths.simulations_root`, `simulation.name`, `simulation.los`.
- Invalid LOS: only `x`, `y`, `z`.
- Transition CSV errors: verify `paths.transitions_csv_root` and CSV columns `kmin,kmax`.
- CLI parsing errors: use only `--config` and repeatable `--set key=value`.
