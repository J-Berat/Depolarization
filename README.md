# Depolarization

Julia workflows for depolarization and synchrotron/Faraday analysis on RAMSES simulation products.

## Quick Start

```bash
# 1) bootstrap local Julia environment (installs + precompiles)
julia --startup-file=no --project=. scripts/precompile.jl

# 2) run one job on a target simulation (with runtime overrides)
julia --startup-file=no --project=. src/code/jobs/run_instrumental_effect_job.jl \
  --config config/default.toml \
  --set paths.simulations_root=/path/to/simu_RAMSES \
  --set simulation.name=d1cf05bx10rms18000nograv1024 \
  --set simulation.los=y

# 3) inspect outputs
ls -R outputs/instrumental_effect/d1cf05bx10rms18000nograv1024/LOSy
```

## Project Layout

- `src/code/lib/DepolLib.jl`: shared config, FITS, LOS, and output helpers
- `src/code/jobs/`: consolidated multi-step pipelines
- `src/code/instrumental_effect/`: instrumental filtering module
- `src/code/fits_io.jl`, `src/code/los_utils.jl`: low-level shared support
- `scripts/precompile.jl`: local dependency install + precompile helper
- `config/default.toml`: runtime configuration
- `docs/functions.md`: function index for the codebase

## Requirements

- Julia 1.10+
- Julia packages: `CairoMakie`, `FITSIO`, `LaTeXStrings`, `StatsBase`, `Makie`, `FFTW`, `Interpolations`

Project-local install (without using the helper script):

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.add(["CairoMakie","FITSIO","LaTeXStrings","StatsBase","Makie","FFTW","Interpolations"]); Pkg.precompile()'
```

## Precompile (Recommended)

Run this before first use (and after major updates) to reduce startup latency:

```bash
julia --startup-file=no --project=. scripts/precompile.jl
```

This script activates the repository root and creates/updates local `Project.toml` and `Manifest.toml` when needed.

Dependencies-only variant:

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

## Configuration

Edit `config/default.toml` before running. Key groups:

- `paths`
- `simulation`
- `tasks`

Most-used keys:

- `paths.simulations_root`
- `paths.outputs_root`
- `paths.transitions_csv_root`
- `simulation.name`
- `simulation.los` (`x`, `y`, `z`)

All scripts that use the shared runtime loader accept:

```bash
julia --startup-file=no --project=. <script>.jl --config config/default.toml --set key=value --set other.key=value
```

- `--config`: TOML file path (default: `config/default.toml`)
- `--set`: repeatable dotted-key override

## Main Job Entrypoints

```bash
julia --startup-file=no --project=. src/code/jobs/run_mach_suite.jl --config config/default.toml
julia --startup-file=no --project=. src/code/jobs/run_reversal_transition_job.jl --config config/default.toml
julia --startup-file=no --project=. src/code/jobs/run_segmentation_pipeline_job.jl --config config/default.toml
julia --startup-file=no --project=. src/code/jobs/run_instrumental_effect_job.jl --config config/default.toml
```

## Run All Jobs (Batch)

```bash
for job in \
  src/code/jobs/run_mach_suite.jl \
  src/code/jobs/run_reversal_transition_job.jl \
  src/code/jobs/run_segmentation_pipeline_job.jl \
  src/code/jobs/run_instrumental_effect_job.jl
do
  julia --startup-file=no --project=. "$job" --config config/default.toml
done
```

Task keys used by these jobs live under:

- `tasks.mach_suite`
- `tasks.reversal_transition_job`
- `tasks.reversals_map`
- `tasks.b_transition`
- `tasks.segmentation_pipeline_job`
- `tasks.cut_transition`
- `tasks.make_chunk`
- `tasks.instrumental_effect`

## Job Outputs

- `run_mach_suite.jl`
  - `<outputs_root>/mach_suite/multi/LOS<los_y>/mach_suite_multi_LOS<los_y>_p10_scatter.pdf`
  - `<outputs_root>/mach_suite/multi/LOS<los_y>/mach_suite_multi_LOS<los_y>_fraction_scatter.pdf`
  - `<outputs_root>/mach_suite/multi/LOS<los_y>/mach_suite_multi_LOS<los_y>_summary.csv`
- `run_reversal_transition_job.jl`
  - `<outputs_root>/reversal_transition_job/<simulation>/LOS<los>/..._overview.pdf`
  - `<outputs_root>/reversal_transition_job/<simulation>/LOS<los>/..._windows.csv`
  - `<outputs_root>/reversal_transition_job/<simulation>/LOS<los>/..._summary.log`
- `run_segmentation_pipeline_job.jl`
  - chunked FITS under `<outputs_root>/segmentation_pipeline_job_chunks/<simulation>/LOS<los>/`
  - cut FITS under `<outputs_root>/segmentation_pipeline_job_cuts/<simulation>/LOS<los>/`
  - `<outputs_root>/segmentation_pipeline_job/<simulation>/LOS<los>/..._intervals.csv`
  - `<outputs_root>/segmentation_pipeline_job/<simulation>/LOS<los>/..._pmax_panels.pdf`
- `run_instrumental_effect_job.jl`
  - per-filter folders under `<outputs_root>/instrumental_effect/<simulation>/LOS<los>/HardBandPass_*`
  - each folder writes `Qnu_filtered.fits`, `Unu_filtered.fits`, and `RMSynthesis/Pphi_max.fits`

## Expected Input Data Layout

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
        WithFaraday/
          Pmax.fits
          Qnu.fits
          Unu.fits
          FDF.fits
          realFDF.fits
          imagFDF.fits
```

Transition interval CSV files are read from `paths.transitions_csv_root`.

## Output Convention

Top-level output directories are standardized to:

```text
<outputs_root>/<task>/<simulation>/LOS<los>/
```

Filename pattern:

```text
<task>_<simulation>_LOS<los>_<artifact>.<ext>
```

Examples:

```text
outputs/mach_suite/multi/LOSy/mach_suite_multi_LOSy_summary.csv
outputs/reversal_transition_job/d1cf05bx10rms18000nograv1024/LOSy/reversal_transition_job_d1cf05bx10rms18000nograv1024_LOSy_windows.csv
```

Note: some pipelines (notably `instrumental_effect` and segmentation subproducts) also create nested files/folders inside the standardized task directory.

## Developer Notes

- Job entrypoints under `src/code/jobs/` should keep the same runtime contract:
  - load config via shared helpers (`--config`, repeatable `--set`)
  - validate inputs early
  - return a structured `Dict` summary
- Shared runtime/path/output helpers belong in `src/code/lib/DepolLib.jl`.
- `src/code/common.jl` is compatibility-only; new code should prefer direct `DepolLib` usage.
- Inside `src/code/instrumental_effect/`, keep boundaries explicit:
  - compute + IO contracts in `pipeline.jl` and `contracts.jl`
  - plotting-only behavior in `plotting.jl`
  - low-level Fourier/filter helpers in `filters.jl`, `spectra.jl`, `io_and_axes.jl`

## Troubleshooting

- Missing FITS files: validate `paths.simulations_root`, `simulation.name`, `simulation.los`, and required files in the layout above.
- LOS errors: only `x`, `y`, `z` are accepted.
- Transition pipeline errors: verify `paths.transitions_csv_root` and CSV columns `kmin,kmax`.
- Shape mismatch errors: confirm all related cubes/maps have compatible dimensions.
