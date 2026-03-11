# Development Guide

## Repository Structure

- `src/code/lib/DepolLib.jl`: shared runtime/config/path helpers and REPL helper text
- `src/code/common.jl`: compatibility shim re-exporting `DepolLib`
- `src/code/jobs/`: executable job entrypoints
- `src/code/instrumental_effect/`: instrumental module internals
- `scripts/precompile.jl`: dependency setup and warmup
- `scripts/run_with_precompile.sh`: convenience wrapper for precompile+job run
- `docs/functions.md`: auto-generated function list
- `test/`: unit/integration-style tests

## Instrumental Pipeline Layout

`src/code/instrumental_effect/pipeline/` is split by responsibility:

- `helpers.jl`: utility helpers and figure save helpers
- `integrity.jl`: input/invariant checks and integrity report writing
- `filter_pass.jl`: filtering and per-filter FITS generation
- `plots.jl`: spectrum/map plotting routines
- `channel_alignment.jl`: canal/B_perp alignment diagnostics and histogram output
- `run_pipeline.jl`: orchestration based on `RunFlags`

`src/code/instrumental_effect/pipeline.jl` is a compatibility include entrypoint.

## Running Tests

Full test suite:

```bash
julia --startup-file=no --project=. test/runtests.jl
```

CI uses the same command on Julia 1.12.

## Regenerating Function Docs

```bash
./scripts/generate_functions_doc.sh
```

This regenerates `docs/functions.md` from all `*.jl` files in `src/code/`.

## Adding a New Job

1. Create `src/code/jobs/run_<name>_job.jl` with:
   - `run_<name>_job(cfg)::Dict{String,Any}`
   - guarded entrypoint:
     - `if abspath(PROGRAM_FILE) == @__FILE__`
     - `run_job_entrypoint("<task>", run_<name>_job)`
2. Reuse shared helpers from `DepolLib` for:
   - config reads (`cfg_get`, `cfg_require`)
   - LOS and paths (`require_los`, `simulation_field_path`, `standard_output_path`)
   - validation (`require_existing_files`, `require_ndims`, `require_same_size`)
3. Add tests in `test/` and include them in `test/runtests.jl`.
4. Update docs:
   - `docs/jobs.md`
   - `README.md` (job lists/order)
   - `docs/functions.md` (regenerate)

## REPL Helper Text

The in-REPL guidance is generated in `DepolLib.help` (`depol_help` and `repl_help` aliases).
If job names/order/config sections change, keep that helper synchronized with docs.
