# Configuration Reference

## Runtime Model

Jobs load configuration through shared helpers in `DepolLib`:

- `load_runtime_config(task; args=...)`: CLI parsing (`--config`, `--set`)
- `build_cfg(task; config_path=..., overrides=...)`: REPL-friendly configuration

Runtime metadata is written under `cfg["runtime"]`:

- `task`
- `config_path`
- `set_args`

## Root Sections

`config/default.toml` is organized around:

- `[paths]`
- `[simulation]`
- `[tasks.*]`

### `[paths]`

- `simulations_root`: root folder of simulation datasets
- `desktop_output_root`: preferred root for generated outputs
- `transitions_csv_root`: folder containing transition interval CSV files

Legacy key:

- `outputs_root`: still supported as fallback, but deprecated

### `[simulation]`

- `name`: simulation folder name
- `los`: line of sight (`x`, `y`, or `z`)

## Task Sections

### `[tasks.ne_dm_em]`

- `enabled`
- `prompt_constants`
- physical constants and geometry: `lbox_pc`, `pc_to_cm`, `zeta`, `Geff`, `phiPAH`, `XC`
- writing/plot style: `overwrite`, `save_plot`, sizes/ticks/colorbar options

### `[tasks.mach_suite]`

- `los_y`
- `los_x`

### `[tasks.reversal_transition_job]`

- `enabled`
- optional controls such as smoothing and thresholds (`smooth_win`, `sign_eps`, `deriv_tol`, `rng_seed`)
- optional phi controls (`phiArray`, `PhiArray`, `phi_min`, `phi_max`, `nphi`)

Related sections used by this job:

- `[tasks.reversals_map]` (`npix`, `lbox_pc`, `b_scale`, `target_nrev`)
- `[tasks.b_transition]` (`nphi`, optional phi defaults)

### `[tasks.canal_metrics]`

- `enabled`
- `decile`
- `b_unit_factor`
- `lbox_pc`
- `nrev_thresh`
- histogram/binning controls: `nbins_pdf`, `nbins_nrev`, `nbins_stats`, `nbins_matrix_x`, `nbins_matrix_y`, `mincount`

### `[tasks.segmentation_pipeline_job]`

- `enabled`

Related sections:

- `[tasks.cut_transition]`: transition CSV name and cut behavior
- `[tasks.make_chunk]`: `chunk_size`

### `[tasks.instrumental_effect]`

- plotting/analysis toggles:
  - `run_pmax_maps`
  - `run_psd`
  - `run_q_u_p_q2`
  - `run_phi_q_u_p`
  - `run_lic`
  - `run_channel_b_alignment`
- channel alignment rendering mode:
  - `channel_alignment_pdf_plain`
- filter scales:
  - `llarge_list`

## CLI Overrides

Any key can be overridden with dotted-path syntax:

```bash
--set paths.desktop_output_root=./outputs
--set simulation.los=z
--set tasks.instrumental_effect.run_psd=false
```

Scalar parsing rules:

- `true`/`false` -> `Bool`
- integer strings -> `Int`
- decimal strings -> `Float64`
- `nothing` -> `nothing`
- otherwise -> `String`

## Validation Notes

Shared helpers enforce:

- LOS in `{x,y,z}`
- required keys for active jobs
- required file existence before heavy computation

Instrumental pipeline adds integrity checks and writes `integrity_report.txt`.
