# Instrumental Effect Pipeline

This document defines the runtime contract for:

- `src/code/jobs/run_instrumental_effect_job.jl`
- `src/code/instrumental_effect/InstrumentalEffect.jl`

## Run Command

```bash
julia --startup-file=no --project=. src/code/jobs/run_instrumental_effect_job.jl \
  --config config/default.toml \
  --set paths.simulations_root=/path/to/simu_RAMSES \
  --set paths.desktop_output_root=./outputs \
  --set simulation.name=d1cf05bx10rms18000nograv1024 \
  --set simulation.los=y
```

Only `--config` and repeatable `--set key=value` are accepted.

## Required Inputs

For `<simulations_root>/<simulation>/<los>/Synchrotron/WithFaraday/`:

- `Qnu.fits`
- `Unu.fits`
- `Pmax.fits`
- `realFDF.fits`
- `imagFDF.fits`

For `<simulations_root>/<simulation>/`:

- `Bx.fits`
- `By.fits`
- `Bz.fits`
- `density.fits`

Missing files fail early with a dedicated error message.

## Config Keys Used

Main routing keys:

- `paths.simulations_root`
- `paths.desktop_output_root` (preferred, legacy fallback: `paths.outputs_root`)
- `simulation.name`
- `simulation.los`

Task keys under `tasks.instrumental_effect`:

- `llarge_list` (default list in code)
- `run_pmax_maps` (bool)
- `run_psd` (bool)
- `run_q_u_p_q2` (bool)
- `run_phi_q_u_p` (bool)
- `run_lic` (bool, currently kept disabled by default behavior)
- `run_channel_b_alignment` (bool)
- `channel_alignment_pdf_plain` (bool, `true` saves a clean PDF without slider UI)

## Output Layout

Base output directory:

```text
<desktop_output_root>/instrumental/<simulation>/LOS<los>/
```

Per filter scale `<Llarge>`:

```text
HardBandPass_remove_L0_to_1pc_and_<Llarge>to50pc/
  Qnu_filtered.fits
  Unu_filtered.fits
  RMSynthesis/
    Pphi_max.fits
```

Global outputs in base directory:

- `figures/pmax_maps.pdf` (if `run_pmax_maps=true`)
- `figures/psd_panels.pdf` and `figures/pmax_kx.pdf` (if `run_psd=true`)
- `figures/q_nu50.pdf`, `u_nu50.pdf`, `p_nu50.pdf`, `q2_nu50.pdf` (if `run_q_u_p_q2=true`)
- `figures/q_phi.pdf`, `u_phi.pdf`, `p_phi.pdf` (if `run_phi_q_u_p=true`)
- `integrity_report.txt` (always written before heavy compute starts)
- `channel_alignment_summary.csv` and `figures/channel_alignment_delta_theta_hist.pdf` (if `run_channel_b_alignment=true`)

## Integrity Report

The pipeline writes `integrity_report.txt` with:

- global status: `pass` or `fail`
- critical failures count
- warning count
- total checks count
- per-check lines with `level`, `passed`, and message

On critical check failure, the run aborts after the report is written.

Examples of checks:

- axis sanity (`nu`, `phi` non-empty/finiteness/ordering)
- shape/layout checks for `Qnu`, `Unu`, `Pmax`
- finite-value checks (no `NaN`/`Inf`)
- physical plausibility checks for frequency and phi ranges

## Structured Logging Fields

Recent logs are emitted with explicit key/value fields so they are easier to parse:

- grid geometry: `nx`, `ny`, `dx_pc`, `dy_pc`, `f_nyquist_rad_pc_inv`
- filter frequencies: `Lcut_small_pc`, `Llarge_pc`, `flo_rad_pc_inv`, `fhi_rad_pc_inv`
- phi channel selection: `iphi`, `phi_rad_m2`
- channel/B alignment completion: `summary_path`, `figure_path`, `los`

## Troubleshooting

- If outputs appear on Desktop unexpectedly, set `paths.desktop_output_root` explicitly.
- If you still use `paths.outputs_root`, migrate to `paths.desktop_output_root`.
- If run aborts quickly, check `integrity_report.txt` first.
- If `simulation.los` is invalid, only `x`, `y`, `z` are accepted.
