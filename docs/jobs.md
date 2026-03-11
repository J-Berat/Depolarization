# Job Reference

## Recommended Order

1. `run_ne_dm_em_job.jl`
2. `run_mach_suite.jl`
3. `run_reversal_transition_job.jl`
4. `run_canal_metrics_job.jl`
5. `run_segmentation_pipeline_job.jl`
6. `run_instrumental_effect_job.jl`

## Shared Conventions

- CLI: `julia --startup-file=no --project=. <job> --config ... --set ...`
- REPL task names:
  - `ne_dm_em`, `mach`, `reversal`, `canal_metrics`, `segmentation`, `instrumental`
- Standard output layout:
  - `<output_root>/<task>/<simulation>/LOS<los>/...`

## `run_ne_dm_em_job.jl`

Purpose:

- compute electron density cube (`ne`) with Wolfire approximation
- integrate LOS maps (`DM`, `EM`, electron column density)

Inputs:

- `<sim>/<density.fits>`
- `<sim>/<temperature.fits>`

Writes in simulation data tree:

- `<sim>/<los>/Synchrotron/ne.fits`
- `<sim>/<los>/Synchrotron/DM.fits`
- `<sim>/<los>/Synchrotron/EM.fits`
- `<sim>/<los>/Synchrotron/NeColumn.fits`

Writes in output tree (if `save_plot=true`):

- `ne_dm_em/.../dm_em_maps.pdf`
- `ne_dm_em/.../dm_heatmap.pdf`
- `ne_dm_em/.../em_heatmap.pdf`
- `ne_dm_em/.../ne_column_density_heatmap.pdf`

Main config:

- `tasks.ne_dm_em.*`

## `run_mach_suite.jl`

Purpose:

- compute mean sonic/Alfven Mach metrics over a fixed default simulation list
- produce scatter plots and summary CSV

Inputs:

- root cubes: `density`, `temperature`, `Bx`, `By`, `Bz`, `Vx`, `Vy`, `Vz`
- `WithFaraday/Pmax.fits` for both LOS from `tasks.mach_suite.los_y` and `tasks.mach_suite.los_x`

Outputs:

- `mach/multi/LOS<los_y>/p10_scatter.pdf`
- `mach/multi/LOS<los_y>/fraction_scatter.pdf`
- `mach/multi/LOS<los_y>/summary.csv`

Main config:

- `tasks.mach_suite.*`

## `run_reversal_transition_job.jl`

Purpose:

- compute reversal map on LOS magnetic component
- select decile pixels from `Pmax`
- derive reversal windows and FDF peak diagnostics
- generate overview plot + CSV/log summaries

Inputs:

- `<sim>/<los>/Synchrotron/WithFaraday/Pmax.fits`
- `<sim>/<los>/Synchrotron/WithFaraday/realFDF.fits`
- `<sim>/<los>/Synchrotron/WithFaraday/imagFDF.fits`
- LOS magnetic cube (`Bx`/`By`/`Bz` selected by `simulation.los`)

Outputs:

- `reversal/<sim>/LOS<los>/overview.pdf`
- `reversal/<sim>/LOS<los>/windows.csv`
- `reversal/<sim>/LOS<los>/summary.log`

Main config:

- `tasks.reversal_transition_job.*`
- `tasks.reversals_map.*`
- `tasks.b_transition.*`

## `run_canal_metrics_job.jl`

Purpose:

- compute LOS coherence maps (`Cb`, `Cphi`, `<|Bpar|>`, `Nrev`)
- compare canal vs non-canal populations (defined by low `Pmax` decile)

Inputs:

- `<sim>/<los>/Synchrotron/WithFaraday/Pmax.fits`
- `<sim>/<los>/Synchrotron/ne.fits`
- `<sim>/Bx.fits`, `<sim>/By.fits`, `<sim>/Bz.fits`

Outputs:

- `canal_metrics/<sim>/LOS<los>/maps_cb_cphi.pdf`
- `canal_metrics/<sim>/LOS<los>/pdfs_2x2.pdf`
- `canal_metrics/<sim>/LOS<los>/cb_ab_matrix.pdf`

Main config:

- `tasks.canal_metrics.*`

## `run_segmentation_pipeline_job.jl`

Purpose:

- split full cubes into regular chunks
- cut cubes around transition intervals from CSV
- produce summary intervals CSV and optional `Pmax` comparison panels

Inputs:

- transition CSV: `paths.transitions_csv_root + tasks.cut_transition.transitions_csv`
- root cubes: `Bx`, `By`, `Bz`, `Vx`, `Vy`, `Vz`, `density`, `temperature`

Outputs:

- chunk tree: `segmentation_pipeline_job_chunks/...`
- cut tree: `segmentation_pipeline_job_cuts/...`
- intervals CSV: `segmentation/.../intervals.csv`
- optional panel plot: `segmentation/.../pmax_panels.pdf`

Main config:

- `tasks.segmentation_pipeline_job.*`
- `tasks.cut_transition.*`
- `tasks.make_chunk.*`

## `run_instrumental_effect_job.jl`

Purpose:

- run spatial filtering on Q/U cubes
- compute RM-synthesis derived maps
- generate spectral diagnostics and channel alignment diagnostics

Inputs:

- `<sim>/<los>/Synchrotron/WithFaraday/Qnu.fits`
- `<sim>/<los>/Synchrotron/WithFaraday/Unu.fits`
- `<sim>/<los>/Synchrotron/WithFaraday/Pmax.fits`
- `<sim>/<los>/Synchrotron/WithFaraday/realFDF.fits`
- `<sim>/<los>/Synchrotron/WithFaraday/imagFDF.fits`
- `<sim>/Bx.fits`, `<sim>/By.fits`, `<sim>/Bz.fits`, `<sim>/density.fits`

Outputs:

- per-filter folders with `Qnu_filtered.fits`, `Unu_filtered.fits`, `RMSynthesis/Pphi_max.fits`
- `figures/*.pdf`
- `integrity_report.txt`
- `channel_alignment_summary.csv`
- `figures/channel_alignment_delta_theta_hist.pdf`

Main config:

- `tasks.instrumental_effect.*`

Detailed guide:

- `instrumental_effect.md`
