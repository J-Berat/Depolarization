# Depolarization

This repository contains Julia scripts for depolarization and synchrotron/Faraday analysis.

## Instrumental effect workflow (modular)

`src/code/Instrumental_Effect.jl` is now a thin runner that calls a reusable module:

- `src/code/instrumental_effect/InstrumentalEffect.jl` (module entrypoint)
- `src/code/instrumental_effect/config.jl` (typed configuration + run flags)
- `src/code/instrumental_effect/io_and_axes.jl` (grid/axes/channel helpers)
- `src/code/instrumental_effect/filters.jl` (instrument Fourier filtering)
- `src/code/instrumental_effect/spectra.jl` (PSD and peak helpers)
- `src/code/instrumental_effect/plotting.jl` (themes/ticks/plot helpers)
- `src/code/instrumental_effect/pipeline.jl` (orchestration)

## Public API

The module exports:

- `InstrumentalConfig`
- `RunFlags`
- `load_config_defaults()`
- `run_filter_pass(cfg)`
- `run_pipeline(cfg; flags=RunFlags())`

## Run

Default runner (same script entrypoint as before):

```bash
julia src/code/Instrumental_Effect.jl
```

## Configure parameters

Example with custom parameters and selective sections:

```julia
include("src/code/instrumental_effect/InstrumentalEffect.jl")
using .InstrumentalEffect

cfg = load_config_defaults()
cfg = InstrumentalConfig(
    Q_in = cfg.Q_in,
    U_in = cfg.U_in,
    Pmax_nofilter_path = cfg.Pmax_nofilter_path,
    Q_in_phi = cfg.Q_in_phi,
    U_in_phi = cfg.U_in_phi,
    base_out = cfg.base_out,
    Bx_in = cfg.Bx_in,
    By_in = cfg.By_in,
    Bz_in = cfg.Bz_in,
    dens_in = cfg.dens_in,
    n = cfg.n,
    m = cfg.m,
    Lbox_pc = cfg.Lbox_pc,
    Lcut_small = 1.5,
    Llarge_list = [80.0, 40.0, 20.0],
    νmin_MHz = cfg.νmin_MHz,
    νmax_MHz = cfg.νmax_MHz,
    Δν_MHz = cfg.Δν_MHz,
    PhiArray = cfg.PhiArray,
    ichan = 40,
    iphi = 25,
)

flags = RunFlags(
    run_pmax_maps = true,
    run_psd = true,
    run_q_u_p_q2 = false,
    run_phi_q_u_p = false,
    run_lic = false,
)

run_pipeline(cfg; flags=flags)
```

## Notes

- Default values in `load_config_defaults()` preserve the previous behavior and paths.
- `run_filter_pass(cfg)` can be called independently if only filtering + RM-synthesis outputs are needed.
