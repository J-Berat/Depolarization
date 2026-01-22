# Depolarization

Analysis scripts for studying Faraday depolarization in synthetic synchrotron maps. The Julia notebooks in `code/` focus on polarization gradients, RM synthesis products, polarization degree, and related diagnostics derived from FITS cubes produced by RAMSES or similar MHD simulations.

## Repository layout

- `code/`: Julia scripts for the analysis and figure generation.
- `File_IO/`: Supporting data or auxiliary I/O (if present in your local copy).

## Requirements

- Julia 1.8+ (recommended).
- Packages referenced by the scripts, commonly:
  - `FITSIO`
  - `CairoMakie`
  - `FFTW`
  - `LaTeXStrings`
  - `Statistics`
  - `Printf`

Install dependencies in a Julia session:

```julia
import Pkg
Pkg.add(["FITSIO", "CairoMakie", "FFTW", "LaTeXStrings", "Statistics"])
```

## Usage

Each script is a standalone entry point. Update the hard-coded paths at the top of a script to point to your local FITS data, then run it with Julia:

```bash
julia code/GradP.jl
julia code/PolarizationDegree.jl
julia code/InstrumentalEffect.jl
```

Common inputs expected by the scripts include:

- Synchrotron Stokes cubes (`I.fits`, `Q.fits`, `U.fits`).
- Faraday synthesis products (e.g., `Qnu.fits`, `Unu.fits`).
- Magnetic field components (`Bx.fits`, `By.fits`, `Bz.fits`).
- Optional density cubes for column density maps (`density.fits`).

Outputs are written alongside the input data directories and include FITS products (e.g., `Pmax.fits`) and publication-ready figures (PDF).

## Notes

- Several scripts assume a 50 pc box sampled on 256×256 pixels. Adjust constants like `Lbox_pc`, `NPIX`, and frequency grids if your simulation differs.
- RM-synthesis parameters (frequency range and channel width) are set inside each script and should match your data.
