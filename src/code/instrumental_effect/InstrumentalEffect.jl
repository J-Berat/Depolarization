module InstrumentalEffect

using FFTW
using FITSIO
using CairoMakie
using LaTeXStrings
using Statistics
using Printf
using Interpolations

include("../fits_io.jl")
include("config.jl")
include("io_and_axes.jl")
include("filters.jl")
include("spectra.jl")
include("plotting.jl")
include("pipeline.jl")

export InstrumentalConfig, RunFlags, load_config_defaults, run_filter_pass, run_pipeline

end
