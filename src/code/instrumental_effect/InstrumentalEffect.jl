module InstrumentalEffect

using FFTW
using FITSIO
using CairoMakie
using LaTeXStrings
using Statistics
using Printf
using Interpolations

include(joinpath(@__DIR__, "..", "fits_io.jl"))
include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "io_and_axes.jl"))
include(joinpath(@__DIR__, "filters.jl"))
include(joinpath(@__DIR__, "spectra.jl"))
include(joinpath(@__DIR__, "plotting.jl"))
include(joinpath(@__DIR__, "contracts.jl"))
include(joinpath(@__DIR__, "pipeline.jl"))

export InstrumentalConfig, RunFlags, load_config_defaults, run_filter_pass, run_pipeline

end
