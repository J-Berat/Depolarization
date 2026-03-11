# Compatibility entrypoint for the instrumental pipeline internals.
include(joinpath(@__DIR__, "pipeline", "helpers.jl"))
include(joinpath(@__DIR__, "pipeline", "integrity.jl"))
include(joinpath(@__DIR__, "pipeline", "filter_pass.jl"))
include(joinpath(@__DIR__, "pipeline", "plots.jl"))
include(joinpath(@__DIR__, "pipeline", "channel_alignment.jl"))
include(joinpath(@__DIR__, "pipeline", "run_pipeline.jl"))
