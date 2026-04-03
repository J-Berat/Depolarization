using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "code", "lib", "DepolLib.jl"))
include(joinpath(REPO_ROOT, "src", "code", "instrumental_effect", "InstrumentalEffect.jl"))

using .DepolLib
using .InstrumentalEffect

include(joinpath(@__DIR__, "test_depollib.jl"))
include(joinpath(@__DIR__, "test_canal_metrics_job.jl"))
include(joinpath(@__DIR__, "test_instrumental_core.jl"))
include(joinpath(@__DIR__, "test_instrumental_analysis_helpers.jl"))
include(joinpath(@__DIR__, "test_filter_pass.jl"))
include(joinpath(@__DIR__, "test_instrumental_job.jl"))
