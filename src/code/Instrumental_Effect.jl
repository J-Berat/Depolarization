include("instrumental_effect/InstrumentalEffect.jl")

using .InstrumentalEffect

cfg = load_config_defaults()
flags = RunFlags()

run_pipeline(cfg; flags=flags)
