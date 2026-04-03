using MOOSE

function main(args::Vector{String})
    isempty(args) && error("Usage: moose_from_config_runner.jl <config.json>")
    MOOSE.MOOSE_from_config(args[1]; quiet=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
