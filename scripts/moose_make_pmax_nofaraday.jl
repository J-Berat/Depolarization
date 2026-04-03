using FITSIO
using MOOSE

function _arg_value(args::Vector{String}, key::String; default=nothing)
    idx = findfirst(==(key), args)
    idx === nothing && return default
    idx == length(args) && error("Missing value after $key")
    return args[idx + 1]
end

function _read_fits(path::AbstractString)
    FITS(path, "r") do f
        return read(f[1])
    end
end

function make_pmax_nofaraday(args::Vector{String})
    q_path = something(_arg_value(args, "--q"), error("Missing --q"))
    u_path = something(_arg_value(args, "--u"), error("Missing --u"))
    outdir = something(_arg_value(args, "--outdir"), error("Missing --outdir"))
    νstart = parse(Float64, something(_arg_value(args, "--nustart"), error("Missing --nustart")))
    νend = parse(Float64, something(_arg_value(args, "--nuend"), error("Missing --nuend")))
    dν = parse(Float64, something(_arg_value(args, "--dnu"), error("Missing --dnu")))
    φmin = parse(Float64, something(_arg_value(args, "--phimin"), error("Missing --phimin")))
    φmax = parse(Float64, something(_arg_value(args, "--phimax"), error("Missing --phimax")))
    dφ = parse(Float64, something(_arg_value(args, "--dphi"), error("Missing --dphi")))

    Qnu = _read_fits(q_path)
    Unu = _read_fits(u_path)
    νArray = collect(νstart:dν:νend)
    ΦArray = collect(φmin:dφ:φmax)

    FDF, realFDF, imagFDF = MOOSE.RMSynthesis(Qnu, Unu, νArray .* 1e6, ΦArray; log_progress=true)
    Pmax = MOOSE.maxCube(FDF)

    MOOSE.WriteData3D(outdir, FDF, "FDF", ΦArray; ensure_path=true)
    MOOSE.WriteData3D(outdir, realFDF, "realFDF", ΦArray; ensure_path=false)
    MOOSE.WriteData3D(outdir, imagFDF, "imagFDF", ΦArray; ensure_path=false)
    MOOSE.WriteData2D(outdir, Pmax, "Pmax"; ensure_path=false)
end

if abspath(PROGRAM_FILE) == @__FILE__
    make_pmax_nofaraday(ARGS)
end
