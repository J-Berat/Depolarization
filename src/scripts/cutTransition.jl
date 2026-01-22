using FITSIO
using Printf
using Statistics

include(joinpath(@__DIR__, "../src/utils/fits_utils.jl"))

# ============================================================
# CUT CUBES BETWEEN (kmin,kmax) INTERVALS FROM A CSV
# + intervals CSV includes distances in pc
# ============================================================

# ------------------------------------------------------------
# USER CHOICES
# ------------------------------------------------------------

const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS       = "y"  # "x","y","z"
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"

# CSV des transitions (décile 1)
const TRANSITIONS_CSV =
"/Users/jb270005/Desktop/Depolarization_canals/sightline_plots/BLOS_decile1_LOSy_transitions.csv"

# OUTDIR (exactement comme tu veux)
const OUT_ROOT =
"/Users/jb270005/Desktop/Depolarization_canals/d1cf05bx10rms18000nograv1024"

# Geometry (for pc conversion)
const NPIX    = 256
const LBOX_PC = 50.0

# Fusion si kmin_next <= kmax_prev + MERGE_GAP_PIX
const MERGE_GAP_PIX = 1

const OVERWRITE = true

# Cubes à découper
const VARS = ["Bx","By","Bz","density","Temperature","Vx","Vy","Vz"]

const VAR_FILEMAP = Dict(
    "Bx" => "Bx.fits",
    "By" => "By.fits",
    "Bz" => "Bz.fits",
    "density" => "density.fits",
    "Temperature" => "Temperature.fits",
    "Vx" => "Vx.fits",
    "Vy" => "Vy.fits",
    "Vz" => "Vz.fits",
)

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------
function read_intervals_from_csv(path::String)
    @assert isfile(path) "Missing CSV: $path"
    lines = readlines(path)
    @assert !isempty(lines) "Empty CSV: $path"

    header = split(strip(lines[1]), ',')
    col = Dict(h => i for (i,h) in enumerate(header))
    @assert haskey(col,"kmin") && haskey(col,"kmax") "CSV must contain kmin,kmax"

    intervals = Tuple{Int,Int}[]
    for ln in lines[2:end]
        s = strip(ln)
        isempty(s) && continue
        parts = split(s, ',')
        push!(intervals,
              (parse(Int, parts[col["kmin"]]),
               parse(Int, parts[col["kmax"]])))
    end
    sort!(intervals, by=x->x[1])
    return intervals
end

function merge_intervals(intervals; gap_pix::Int=1)
    isempty(intervals) && return Tuple{Int,Int}[]
    out = Tuple{Int,Int}[]
    kL, kR = intervals[1]
    for (L,R) in intervals[2:end]
        if L <= kR + gap_pix
            kR = max(kR, R)
        else
            push!(out, (kL, kR))
            kL, kR = L, R
        end
    end
    push!(out, (kL, kR))
    return out
end

los_axis(los::String) = los=="x" ? 1 : los=="y" ? 2 : los=="z" ? 3 : error("LOS must be x/y/z")

function cut_cube(A::Array{T,3}, ax::Int, kmin::Int, kmax::Int) where {T}
    if ax == 1
        return Array(@view A[kmin:kmax, :, :])
    elseif ax == 2
        return Array(@view A[:, kmin:kmax, :])
    else
        return Array(@view A[:, :, kmin:kmax])
    end
end

# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

SIMU_DIR = joinpath(SIMU_ROOT, SIMU_NAME)
@assert isdir(SIMU_DIR) "Missing simulation dir: $SIMU_DIR"

# Distance grid (pc) matching your other scripts
dist = collect(range(0.0, LBOX_PC; length=NPIX))

intervals_raw = read_intervals_from_csv(TRANSITIONS_CSV)
@assert !isempty(intervals_raw) "No intervals found in CSV: $TRANSITIONS_CSV"

intervals = merge_intervals(intervals_raw; gap_pix=MERGE_GAP_PIX)

@info "Intervals" raw=length(intervals_raw) merged=length(intervals)

# Output folder
out_dir = joinpath(OUT_ROOT, "cuts_LOS_" * LOS)
mkpath(out_dir)

# ------------------------------------------------------------
# CSV of merged intervals + distances in pc
# ------------------------------------------------------------
intervals_csv = joinpath(out_dir, "intervals_kmin_kmax_pc.csv")
open(intervals_csv, "w") do io
    println(io, "segment,kmin,kmax,smin_pc,smax_pc,width_pc,npix")
    for (seg,(kmin,kmax)) in enumerate(intervals)
        @assert 1 ≤ kmin ≤ kmax ≤ NPIX "kmin/kmax out of bounds: ($kmin,$kmax) with NPIX=$NPIX"
        smin = dist[kmin]
        smax = dist[kmax]
        width = smax - smin
        npix = kmax - kmin + 1
        println(io, @sprintf("%d,%d,%d,%.6f,%.6f,%.6f,%d", seg, kmin, kmax, smin, smax, width, npix))
    end
end
@info "Saved merged intervals CSV with pc distances" path=intervals_csv

# ------------------------------------------------------------
# CUT & WRITE
# ------------------------------------------------------------
ax = los_axis(LOS)

for (seg,(kmin,kmax)) in enumerate(intervals)
    segdir = joinpath(out_dir, @sprintf("seg_%02d_k%04d_%04d", seg, kmin, kmax))
    mkpath(segdir)

    @info "Cutting segment" seg=seg kmin=kmin kmax=kmax

    for var in VARS
        fname = VAR_FILEMAP[var]
        inpath = joinpath(SIMU_DIR, fname)
        if !isfile(inpath)
            @warn "Missing file, skipped" var=var path=inpath
            continue
        end

        A = read_FITS(inpath)
        @assert ndims(A) == 3 "Expected 3D cube for $var, got size=$(size(A))"

        Acut = cut_cube(A, ax, kmin, kmax)

        outpath = joinpath(segdir, fname)
        write_FITS(outpath, Acut; overwrite=OVERWRITE)
    end
end

@info "DONE — cubes written to" out_dir
