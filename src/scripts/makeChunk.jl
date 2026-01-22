using FITSIO
using Printf

include(joinpath(@__DIR__, "../src/utils/fits_utils.jl"))

# -----------------------------------------------------------
# PARAMETERS
# -----------------------------------------------------------

const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"
const OUT_ROOT  = "/Users/jb270005/Desktop/Depolarization_canals"

# Choose the chunk size here (e.g. 32, 64, 128)
const CHUNK_SIZE = 32

# Simulation list
SimuList = [
    "d1cf00bx10rms18000nograv1024",
    "d1cf02bx10rms09000nograv1024",
    "d1cf02bx10rms18000nograv1024",
    "d1cf02bx10rms36000nograv1024",
    "d1cf02bx10rms72000nograv1024",
    "d1cf05bx10rms09000nograv1024",
    "d1cf05bx10rms18000nograv1024", # simulation 7
    "d1cf05bx10rms36000nograv1024",
    "d1cf05bx10rms72000nograv1024",
    "d1cf08bx10rms09000nograv1024",
    "d1cf08bx10rms18000nograv1024",
    "d1cf08bx10rms36000nograv1024",
    "d1cf08bx10rms72000nograv1024",
    "d1cf10bx10rms18000nograv1024",
]

# Fields to read directly from FITS
const FIELDS = [
    "Bx", "By", "Bz",
    "Vx", "Vy", "Vz",
    "density", "temperature",
]

# -----------------------------------------------------------
# GENERIC SPLIT FUNCTION
# -----------------------------------------------------------

function split_along_dim(
    data::AbstractArray,
    outdir::String;
    base::String,
    dim::Int,
    chunk_size::Int,
)
    n = size(data, dim)
    @assert n % chunk_size == 0 "Dimension $dim size=$n is not a multiple of chunk_size=$chunk_size"

    n_chunks = div(n, chunk_size)
    mkpath(outdir)

    for i in 1:n_chunks
        i1 = (i - 1) * chunk_size + 1
        i2 = i * chunk_size

        ranges = ntuple(d -> d == dim ? (i1:i2) : Colon(), ndims(data))
        slice  = Array(@view data[ranges...])  # contiguous for FITSIO

        chunk_dir = joinpath(outdir, @sprintf("chunk%02d", i))
        mkpath(chunk_dir)

        outfile = joinpath(chunk_dir, base * ".fits")
        println("    Writing $outfile")

        FITS(outfile, "w") do f
            write(f, slice)  # minimal FITS header, no copy
        end
    end
end

# -----------------------------------------------------------
# MAIN — ONLY SIMULATION #7
# -----------------------------------------------------------

simu_name = SimuList[7]
println("=== Processing simulation #7: $simu_name ===")

simu_path = joinpath(SIMU_ROOT, simu_name)

# This creates a size-specific folder layer, e.g. "64pix"
size_tag = @sprintf("%dpix", CHUNK_SIZE)

# Two separate output trees under the size folder:
# - split_x : chunks along dimension 1 -> (CHUNK_SIZE, 256, 256)
# - split_y : chunks along dimension 2 -> (256, CHUNK_SIZE, 256)
# - split_z : chunks along dimension 3 -> (256, 256, CHUNK_SIZE)
out_x = joinpath(OUT_ROOT, simu_name, size_tag, "split_x")
out_y = joinpath(OUT_ROOT, simu_name, size_tag, "split_y")
out_z = joinpath(OUT_ROOT, simu_name, size_tag, "split_z")

mkpath(out_x)
mkpath(out_y)
mkpath(out_z)

for field in FIELDS
    println("  Field: $field")

    infile = joinpath(simu_path, field * ".fits")
    @assert isfile(infile) "Missing file: $infile"

    data = read_FITS(infile)

    # Sanity check (expected 256×256×256)
    @assert size(data) == (256, 256, 256) "Unexpected cube size for $field: got $(size(data))"

    # Split along first dimension (x)
    split_along_dim(data, out_x; base=field, dim=1, chunk_size=CHUNK_SIZE)

    # Split along second dimension (y)
    split_along_dim(data, out_y; base=field, dim=2, chunk_size=CHUNK_SIZE)

    # Split along third dimension (z)
    split_along_dim(data, out_z; base=field, dim=3, chunk_size=CHUNK_SIZE)
end

println("Done.")
