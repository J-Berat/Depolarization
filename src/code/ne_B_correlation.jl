using FITSIO
using Statistics
using CairoMakie
using LaTeXStrings

# ============================================================
# FITS READER 
# ============================================================
function read_FITS_file(file)
    read(FITS(file)[1])
end

# ============================================================
# PATHS & FILES
# ============================================================
const DATA_DIR = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024"

const FILE_BX = "Bx.fits"
const FILE_BY = "By.fits"
const FILE_BZ = "Bz.fits"

const NE_SUBDIR = "Synchrotron"
const FILE_NE   = "ne.fits"

# ============================================================
# PARAMETERS
# ============================================================
const LOS = "x"          # "x", "y", or "z"
const LOS_DIM = 1        # dimension index corresponding to LOS

const WHAT = :rF         # :CF or :rF
const VMIN = nothing     # e.g. -0.3 or nothing
const VMAX = nothing     # e.g.  0.3 or nothing

const LBOX_PC = 50.0
const OUT = "FaradayCov_$(WHAT)_LOS_$(LOS).pdf"

# ============================================================
# HELPERS
# ============================================================
@inline function Bparallel(Bx, By, Bz, los::AbstractString)
    los == "x" && return Bx
    los == "y" && return By
    los == "z" && return Bz
    error("los must be \"x\", \"y\", or \"z\"")
end

@inline function ne_path(base::AbstractString, los::AbstractString)
    joinpath(base, los, NE_SUBDIR, FILE_NE)
end

"""
    faraday_cov(ne, Bpar; dims=:)

Returns:
  C_F = ⟨ne B∥⟩ − ⟨ne⟩⟨B∥⟩
  r_F = C_F / (σ_ne σ_B∥)
"""
function faraday_cov(ne::AbstractArray, Bpar::AbstractArray; dims=:)
    @assert size(ne) == size(Bpar)

    μ_ne  = mean(ne; dims=dims)
    μ_B   = mean(Bpar; dims=dims)
    μ_neB = mean(ne .* Bpar; dims=dims)

    C_F = μ_neB .- μ_ne .* μ_B

    σ_ne = std(ne; dims=dims)
    σ_B  = std(Bpar; dims=dims)
    r_F  = C_F ./ (σ_ne .* σ_B)

    return C_F, r_F
end

# robust extrema ignoring NaN / Inf
function finite_extrema(A)
    vmin =  Inf
    vmax = -Inf
    @inbounds for x in A
        if isfinite(x)
            vmin = min(vmin, x)
            vmax = max(vmax, x)
        end
    end
    return vmin, vmax
end

# ============================================================
# LOAD DATA
# ============================================================
Bx = read_FITS_file(joinpath(DATA_DIR, FILE_BX))
By = read_FITS_file(joinpath(DATA_DIR, FILE_BY))
Bz = read_FITS_file(joinpath(DATA_DIR, FILE_BZ))
ne = read_FITS_file(ne_path(DATA_DIR, LOS))

@assert size(Bx) == size(By) == size(Bz) == size(ne) "Cube size mismatch"

Bpar = Bparallel(Bx, By, Bz, LOS)

# ============================================================
# FARADAY-RELEVANT COVARIANCE (2D MAP)
# ============================================================
C_map, r_map = faraday_cov(ne, Bpar; dims=LOS_DIM)

C_map = dropdims(C_map; dims=LOS_DIM)
r_map = dropdims(r_map; dims=LOS_DIM)

Z = (WHAT == :CF) ? C_map : r_map

# ============================================================
# AXES IN pc
# ============================================================
nx, ny = size(Z)
x_pc = range(0, LBOX_PC; length=nx)
y_pc = range(0, LBOX_PC; length=ny)

ticks = 0:10:LBOX_PC

# ============================================================
# COLORRANGE (ROBUST)
# ============================================================
vmin_auto, vmax_auto = finite_extrema(Z)

if WHAT == :rF && VMIN === nothing && VMAX === nothing
    m = max(abs(vmin_auto), abs(vmax_auto))
    cr = (Float32(-m), Float32(m))      # symmetric, recommended
else
    vmin = (VMIN === nothing) ? vmin_auto : VMIN
    vmax = (VMAX === nothing) ? vmax_auto : VMAX
    cr = (Float32(vmin), Float32(vmax))
end

# ============================================================
# PLOT (LATEX THEME, BIG LABELS, NO TITLE)
# ============================================================
with_theme(theme_latexfonts()) do

    fig = Figure(size=(900, 800))

    ax = Axis(
        fig[1, 1],
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",

        xlabelsize = 28,
        ylabelsize = 28,

        xticklabelsize = 22,
        yticklabelsize = 22,

        xticks = (ticks, string.(Int.(ticks))),
        yticks = (ticks, string.(Int.(ticks))),

        xgridvisible = false,
        ygridvisible = false,
    )

    hm = heatmap!(
        ax,
        x_pc,
        y_pc,
        Z;
        colorrange = cr,
        colormap = :magma
    )

    Colorbar(
        fig[1, 2],
        hm;
        label = (WHAT == :CF ? L"C_F" : L"r_F"),
        labelsize = 26,
        ticklabelsize = 22,
        width = 22
    )

    display(fig)      # <<< interactive display
    save(OUT, fig)    # <<< PDF export
end

println("Saved → $OUT")