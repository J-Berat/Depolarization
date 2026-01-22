using FITSIO
using CairoMakie
using Statistics
using LaTeXStrings
include("src/io/fits_io.jl")

# ------------------------------------------------------------
# USER CHOICES
# ------------------------------------------------------------
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"
const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS = "y"

const NPIX = 256
const LBOX_PC = 50.0

const CMAP = :viridis
const VMIN = 0.0
const VMAX = 1.0

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
const SYNCHRO_DIR = joinpath(
    SIMU_ROOT, SIMU_NAME, LOS, "Synchrotron", "WithFaraday"
)

const I_FILE = joinpath(SYNCHRO_DIR, "I.fits")
const Q_FILE = joinpath(SYNCHRO_DIR, "Q.fits")
const U_FILE = joinpath(SYNCHRO_DIR, "U.fits")

const OUT_DIR = joinpath(SYNCHRO_DIR, "Plots")
const OUT_PDF = joinpath(OUT_DIR, "p_degree_$(LOS).pdf")

mkpath(OUT_DIR)

# ------------------------------------------------------------
# DEGREE OF POLARIZATION
# p = sqrt(Q^2 + U^2) / I
# ------------------------------------------------------------
function degree_polarization(I, Q, U)
    @assert size(I) == size(Q) == size(U)
    P = sqrt.(Q.^2 .+ U.^2)
    return P ./ I
end

# ------------------------------------------------------------
# Reduce cube → 2D map
# Convention: (nchan, ny, nx) → max over channels
# ------------------------------------------------------------
function to_map(A)
    if ndims(A) == 2
        return A
    elseif ndims(A) == 3
        return dropdims(maximum(A; dims=1), dims=1)
    else
        error("Unsupported array dimensions")
    end
end

# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------
function plot_p_map(p2)
    fig = Figure(size=(950, 800))

    ax = Axis(
        fig[1, 1],
        title  = "LOS: $(LOS)",
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]"
    )

    x = range(0, LBOX_PC, length=size(p2, 2))
    y = range(0, LBOX_PC, length=size(p2, 1))

    hm = heatmap!(
        ax, x, y, p2;
        colormap = CMAP,
        colorrange = (VMIN, VMAX)
    )

    Colorbar(fig[1, 2], hm, label = L"p = \sqrt{Q^2+U^2}/I")

    ax.xticks = 0:10:LBOX_PC
    ax.yticks = 0:10:LBOX_PC

    fig
end

# ------------------------------------------------------------
# RUN
# ------------------------------------------------------------
I = read_FITS(I_FILE)
Q = read_FITS(Q_FILE)
U = read_FITS(U_FILE)

p  = degree_polarization(I, Q, U)
p2 = to_map(p)

fig = plot_p_map(p2)
save(OUT_PDF, fig)

@info "PDF saved → $OUT_PDF"
