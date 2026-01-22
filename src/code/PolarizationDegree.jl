using FITSIO
using CairoMakie
using Statistics
using LaTeXStrings
include(joinpath(@__DIR__, "../constants.jl"))
include(joinpath(@__DIR__, "../io/fits_io.jl"))

# ------------------------------------------------------------
# USER CHOICES
# ------------------------------------------------------------
const SIMU_ROOT = DepolarizationConstants.PolarizationDegree.SIMU_ROOT
const SIMU_NAME = DepolarizationConstants.PolarizationDegree.SIMU_NAME
const LOS = DepolarizationConstants.PolarizationDegree.LOS

const NPIX = DepolarizationConstants.PolarizationDegree.NPIX
const LBOX_PC = DepolarizationConstants.PolarizationDegree.LBOX_PC

const CMAP = DepolarizationConstants.PolarizationDegree.CMAP
const VMIN = DepolarizationConstants.PolarizationDegree.VMIN
const VMAX = DepolarizationConstants.PolarizationDegree.VMAX

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
const SYNCHRO_DIR = DepolarizationConstants.PolarizationDegree.SYNCHRO_DIR
const I_FILE = DepolarizationConstants.PolarizationDegree.I_FILE
const Q_FILE = DepolarizationConstants.PolarizationDegree.Q_FILE
const U_FILE = DepolarizationConstants.PolarizationDegree.U_FILE
const OUT_DIR = DepolarizationConstants.PolarizationDegree.OUT_DIR
const OUT_PDF = DepolarizationConstants.PolarizationDegree.OUT_PDF

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
