using FITSIO
using Statistics
using StatsBase
using CairoMakie
using LaTeXStrings
include("Desktop/Depolarization/src/io/fits_io.jl")

# ============================================================
# USER PARAMETERS
# ============================================================
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024"

const LOS = "y"

const LBOX_PC = 50.0
const NPIX = 256
const PC_PER_PIX = LBOX_PC / NPIX

const VMIN = 0.0
const VMAX = 12.0
const COLOR_RANGE = (VMIN, VMAX)

const LABEL_SIZE = 26
const TICKLABEL_SIZE = 22
const TITLE_SIZE = 26
const LEGEND_SIZE = 22
const COLORBAR_SIZE = 24

const SAVE_PDF = true
const OUT_PDF = joinpath(SIMU_ROOT, "Pmax_deciles_histograms_LOS_$LOS.pdf")

# ============================================================
# FILE PATHS
# ============================================================

# LOS-dependent products
pmax_file = joinpath(SIMU_ROOT, LOS, "Synchrotron", "WithFaraday", "Pmax.fits")
ne_file   = joinpath(SIMU_ROOT, LOS, "Synchrotron", "ne.fits")

# Physical cubes (assumed at root)
T_file  = joinpath(SIMU_ROOT, "temperature.fits")
Bx_file = joinpath(SIMU_ROOT, "Bx.fits")
By_file = joinpath(SIMU_ROOT, "By.fits")
Bz_file = joinpath(SIMU_ROOT, "Bz.fits")

# --- Safety checks
for (lab, f) in [
    ("Pmax", pmax_file), ("ne", ne_file), ("T", T_file),
    ("Bx", Bx_file), ("By", By_file), ("Bz", Bz_file)
]
    isfile(f) || error("Missing file for $lab: $f")
end

# ============================================================
# READ FITS
# ============================================================
Pmax = read_FITS(pmax_file)

ne = read_FITS(ne_file)
T  = read_FITS(T_file)

Bx = read_FITS(Bx_file)
By = read_FITS(By_file)
Bz = read_FITS(Bz_file)

# ============================================================
# DEFINE REGIONS (1st & 10th deciles of Pmax)
# ============================================================
p10 = quantile(vec(Pmax), 0.10)
p90 = quantile(vec(Pmax), 0.90)

mask_low  = Pmax .<= p10
mask_high = Pmax .>= p90

# ============================================================
# EXTRACTION FUNCTION (apply 2D mask along LOS)
# ============================================================
function extract_masked_los(cube, mask2D)
    nx, ny, nz = size(cube)
    vals = Float64[]
    for i in 1:nx, j in 1:ny
        if mask2D[i, j]
            append!(vals, cube[i, j, :])
        end
    end
    return vals
end

# ============================================================
# EXTRACT VALUES
# ============================================================
Bmod = 1000 .* sqrt.(Bx.^2 .+ By.^2 .+ Bz.^2)   # µG

ne_low  = extract_masked_los(ne,   mask_low)
T_low   = extract_masked_los(T,    mask_low)
B_low   = extract_masked_los(Bmod, mask_low)

ne_high = extract_masked_los(ne,   mask_high)
T_high  = extract_masked_los(T,    mask_high)
B_high  = extract_masked_los(Bmod, mask_high)

# ============================================================
# HISTOGRAMS (COUNTS)
# ============================================================
nbins = 80

h_ne_low  = fit(Histogram, ne_low,  nbins=nbins, closed=:left)
h_T_low   = fit(Histogram, T_low,   nbins=nbins, closed=:left)
h_B_low   = fit(Histogram, B_low,   nbins=nbins, closed=:left)

h_ne_high = fit(Histogram, ne_high, nbins=nbins, closed=:left)
h_T_high  = fit(Histogram, T_high,  nbins=nbins, closed=:left)
h_B_high  = fit(Histogram, B_high,  nbins=nbins, closed=:left)

# ============================================================
# TICKS IN pc (0–50 by 10)
# ============================================================
pc_ticks  = 0:10:50
pix_ticks = pc_ticks ./ PC_PER_PIX

# ============================================================
# PLOT + EXPORT
# ============================================================
with_theme(theme_latexfonts()) do

    fig = Figure(size = (1300, 560))

    # --------------------------------------------------------
    # (A) Pmax heatmap + contours
    # --------------------------------------------------------
    axH = Axis(
        fig[1:3, 1],
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        xlabelsize = LABEL_SIZE,
        ylabelsize = LABEL_SIZE,
        xticks = (pix_ticks, string.(pc_ticks)),
        yticks = (pix_ticks, string.(pc_ticks)),
        xticklabelsize = TICKLABEL_SIZE,
        yticklabelsize = TICKLABEL_SIZE,
        xgridvisible = false,
        ygridvisible = false,
    )

    hm = heatmap!(
        axH,
        Pmax;
        colormap   = :magma,
        colorrange = COLOR_RANGE,
    )

    contour!(axH, mask_low;  levels=[0.5], color=:black, linewidth=3)
    contour!(axH, mask_high; levels=[0.5], color=:red,   linewidth=3)

    Colorbar(
        fig[1:3, 2],
        hm,
        label = L"$P_{\mathrm{max}}\,[\mathrm{K}]$",
        labelsize = COLORBAR_SIZE,
        ticklabelsize = TICKLABEL_SIZE,
    )

    # --------------------------------------------------------
    # (B) HISTOGRAMS — counts
    # --------------------------------------------------------
    ax1 = Axis(
        fig[1, 3],
        xlabel = L"$n_e\,[\mathrm{cm}^{-3}]$",
        ylabel = "Counts",
        title  = "Decile comparison",
        titlesize = TITLE_SIZE,
        titlefont = :regular,
        xlabelsize = LABEL_SIZE,
        ylabelsize = LABEL_SIZE,
        xticklabelsize = TICKLABEL_SIZE,
        yticklabelsize = TICKLABEL_SIZE,
        xgridvisible = false,
        ygridvisible = false,
    )

    ax2 = Axis(
        fig[2, 3],
        xlabel = L"$T\,[\mathrm{K}]$",
        ylabel = "Counts",
        titlefont = :regular,
        xlabelsize = LABEL_SIZE,
        ylabelsize = LABEL_SIZE,
        xticklabelsize = TICKLABEL_SIZE,
        yticklabelsize = TICKLABEL_SIZE,
        xgridvisible = false,
        ygridvisible = false,
    )

    ax3 = Axis(
        fig[3, 3],
        xlabel = L"$|\mathbf{B}|\,[\mu\mathrm{G}]$",
        ylabel = "Counts",
        titlefont = :regular,
        xlabelsize = LABEL_SIZE,
        ylabelsize = LABEL_SIZE,
        xticklabelsize = TICKLABEL_SIZE,
        yticklabelsize = TICKLABEL_SIZE,
        xgridvisible = false,
        ygridvisible = false,
    )

    stairs!(ax1, h_ne_low.edges[1][1:end-1],  h_ne_low.weights,  color=:black, linewidth=3)
    stairs!(ax1, h_ne_high.edges[1][1:end-1], h_ne_high.weights, color=:red,   linewidth=3)

    stairs!(ax2, h_T_low.edges[1][1:end-1],   h_T_low.weights,   color=:black, linewidth=3)
    stairs!(ax2, h_T_high.edges[1][1:end-1],  h_T_high.weights,  color=:red,   linewidth=3)

    stairs!(ax3, h_B_low.edges[1][1:end-1],   h_B_low.weights,   color=:black, linewidth=3)
    stairs!(ax3, h_B_high.edges[1][1:end-1],  h_B_high.weights,  color=:red,   linewidth=3)

    # legend
    l1 = lines!(ax1, [NaN], [NaN], color=:black, linewidth=3)
    l2 = lines!(ax1, [NaN], [NaN], color=:red,   linewidth=3)

    axislegend(
        ax1,
        [l1, l2],
        ["1st decile", "10th decile"],
        position = :rt,
        framevisible = false,
        labelsize = LEGEND_SIZE,
    )

    # --------------------------------------------------------
    # Layout
    # --------------------------------------------------------
    rowgap!(fig.layout, 12)
    colgap!(fig.layout, 16)

    colsize!(fig.layout, 1, Relative(0.42))
    colsize!(fig.layout, 2, Auto(0.0))
    colsize!(fig.layout, 3, Relative(0.58))

    # --------------------------------------------------------
    # EXPORT
    # --------------------------------------------------------
    if SAVE_PDF
        save(OUT_PDF, fig; pt_per_unit = 1)
        @info "PDF saved to $OUT_PDF"
    end

    fig
end
