using FITSIO
using CairoMakie
using Colors
using LaTeXStrings
include("src/io/fits_io.jl")

# =========================
# USER PARAMETERS
# =========================
LOS  = "z"   # <--- change here: "x", "y", or "z"
simu = "d1cf05bx10rms18000nograv1024"

# Physical size (pc) and pixel size
Lbox_pc = 50.0
npix = 256
pix_to_pc = Lbox_pc / npix

# Axis ticks in pc (0 -> 50)
pc_ticks = 0:10:50
pc_ticklabels = string.(pc_ticks)
pix_ticks = pc_ticks ./ pix_to_pc  # positions in pixel coordinates

# =========================
# Colormap (custom plasma)
# =========================
custom_plasma = cgrad(
    [
        colorant"#3B0F70",
        colorant"#8C2981",
        colorant"#E16462",
        colorant"#FCA636",
        colorant"#FCDEA2",
    ],
    categorical = false
)

# =========================
# File paths
# =========================
base_split = joinpath(
    "/Users/jb270005/Desktop/Depolarization_canals",
    simu,
    "split_$(LOS)"
)

chunk_files = [
    joinpath(base_split, "chunk01", LOS, "Synchrotron", "WithFaraday", "Pmax.fits"),
    joinpath(base_split, "chunk02", LOS, "Synchrotron", "WithFaraday", "Pmax.fits"),
    joinpath(base_split, "chunk03", LOS, "Synchrotron", "WithFaraday", "Pmax.fits"),
    joinpath(base_split, "chunk04", LOS, "Synchrotron", "WithFaraday", "Pmax.fits"),
]

# If you also have the unsplit full Pmax in simu_RAMSES:
single_file = joinpath(
    "/Users/jb270005/Desktop/simu_RAMSES",
    simu, LOS, "Synchrotron", "WithFaraday", "Pmax.fits"
)

# =========================
# Titles (physical depth intervals)
# =========================
titles_pc = [
    L"0\text{-}12.5~\mathrm{pc}",
    L"12.5\text{-}25~\mathrm{pc}",
    L"25\text{-}37.5~\mathrm{pc}",
    L"37.5\text{-}50~\mathrm{pc}",
]

# =========================
# Read data
# =========================
P_chunks = [read_FITS(f) for f in chunk_files]
P_single = read_FITS(single_file)

# =========================
# FIGURES: 4 independent chunk plots
# =========================
with_theme(theme_latexfonts()) do
    for i in 1:4
        fig = Figure(size = (850, 700))

        ax = Axis(
            fig[1, 1],
            title = titles_pc[i],
            xlabel = L"\textrm{Distance}~[\mathrm{pc}]",
            ylabel = L"\textrm{Distance}~[\mathrm{pc}]",
            xgridvisible = false,
            ygridvisible = false,
            backgroundcolor = :white,
            xticks = (pix_ticks, pc_ticklabels),
            yticks = (pix_ticks, pc_ticklabels),
            titlesize = 28,
            xticklabelsize = 22,
            yticklabelsize = 22,
            xlabelsize = 24,
            ylabelsize = 24,
        )

        hm = heatmap!(ax, P_chunks[i]; colormap = custom_plasma)
        Colorbar(fig[1, 2], hm; label = L"[K]", labelsize = 24, ticklabelsize = 20)

        colgap!(fig.layout, 20)

        save(
            "/Users/jb270005/Desktop/Depolarization_canals/Pmax_chunk$(lpad(i,2,'0'))_split_$(LOS).pdf",
            fig;
            pt_per_unit = 1
        )
    end
end

# =========================
# FIGURE: single Pmax
# =========================
with_theme(theme_latexfonts()) do
    fig = Figure(size = (850, 700))

    ax = Axis(
        fig[1, 1],
        xlabel = L"\textrm{Distance}~[\mathrm{pc}]",
        ylabel = L"\textrm{Distance}~[\mathrm{pc}]",
        xgridvisible = false,
        ygridvisible = false,
        backgroundcolor = :white,
        xticks = (pix_ticks, pc_ticklabels),
        yticks = (pix_ticks, pc_ticklabels),
        xticklabelsize = 22,
        yticklabelsize = 22,
        xlabelsize = 24,
        ylabelsize = 24,
    )

    hm = heatmap!(ax, P_single; colormap = custom_plasma)
    Colorbar(fig[1, 2], hm; label = L"[K]", labelsize = 24, ticklabelsize = 20)

    colgap!(fig.layout, 20)

    save(
        "/Users/jb270005/Desktop/Depolarization_canals/Pmax_single_$(LOS).pdf",
        fig;
        pt_per_unit = 1
    )

    fig
end
