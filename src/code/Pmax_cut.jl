using FITSIO
using CairoMakie
using LaTeXStrings
include("Desktop/Depolarization/src/io/fits_io.jl")

# ------------------------------------------------------------
# INPUTS
# ------------------------------------------------------------
paths = [
    "/Users/jb270005/Desktop/Depolarization_canals/d1cf05bx10rms18000nograv1024/cuts_from_transitions/LOS_y/w01_k108_160/y/Synchrotron/WithFaraday",
    "/Users/jb270005/Desktop/Depolarization_canals/d1cf05bx10rms18000nograv1024/cuts_from_transitions/LOS_y/w02_k180_198/y/Synchrotron/WithFaraday",
    "/Users/jb270005/Desktop/Depolarization_canals/d1cf05bx10rms18000nograv1024/cuts_from_transitions/LOS_y/w03_k207_219/y/Synchrotron/WithFaraday",
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday",
]
pmax_files = joinpath.(paths, "Pmax.fits")

# ------------------------------------------------------------
# USER CHOICE: COLORRANGE MODE
# ------------------------------------------------------------
# true  -> same colorrange (and thus same colorbar limits) as the full cube for all panels
# false -> individual autoscaled colorrange per panel
const USE_FULLCUBE_COLORRANGE = false

# ------------------------------------------------------------
# MARKER STYLE
# ------------------------------------------------------------
const MAP_MARK_COLOR = (:white, 0.6)
const MAP_MARK_STYLE = :dash
const MAP_MARK_LW = 2.5

const i1 = 159
const j1 = 206

# ------------------------------------------------------------
# PHYSICAL AXIS SETUP
# ------------------------------------------------------------
const LBOX_PC = 50.0
ticks_pc = 0:10:50
ticks = (collect(ticks_pc), string.(ticks_pc))

# ------------------------------------------------------------
# TYPOGRAPHY
# ------------------------------------------------------------
const TITLE_SIZE = 30
const XLABEL_SIZE = 34
const YLABEL_SIZE = 34
const XTICKLABEL_SIZE = 24
const YTICKLABEL_SIZE = 24
const CB_LABEL_SIZE = 30
const CB_TICKLABEL_SIZE = 24
const CB_WIDTH = 20

function finite_minmax(A)
    mn = Inf
    mx = -Inf
    @inbounds for v in A
        if isfinite(v)
            if v < mn; mn = v; end
            if v > mx; mx = v; end
        end
    end
    return (mn, mx)
end

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
maps = [Float32.(read_FITS(file)) for file in pmax_files]
xs = [range(0.0, LBOX_PC; length=size(M, 2)) for M in maps]
ys = [range(0.0, LBOX_PC; length=size(M, 1)) for M in maps]

titles = [
    L"10~\mathrm{pc}",
    L"3.6~\mathrm{pc}",
    L"2.4~\mathrm{pc}",
    L"\text{Full cube}",
]

# full cube is panel 4
full_colorrange = finite_minmax(maps[4])

# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------
with_theme(theme_latexfonts()) do
    fig = Figure(size = (1300, 1000))

    for (k, (M, x, y)) in enumerate(zip(maps, xs, ys))
        r = (k ≤ 2) ? 1 : 2
        c = (k % 2 == 1) ? 1 : 2

        gl = fig[r, c] = GridLayout()

        ax = Axis(
            gl[1, 1],
            title = titles[k],
            titlesize = TITLE_SIZE,

            xlabel = L"\text{Distance [pc]}",
            ylabel = L"\text{Distance [pc]}",
            xlabelsize = XLABEL_SIZE,
            ylabelsize = YLABEL_SIZE,

            xticks = ticks,
            yticks = ticks,
            xticklabelsize = XTICKLABEL_SIZE,
            yticklabelsize = YTICKLABEL_SIZE,
        )

        hm = if USE_FULLCUBE_COLORRANGE
            heatmap!(ax, x, y, M; colormap=:magma, colorrange=full_colorrange)
        else
            heatmap!(ax, x, y, M; colormap=:magma)  # autoscale per panel
        end

        xi = x[i1]
        yj = y[j1]
        vlines!(ax, [xi]; color=MAP_MARK_COLOR, linestyle=MAP_MARK_STYLE, linewidth=MAP_MARK_LW)
        hlines!(ax, [yj]; color=MAP_MARK_COLOR, linestyle=MAP_MARK_STYLE, linewidth=MAP_MARK_LW)

        Colorbar(
            gl[1, 2],
            hm;
            label = L"P_{\max}\,[K]",
            labelsize = CB_LABEL_SIZE,
            ticklabelsize = CB_TICKLABEL_SIZE,
            width = CB_WIDTH,
        )
    end

    fig
end
