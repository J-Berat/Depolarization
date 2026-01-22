# ============================================================
# B(n) — line plot only, LaTeX fonts, log-log or linear
# ============================================================

using FITSIO;
using CairoMakie;
using Statistics;
using LaTeXStrings;
include("src/io/fits_io.jl");

# -------------------------
# USER CHOICES
# -------------------------
const DATA_DIR = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024";

const FILE_BX = "Bx.fits";
const FILE_BY = "By.fits";
const FILE_BZ = "Bz.fits";
const FILE_N = "density.fits";

const NBINS = 30;
const MIN_PER_BIN = 50;
const BINNING_IN_LOG = true;
const PLOT_LOGLOG = true;

const SAVE_FIG = true;
const OUT_FIG = "B_vs_n_line.pdf";

# -------------------------
# LOAD DATA
# -------------------------
Bx = read_FITS(joinpath(DATA_DIR, FILE_BX));
By = read_FITS(joinpath(DATA_DIR, FILE_BY));
Bz = read_FITS(joinpath(DATA_DIR, FILE_BZ));
n  = read_FITS(joinpath(DATA_DIR, FILE_N ));

# -------------------------
# BUILD |B| AND FLATTEN
# -------------------------
B  = sqrt.(Bx.^2 .+ By.^2 .+ Bz.^2);

nv = vec(n);
Bv = vec(B);

mask = isfinite.(nv) .& isfinite.(Bv) .& (nv .> 0) .& (Bv .> 0);
nv = nv[mask];
Bv = Bv[mask];

# -------------------------
# BINNING (median per bin)
# -------------------------
xbin  = BINNING_IN_LOG ? log10.(nv) : nv;
edges = range(minimum(xbin), maximum(xbin), length=NBINS+1);

x_center = Float64[];
B_med    = Float64[];

for i in 1:NBINS
    sel = (xbin .>= edges[i]) .& (xbin .< edges[i+1]);
    count(sel) < MIN_PER_BIN && continue;

    xc = 0.5 * (edges[i] + edges[i+1]);
    x_phys = BINNING_IN_LOG ? 10.0^xc : xc;

    push!(x_center, x_phys);
    push!(B_med, median(Bv[sel]));
end;

# -------------------------
# PREPARE PLOT ARRAYS
# -------------------------
xp = PLOT_LOGLOG ? log10.(x_center) : x_center;
yp = PLOT_LOGLOG ? log10.(B_med)    : B_med;

xlab = PLOT_LOGLOG ?
    L"\log_{10}\!\left(n\;[\mathrm{cm^{-3}}]\right)" :
    L"n\;[\mathrm{cm^{-3}}]";

ylab = PLOT_LOGLOG ?
    L"\log_{10}\!\left(\left|\mathbf{B}\right|\;[\mu\mathrm{G}]\right)" :
    L"\left|\mathbf{B}\right|\;[\mu\mathrm{G}]";

# -------------------------
# PLOT
# -------------------------
with_theme(theme_latexfonts()) do

    fig = Figure(size=(720, 520));

    ax = Axis(fig[1, 1],
        xlabel = xlab,
        ylabel = ylab,
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 0,   # ticks outward
        ytickalign = 0
    );

    lines!(ax, xp, yp; linewidth=2);

    SAVE_FIG && save(joinpath(DATA_DIR, OUT_FIG), fig);

    fig
end;
