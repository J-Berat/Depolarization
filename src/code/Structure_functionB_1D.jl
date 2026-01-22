using FITSIO
using Statistics
using CairoMakie
using LaTeXStrings
include("Desktop/Depolarization/src/io/fits_io.jl")

# ------------------------------------------------------------
# Extract B_parallel LOS signal (in μG)
# cube axes are (x, y, z)
# ------------------------------------------------------------
function los_signal_Bparallel_uG(Bx, By, Bz, LOS::AbstractString, i_map::Int, j_map::Int)
    if LOS == "x"
        # depth=x, map=(y,z)
        v = vec(@view Bx[:, i_map, j_map])
    elseif LOS == "y"
        # depth=y, map=(x,z)
        v = vec(@view By[i_map, :, j_map])
    elseif LOS == "z"
        # depth=z, map=(x,y)
        v = vec(@view Bz[i_map, j_map, :])
    else
        error("LOS must be \"x\", \"y\", or \"z\"")
    end
    return 1000.0 .* v  # μG
end

# ------------------------------------------------------------
# 1D structure function: S_p(ℓ)=<|v(k+ℓ)-v(k)|^p>
# ------------------------------------------------------------
function structure_function_1d(v::AbstractVector; p::Int=2, maxlag=nothing)
    n = length(v)
    maxlag === nothing && (maxlag = n ÷ 2)
    maxlag = min(maxlag, n - 1)

    lags = collect(1:maxlag)
    Sp = Vector{Float64}(undef, maxlag)

    for ℓ in lags
        v1 = @view v[(1+ℓ):n]
        v0 = @view v[1:(n-ℓ)]
        dv = v1 .- v0
        Sp[ℓ] = mean(abs.(dv) .^ p)
    end

    return lags, Sp
end

# ------------------------------------------------------------
# USER PARAMETERS
# ------------------------------------------------------------
SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024"

LOS   = "y"      # "x", "y", or "z"
i_map = 30      # pixel in plane-of-sky
j_map = 128

# NOTE pixel meaning:
# LOS="x" => (i_map,j_map) = (y,z)
# LOS="y" => (i_map,j_map) = (x,z)
# LOS="z" => (i_map,j_map) = (x,y)

# Structure function settings
p_order = 2
maxlag  = 120

# Physical sizes (pc)
Lmap_pc = 50.0                 # for Pmax axes (0..50 pc)
tick_pc = 0:10:50
Lbox_pc = 50.0                 # LOS depth size in pc (set to nothing to keep ℓ in cells)
use_pc_los = true

# Pmax colormap + range (USER)
VMIN = 0.0
VMAX = 12.0
CMAP = :magma

# Font sizes
labsize  = 26
ticksize = 20

# ------------------------------------------------------------
# PATHS (auto from LOS)
# ------------------------------------------------------------
Bx_path = joinpath(SIMU_ROOT, "Bx.fits")
By_path = joinpath(SIMU_ROOT, "By.fits")
Bz_path = joinpath(SIMU_ROOT, "Bz.fits")
PMAX_PATH = joinpath(SIMU_ROOT, LOS, "Synchrotron", "WithFaraday", "Pmax.fits")

# ------------------------------------------------------------
# READ DATA
# ------------------------------------------------------------
Bx = read_FITS(Bx_path)
By = read_FITS(By_path)
Bz = read_FITS(Bz_path)
@assert ndims(Bx) == 3 && size(Bx) == size(By) == size(Bz)

Pmax = read_FITS(PMAX_PATH)
@assert ndims(Pmax) == 2 "Pmax must be a 2D map"

# ------------------------------------------------------------
# COMPUTE STRUCTURE FUNCTION
# ------------------------------------------------------------
Bpar = los_signal_Bparallel_uG(Bx, By, Bz, LOS, i_map, j_map)
lags, Sp = structure_function_1d(Bpar; p=p_order, maxlag=maxlag)

n_depth = length(Bpar)
dℓ = (Lbox_pc === nothing) ? 1.0 : (Lbox_pc / n_depth)
xvals = (use_pc_los && Lbox_pc !== nothing) ? (lags .* dℓ) : lags
xlabel_sf = (use_pc_los && Lbox_pc !== nothing) ? L"$\ell\ \mathrm{[pc]}$" : L"$\ell\ \mathrm{[cell]}$"

# ------------------------------------------------------------
# COORDINATES FOR PMAX AXES (Distance [pc])
# ------------------------------------------------------------
Nx, Ny = size(Pmax)
xpc = range(0.0, Lmap_pc; length=Nx)
ypc = range(0.0, Lmap_pc; length=Ny)

# Crosshair at chosen pixel (assumes Pmax indexing matches your i_map/j_map convention)
x0 = xpc[i_map]
y0 = ypc[j_map]
cross_col = RGBAf(1, 1, 1, 0.8)

# ------------------------------------------------------------
# PLOT (LaTeX theme)
# ------------------------------------------------------------
with_theme(theme_latexfonts()) do
    fig = Figure(size=(1500, 520))

    # Left: Pmax heatmap
    axL = Axis(
        fig[1, 1],
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        xlabelsize = labsize,
        ylabelsize = labsize,
        xticklabelsize = ticksize,
        yticklabelsize = ticksize,
        xgridvisible = false,
        ygridvisible = false,
        xticks = (collect(tick_pc), string.(collect(tick_pc))),
        yticks = (collect(tick_pc), string.(collect(tick_pc)))
    )

    hm = heatmap!(
        axL, xpc, ypc, Pmax;
        colormap = CMAP,
        colorrange = (VMIN, VMAX)
    )

    vlines!(axL, [x0]; color=cross_col, linewidth=2.0, linestyle=:dash)
    hlines!(axL, [y0]; color=cross_col, linewidth=2.0, linestyle=:dash)
    xlims!(axL, 0, Lmap_pc)
    ylims!(axL, 0, Lmap_pc)

    # Colorbar for Pmax
    Colorbar(
        fig[1, 2],
        hm;
        label = L"$P_{\max}$",
        labelsize = 22,
        ticklabelsize = 18
    )

    # Right: structure function
    axR = Axis(
        fig[1, 3],
        xlabel = xlabel_sf,
        ylabel = L"$S_2(\ell)\ [\mu\mathrm{G}^2]$",
        xlabelsize = labsize,
        ylabelsize = labsize,
        xticklabelsize = ticksize,
        yticklabelsize = ticksize,
        xgridvisible = false,
        ygridvisible = false
    )

    lines!(axR, xvals, Sp, linewidth=2.8)
    scatter!(axR, xvals, Sp, markersize=6)

    colgap!(fig.layout, 1, 14)
    colgap!(fig.layout, 2, 28)

    fig
end
