using CairoMakie
using Makie
using Statistics
using StatsBase
using FITSIO

# ============================================================
# USER CHOICES
# ============================================================
SimuList = [
    "d1cf00bx10rms18000nograv1024",
    "d1cf02bx10rms09000nograv1024",
    "d1cf02bx10rms18000nograv1024",
    "d1cf02bx10rms36000nograv1024",
    "d1cf02bx10rms72000nograv1024",
    "d1cf05bx10rms09000nograv1024",
    "d1cf05bx10rms18000nograv1024",
    "d1cf05bx10rms36000nograv1024",
    "d1cf05bx10rms72000nograv1024",
    "d1cf08bx10rms09000nograv1024",
    "d1cf08bx10rms18000nograv1024",
    "d1cf08bx10rms36000nograv1024",
    "d1cf08bx10rms72000nograv1024",
    "d1cf10bx10rms18000nograv1024"
]

SIMU_IDX = 7          # <<< choose simulation
LOS      = "y"        # <<< DEFAULT LOS
@assert LOS == "x" || LOS == "y"

SIMU = SimuList[SIMU_IDX]
BASE = "/Users/jb270005/Desktop/simu_RAMSES"

# ============================================================
# PARAMETERS
# ============================================================
Lbox_pc        = 50.0
mask_q         = 0.01
canal_decile   = 0.10
use_same_thr   = true

# Ratio
R_colorrange = (0.0, 1.0)

# >>> USER-DEFINED Pmax limits (SET HERE) <<<
PnoRM_VMIN = 0.0      # [K]
PnoRM_VMAX = 20  # [K]
PRM_VMIN   = 0.0      # [K]
PRM_VMAX   = 10   # [K]

# Font sizes
FS_TITLE = 30
FS_LABEL = 30
FS_TICK  = 24
FS_CBAR  = 26

# ============================================================
# PATHS
# ============================================================
p_withRM = joinpath(BASE, SIMU, LOS, "Synchrotron", "WithFaraday", "Pmax.fits")
p_noRM   = joinpath(BASE, SIMU, LOS, "Synchrotron", "noFaraday", "Pnumax.fits")

# ============================================================
# FITS I/O
# ============================================================
readfits(path::AbstractString) = Float32.(FITS(path, "r") do f
    read(f[1])
end)

# ============================================================
# LOAD
# ============================================================
P0 = readfits(p_noRM)    # PnoRMmax
P1 = readfits(p_withRM)  # PRMmax
@assert size(P0) == size(P1)

ny, nx = size(P0)

x_pc = range(0.0, Lbox_pc; length=nx)
y_pc = range(0.0, Lbox_pc; length=ny)
ticks_pc = 0:10:50
ticklabs = string.(ticks_pc)

# ============================================================
# MAPS
# ============================================================
# Ratio
eps = quantile(vec(P0), mask_q)
R = fill(NaN32, size(P0))
@inbounds for j in 1:ny, i in 1:nx
    if P0[j, i] ≥ eps
        R[j, i] = P1[j, i] / P0[j, i]
    end
end

# Canals
thr0 = quantile(vec(P0), canal_decile)
thr1 = quantile(vec(P1), canal_decile)

C0 = P0 .< thr0
C1 = use_same_thr ? (P1 .< thr0) : (P1 .< thr1)

C0f = Float32.(C0)
C1f = Float32.(C1)

cmap01 = cgrad([:white, :black]; categorical=true)

function heatmap_opt!(ax, x, y, A; colormap, colorrange=nothing)
    colorrange === nothing ?
        heatmap!(ax, x, y, A; colormap) :
        heatmap!(ax, x, y, A; colormap, colorrange)
end

# ============================================================
# OUTPUT FILES
# ============================================================
out_ratio  = "Pmax_ratio_RM_over_noRM__$(SIMU)__LOS_$(LOS).pdf"
out_canals = "Pmax_canals_C0_C1__$(SIMU)__LOS_$(LOS).pdf"
out_pmax   = "Pmax_maps_noRM_vs_RM__$(SIMU)__LOS_$(LOS).pdf"

# ============================================================
# PLOTTING
# ============================================================
with_theme(theme_latexfonts()) do

    # ---------------------- Figure 1 : Ratio
    figR = Figure(size = (980, 820))
    axR = Axis(
        figR[1,1],
        title = L"$P_{\max}^{\mathrm{RM}} / P_{\max}^{\mathrm{noRM}}$",
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        aspect = DataAspect(),
        xticks = (ticks_pc, ticklabs),
        yticks = (ticks_pc, ticklabs),
        titlesize = FS_TITLE,
        xlabelsize = FS_LABEL,
        ylabelsize = FS_LABEL,
        xticklabelsize = FS_TICK,
        yticklabelsize = FS_TICK,
    )

    hmR = heatmap_opt!(axR, x_pc, y_pc, R; colormap=:viridis, colorrange=R_colorrange)
    Colorbar(figR[1,2], hmR,
        label = L"$P_{\max}^{\mathrm{RM}} / P_{\max}^{\mathrm{noRM}}$",
        labelsize = FS_CBAR,
        ticklabelsize = FS_TICK)

    save(out_ratio, figR)
    display(figR)

    # ---------------------- Figure 2 : Canals
    figC = Figure(size = (1400, 600))

    axC0 = Axis(
        figC[1,1],
        title = L"$C_0:\; P_{\max}^{\mathrm{noRM}} < Q_{0.1}$",
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        aspect = DataAspect(),
        xticks = (ticks_pc, ticklabs),
        yticks = (ticks_pc, ticklabs),
        titlesize = FS_TITLE,
        xlabelsize = FS_LABEL,
        ylabelsize = FS_LABEL,
        xticklabelsize = FS_TICK,
        yticklabelsize = FS_TICK,
    )

    axC1 = Axis(
        figC[1,2],
        title = use_same_thr ?
            L"$C_1:\; P_{\max}^{\mathrm{RM}} < Q_{0.1}(P_{\max}^{\mathrm{noRM}})$" :
            L"$C_1:\; P_{\max}^{\mathrm{RM}} < Q_{0.1}$",
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        aspect = DataAspect(),
        xticks = (ticks_pc, ticklabs),
        yticks = (ticks_pc, ticklabs),
        titlesize = FS_TITLE,
        xlabelsize = FS_LABEL,
        ylabelsize = FS_LABEL,
        xticklabelsize = FS_TICK,
        yticklabelsize = FS_TICK,
    )

    hmC0 = heatmap!(axC0, x_pc, y_pc, C0f; colormap=cmap01, colorrange=(0,1))
    hmC1 = heatmap!(axC1, x_pc, y_pc, C1f; colormap=cmap01, colorrange=(0,1))

    Colorbar(figC[1,3], hmC1,
        label = L"\mathrm{mask}",
        ticks = ([0,1], ["0","1"]),
        labelsize = FS_CBAR,
        ticklabelsize = FS_TICK)

    save(out_canals, figC)
    display(figC)

    # ---------------------- Figure 3 : Pmax maps (NO TITLES)
    figP = Figure(size = (1600, 700))

    axP0 = Axis(
        figP[1,1],
        title = "",                      # no title
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        aspect = DataAspect(),
        xticks = (ticks_pc, ticklabs),
        yticks = (ticks_pc, ticklabs),
        xlabelsize = FS_LABEL,
        ylabelsize = FS_LABEL,
        xticklabelsize = FS_TICK,
        yticklabelsize = FS_TICK,
    )

    axP1 = Axis(
        figP[1,3],
        title = "",                      # no title
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        aspect = DataAspect(),
        xticks = (ticks_pc, ticklabs),
        yticks = (ticks_pc, ticklabs),
        xlabelsize = FS_LABEL,
        ylabelsize = FS_LABEL,
        xticklabelsize = FS_TICK,
        yticklabelsize = FS_TICK,
    )

    hmP0 = heatmap!(
        axP0, x_pc, y_pc, P0;
        colormap=:magma,
        colorrange=(PnoRM_VMIN, PnoRM_VMAX),
    )

    hmP1 = heatmap!(
        axP1, x_pc, y_pc, P1;
        colormap=:magma,
        colorrange=(PRM_VMIN, PRM_VMAX),
    )

    Colorbar(figP[1,2], hmP0,
        label = L"$P_{\mathrm{noRM,max}}\;[\mathrm{K}]$",
        labelsize = FS_CBAR,
        ticklabelsize = FS_TICK)

    Colorbar(figP[1,4], hmP1,
        label = L"$P_{\mathrm{RM,max}}\;[\mathrm{K}]$",
        labelsize = FS_CBAR,
        ticklabelsize = FS_TICK)

    save(out_pmax, figP)
    display(figP)
end
