using CairoMakie
using LaTeXStrings
using Statistics

include("fits_io.jl")

# ============================================================
# PARAMETERS
# ============================================================
simu = "d1cf05bx10rms18000nograv1024"
LOS  = "x"

base_dir = "/Users/jb270005/Desktop/simu_RAMSES"
syn_dir  = joinpath(base_dir, simu, LOS, "Synchrotron", "WithFaraday")

pmax_file = joinpath(syn_dir, "Pmax.fits")
tnu_file  = joinpath(syn_dir, "Tnu.fits")

out_pdf = joinpath(syn_dir, "Pmax_and_beta.pdf")

# Geometry
Lbox_pc = 50.0

# Frequencies
ν0 = 120e6          # Hz
Δν = 98e3           # Hz

# Color ranges (nothing = automatic)
VMIN_PMAX = 0.0
VMAX_PMAX = 12.0
VMIN_BETA = nothing
VMAX_BETA = nothing

# Colormap
CMAP = :magma

# ============================================================
# COMPUTE β = d log T / d log ν
# ============================================================
function beta_map_from_cube(Tcube::Array{Float32,3}, nu::Vector{Float64})
    nx, ny, nν = size(Tcube)
    @assert length(nu) == nν "Frequency dimension mismatch"

    x = log.(nu)
    mx = mean(x)
    varx = sum((x .- mx).^2)

    β = Array{Float32}(undef, nx, ny)

    @inbounds for i in 1:nx, j in 1:ny
        y = log.(Float64.(view(Tcube, i, j, :)))
        my = mean(y)
        covxy = sum((x .- mx) .* (y .- my))
        β[i, j] = Float32(covxy / varx)
    end

    return β
end

# ============================================================
# HEATMAP WRAPPER 
# ============================================================
function hm!(ax, x, y, Z; cmap, vmin, vmax)
    if vmin === nothing || vmax === nothing
        heatmap!(ax, x, y, Z; colormap=cmap)
    else
        heatmap!(ax, x, y, Z;
            colormap=cmap,
            colorrange=(Float32(vmin), Float32(vmax))
        )
    end
end

Pmap  = Float32.(read_FITS(pmax_file))
Tcube = Float32.(read_FITS(tnu_file))

@assert ndims(Tcube) == 3 "Tnu.fits must be a 3D cube"

nx, ny, nν = size(Tcube)
nu = Float64.(ν0 .+ (0:nν-1) .* Δν)

println("Tcube size = ", size(Tcube))
println("ν range    = ", nu[1]/1e6, " → ", nu[end]/1e6, " MHz")

beta_map = beta_map_from_cube(Tcube, nu)

# Spatial axes
ticks_pc = 0:10:50
x = range(0, Lbox_pc; length=nx)
y = range(0, Lbox_pc; length=ny)

# ============================================================
# FIGURE
# ============================================================
with_theme(theme_latexfonts()) do
    fig = Figure(size=(1200, 550))

    labsize  = 26
    ticksize = 22

    # --- Pmax ---
    ax1 = Axis(fig[1,1],
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        title  = "",
        xticks = (ticks_pc, string.(ticks_pc)),
        yticks = (ticks_pc, string.(ticks_pc)),
        xgridvisible = false,
        ygridvisible = false
    )
    ax1.xlabelsize = labsize
    ax1.ylabelsize = labsize
    ax1.xticklabelsize = ticksize
    ax1.yticklabelsize = ticksize

    h1 = hm!(ax1, x, y, Pmap;
        cmap = CMAP,
        vmin = VMIN_PMAX,
        vmax = VMAX_PMAX
    )

    cb1 = Colorbar(
        fig[1,2], h1,
        label = L"P_{\max}\,[\mathrm{K}]",
        width = 18
    )
    cb1.labelsize = labsize
    cb1.ticklabelsize = ticksize

    # --- beta ---
    ax2 = Axis(fig[1,3],
        xlabel = "Distance [pc]",
        ylabel = "Distance [pc]",
        title  = "",
        xticks = (ticks_pc, string.(ticks_pc)),
        yticks = (ticks_pc, string.(ticks_pc)),
        xgridvisible = false,
        ygridvisible = false
    )
    ax2.xlabelsize = labsize
    ax2.ylabelsize = labsize
    ax2.xticklabelsize = ticksize
    ax2.yticklabelsize = ticksize

    h2 = hm!(ax2, x, y, beta_map;
        cmap = CMAP,
        vmin = VMIN_BETA,
        vmax = VMAX_BETA
    )

    cb2 = Colorbar(
        fig[1,4], h2,
        label = L"\beta",
        width = 18
    )
    cb2.labelsize = labsize
    cb2.ticklabelsize = ticksize

    colgap!(fig.layout, 14)

    save(out_pdf, fig)
    display(fig)
end

println("Saved → ", out_pdf)
