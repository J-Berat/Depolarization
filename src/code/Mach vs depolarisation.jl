############################################################
# Mach vs depolarisation – 4-panel layout
#
# Top row    : LOS = y,   P10 in [1, 6]   (NO xlabels / xticklabels)
# Bottom row : LOS = x,   P10 in [0, 0.5]
#
# Columns:
#   left  -> <Ms>_3D
#   right -> <MA>_3D
#
# Colours by rms
# No grids anywhere
############################################################

using Statistics
using CairoMakie
using LaTeXStrings

include("fits_io.jl")

# ----------------------------------------------------------
# Parameters
# ----------------------------------------------------------
ROOT  = "/Users/jb270005/Desktop/simu_RAMSES"
LOS_Y = "y"
LOS_X = "x"

USE_FINITE_ONLY = true

mu = 1.27
γ  = 5/3
kB = 1.380649e-16
mH = 1.6735575e-24

V_UNIT = :kms
B_MG_TO_UG = 1000.0

# ----------------------------------------------------------
# Simulations
# ----------------------------------------------------------
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

# ----------------------------------------------------------
# Helpers
# ----------------------------------------------------------
phys_dir(s) = joinpath(ROOT, s)
pmax_path(s; los) =
    joinpath(ROOT, s, los, "Synchrotron", "WithFaraday", "Pmax.fits")

to_cgs(V) = V_UNIT == :kms ? (V .* 1e5) : V
cs(T) = sqrt.(γ * kB .* T ./ (mu * mH))
rho(n) = (mu * mH) .* n

depol_P10(P) = quantile(USE_FINITE_ONLY ? P[isfinite.(P)] : vec(P), 0.10)

rms_label(s) = match(r"rms(.*?)nograv1024", s).captures[1]

# ----------------------------------------------------------
# Colours
# ----------------------------------------------------------
rms_colors = Dict(
    "09000" => :royalblue,
    "18000" => :darkorange,
    "36000" => :seagreen,
    "72000" => :firebrick
)

# ----------------------------------------------------------
# Compute
# ----------------------------------------------------------
Mach3D, MachA3D = Float64[], Float64[]
P10_y, P10_x    = Float64[], Float64[]
RMSlbl          = String[]

for simu in SimuList
    dir = phys_dir(simu)

    n  = Float32.(read_FITS(joinpath(dir, "density.fits")))
    T  = Float32.(read_FITS(joinpath(dir, "temperature.fits")))

    Bx = Float32.(read_FITS(joinpath(dir, "Bx.fits"))) .* Float32(B_MG_TO_UG) .* 1f-6
    By = Float32.(read_FITS(joinpath(dir, "By.fits"))) .* Float32(B_MG_TO_UG) .* 1f-6
    Bz = Float32.(read_FITS(joinpath(dir, "Bz.fits"))) .* Float32(B_MG_TO_UG) .* 1f-6

    Vx = to_cgs(Float32.(read_FITS(joinpath(dir, "Vx.fits"))))
    Vy = to_cgs(Float32.(read_FITS(joinpath(dir, "Vy.fits"))))
    Vz = to_cgs(Float32.(read_FITS(joinpath(dir, "Vz.fits"))))

    vmag = sqrt.(Vx.^2 .+ Vy.^2 .+ Vz.^2)
    Ms   = mean(vmag ./ cs(T))

    Bmag = sqrt.(Bx.^2 .+ By.^2 .+ Bz.^2)
    vA   = Bmag ./ sqrt.(4f0 * Float32(pi) .* rho(n))
    MA   = mean(vmag ./ vA)

    push!(Mach3D, Ms)
    push!(MachA3D, MA)
    push!(P10_y, depol_P10(Float32.(read_FITS(pmax_path(simu; los=LOS_Y)))))
    push!(P10_x, depol_P10(Float32.(read_FITS(pmax_path(simu; los=LOS_X)))))
    push!(RMSlbl, rms_label(simu))
end

# ----------------------------------------------------------
# Plot (2x2 layout)
# ----------------------------------------------------------
with_theme(theme_latexfonts()) do
    fig = Figure(size = (1100, 900))
    gl  = GridLayout(fig[1, 1], colgap = 0, rowgap = 0)

    # --- Top row (LOS y) : NO xlabels / xticklabels
    ax_y_Ms = Axis(
        gl[1, 1],
        ylabel = L"$P_{10}(P_{\max})$",
        xlabel = "",
        xticksvisible = false,
        xticklabelsvisible = false,
        xlabelsize = 34,
        ylabelsize = 34,
        yticklabelsize = 26,
        xgridvisible = false,
        ygridvisible = false
    )

    ax_y_MA = Axis(
        gl[1, 2],
        xlabel = "",
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false
    )

    # --- Bottom row (LOS x) : xlabels ON
    ax_x_Ms = Axis(
        gl[2, 1],
        ylabel = L"$P_{10}(P_{\max})$",
        xlabel = L"$\langle \mathcal{M}_s \rangle_{3D}$",
        xlabelsize = 34,
        ylabelsize = 34,
        xticklabelsize = 26,
        yticklabelsize = 26,
        xgridvisible = false,
        ygridvisible = false
    )

    ax_x_MA = Axis(
        gl[2, 2],
        xlabel = L"$\langle \mathcal{M}_A \rangle_{3D}$",
        xlabelsize = 34,
        xticklabelsize = 26,
        yticksvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false
    )

    # --- Link Y by row
    linkyaxes!(ax_y_Ms, ax_y_MA)
    linkyaxes!(ax_x_Ms, ax_x_MA)

    ylims!(ax_y_Ms, 1, 6)
    ylims!(ax_x_Ms, 0, 0.5)

    # --- X limits
    xlims!(ax_y_Ms, 0, 3.5)
    xlims!(ax_x_Ms, 0, 3.5)
    xlims!(ax_y_MA, 0, 1.0)
    xlims!(ax_x_MA, 0, 1.0)

    ax_y_MA.xticks = 0:0.1:1.0
    ax_x_MA.xticks = 0:0.1:1.0

    # --- Scatter
    for (rms, col) in rms_colors
        idx = findall(RMSlbl .== rms)

        scatter!(ax_y_Ms, Mach3D[idx],  P10_y[idx], color=col, markersize=14)
        scatter!(ax_y_MA, MachA3D[idx], P10_y[idx], color=col, markersize=14, label=rms)

        scatter!(ax_x_Ms, Mach3D[idx],  P10_x[idx], color=col, markersize=14)
        scatter!(ax_x_MA, MachA3D[idx], P10_x[idx], color=col, markersize=14)
    end

    axislegend(ax_y_MA; title="rms", position=:rt, labelsize=24, titlesize=26)

    display(fig)
    save("mach_vs_depol_4panels_LOSxy.pdf", fig)
end
