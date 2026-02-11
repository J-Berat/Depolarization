############################################################
# Mach vs depolarisation – 4-panel layout
#
# Y is now the FRACTION of pixels below a FIXED threshold Pth
# Pth is defined as the GLOBAL 1st decile of Pmax across all sims
# (computed separately for each LOS, so "fixed per LOS").
#
# Top row    : LOS = y,   frac in [0, 1]
# Bottom row : LOS = x,   frac in [0, 1]
#
# Columns:
#   left  -> <Ms>_3D
#   right -> <MA>_3D
#
# Colours by rms
# No grids anywhere
# Remove xlabels/xticklabels from first row (as requested previously)
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

function rms_label(s::String)
    m = match(r"rms(.*?)nograv1024", s)
    return m === nothing ? s : m.captures[1]
end

function finite_vec(A)
    v = vec(A)
    return USE_FINITE_ONLY ? v[isfinite.(v)] : v
end

function frac_below_threshold(A, thr)
    v = finite_vec(A)
    return count(<=(thr), v) / length(v)
end

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
# 1) Compute global fixed thresholds Pth (one per LOS)
#    Pth_y = global 10th percentile of Pmax over ALL simulations (LOS=y)
#    Pth_x = global 10th percentile of Pmax over ALL simulations (LOS=x)
# ----------------------------------------------------------
function global_P10_threshold(SimuList, los)
    allvals = Float32[]
    for simu in SimuList
        P = Float32.(read_FITS(pmax_path(simu; los=los)))
        append!(allvals, finite_vec(P))
    end
    return quantile(allvals, 0.10)
end

Pth_y = global_P10_threshold(SimuList, LOS_Y)
Pth_x = global_P10_threshold(SimuList, LOS_X)

@info "Fixed threshold Pth (LOS=y) = $Pth_y"
@info "Fixed threshold Pth (LOS=x) = $Pth_x"

# ----------------------------------------------------------
# 2) Compute Mach numbers and depolarised fractions per sim
# ----------------------------------------------------------
Mach3D, MachA3D = Float64[], Float64[]
fdep_y, fdep_x  = Float64[], Float64[]
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

    P_y = Float32.(read_FITS(pmax_path(simu; los=LOS_Y)))
    P_x = Float32.(read_FITS(pmax_path(simu; los=LOS_X)))

    push!(Mach3D, Ms)
    push!(MachA3D, MA)
    push!(fdep_y, frac_below_threshold(P_y, Pth_y))
    push!(fdep_x, frac_below_threshold(P_x, Pth_x))
    push!(RMSlbl, rms_label(simu))
end

# ----------------------------------------------------------
# 3) Plot (2x2 layout): Y is a fraction in [0,1]
# ----------------------------------------------------------
with_theme(theme_latexfonts()) do
    fig = Figure(size = (1100, 900))
    gl  = GridLayout(fig[1, 1], colgap = 0, rowgap = 0)

    # --- Top row (LOS y): NO xlabels/xticklabels
    ax_y_Ms = Axis(
        gl[1, 1],
        ylabel = L"$f\,(P_{\max}\le P_{\mathrm{th}})$",
        xlabel = "",
        xticksvisible = false,
        xticklabelsvisible = false,
        xlabelsize = 34, ylabelsize = 34,
        yticklabelsize = 26,
        xgridvisible = false, ygridvisible = false
    )

    ax_y_MA = Axis(
        gl[1, 2],
        xlabel = "",
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false
    )

    # --- Bottom row (LOS x): xlabels ON
    ax_x_Ms = Axis(
        gl[2, 1],
        ylabel = L"$f\,(P_{\max}\le P_{\mathrm{th}})$",
        xlabel = L"$\langle \mathcal{M}_s \rangle_{3D}$",
        xlabelsize = 34, ylabelsize = 34,
        xticklabelsize = 26, yticklabelsize = 26,
        xgridvisible = false, ygridvisible = false
    )

    ax_x_MA = Axis(
        gl[2, 2],
        xlabel = L"$\langle \mathcal{M}_A \rangle_{3D}$",
        xlabelsize = 34, xticklabelsize = 26,
        yticksvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false
    )

    # Link Y by row + set [0,1]
    linkyaxes!(ax_y_Ms, ax_y_MA)
    linkyaxes!(ax_x_Ms, ax_x_MA)
    ylims!(ax_y_Ms, 0, 0.5)
    ylims!(ax_x_Ms, 0, 0.5)

    # X limits
    xlims!(ax_y_Ms, 0, 3.5)
    xlims!(ax_x_Ms, 0, 3.5)
    xlims!(ax_y_MA, 0, 1.0)
    xlims!(ax_x_MA, 0, 1.0)
    ax_y_MA.xticks = 0:0.1:1.0
    ax_x_MA.xticks = 0:0.1:1.0

    # Scatter
    for (rms, col) in rms_colors
        idx = findall(RMSlbl .== rms)

        scatter!(ax_y_Ms, Mach3D[idx],  fdep_y[idx], color=col, markersize=14)
        scatter!(ax_y_MA, MachA3D[idx], fdep_y[idx], color=col, markersize=14, label=rms)

        scatter!(ax_x_Ms, Mach3D[idx],  fdep_x[idx], color=col, markersize=14)
        scatter!(ax_x_MA, MachA3D[idx], fdep_x[idx], color=col, markersize=14)
    end

    axislegend(ax_y_MA; title="rms", position=:rt, labelsize=24, titlesize=26)

    display(fig)
    save("mach_vs_depol_fraction_fixed_threshold.pdf", fig)
end
