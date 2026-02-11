using FITSIO
using CairoMakie
using Random
using Statistics
using LaTeXStrings

# ============================================================
# Simulation list
# ============================================================
const SimuList = [
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

# ============================================================
# FITS I/O
# ============================================================
readfits(path::AbstractString) = FITS(path, "r") do f
    read(f[1])
end

base_dir(s) = joinpath("/Users/jb270005/Desktop/simu_RAMSES", s)
Bx_path(s)  = joinpath(base_dir(s), "Bx.fits")
By_path(s)  = joinpath(base_dir(s), "By.fits")
Bz_path(s)  = joinpath(base_dir(s), "Bz.fits")
ne_path(s, los) = joinpath(base_dir(s), lowercase(los), "Synchrotron", "ne.fits")

# ============================================================
# Simulation selector
# ============================================================
simu(x) = x isa Integer ? SimuList[x] :
          x isa AbstractString && all(isdigit, x) ? SimuList[parse(Int, x)] :
          x

# ============================================================
# ψ cube (LOS → dim 3)
# ============================================================
function psi_los(Bx, By, Bz; los="z")
    los = lowercase(los)
    los == "z" && return atan.(By, Bx)
    los == "x" && return permutedims(atan.(Bz, By), (2, 3, 1))
    los == "y" && return permutedims(atan.(Bz, Bx), (1, 3, 2))
    error("LOS must be x, y or z")
end

# ============================================================
# B_LOS (µG) reordered to (POS1, POS2, LOS)
# ============================================================
function Blos_los(Bx, By, Bz; los="z")
    los = lowercase(los)
    los == "z" && return 1000 .* Bz
    los == "x" && return 1000 .* permutedims(Bx, (2, 3, 1))
    los == "y" && return 1000 .* permutedims(By, (1, 3, 2))
    error("LOS must be x, y or z")
end

# ============================================================
# ne reordered to (POS1, POS2, LOS)
# ============================================================
ne_los(ne; los="z") =
    los == "z" ? ne :
    los == "x" ? permutedims(ne, (2, 3, 1)) :
    los == "y" ? permutedims(ne, (1, 3, 2)) :
    error("LOS must be x, y or z")

# ============================================================
# Circular mean along LOS (deg)
# ============================================================
mean_psi_deg(ψ) =
    dropdims(
        rad2deg.(atan.(mean(sin.(ψ), dims=3),
                        mean(cos.(ψ), dims=3))),
        dims=3
    )

# ============================================================
# Plot
# ============================================================
function plot_meanpsi_sightline(
        s;
        los="z",
        seed=123,
        Lbox_pc=50.0,
        colormap=:twilight
    )

    Random.seed!(seed)

    Bx = readfits(Bx_path(s))
    By = readfits(By_path(s))
    Bz = readfits(Bz_path(s))
    ne = readfits(ne_path(s, los))

    ψ    = psi_los(Bx, By, Bz; los=los)
    ψ̄   = mean_psi_deg(ψ)
    BLOS = Blos_los(Bx, By, Bz; los=los)
    ne3  = ne_los(ne; los=los)

    n1, n2, nlos = size(ψ)
    i, j = rand(1:n1), rand(1:n2)

    dlos    = range(0, Lbox_pc; length=nlos)
    ψ_sl    = rad2deg.(ψ[i, j, :])
    BLOS_sl = BLOS[i, j, :]
    ne_sl   = ne3[i, j, :]

    xmap = range(0, Lbox_pc; length=n1)
    ymap = range(0, Lbox_pc; length=n2)
    ticks = 0:10:50  # force 0..50

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1250, 950))

        labelsize = 28
        ticksize  = 24
        titlesize_map = 36
        ticksize_cbar = 26
        legendsize = 30

        # =============================
        # Mean orientation map
        # =============================
        gl = GridLayout()
        fig[1, 1] = gl
        colgap!(gl, 8)

        axm = Axis(
            gl[1, 1],
            xlabel="Distance [pc]",
            ylabel="Distance [pc]",
            title=L"\vec{B}_{\mathrm{LOS}}\ \text{mean orientation}",
            titlesize=titlesize_map,
            titlefont=:regular,
            xlabelsize=labelsize,
            ylabelsize=labelsize,
            xticklabelsize=ticksize,
            yticklabelsize=ticksize
        )

        hm = heatmap!(axm, xmap, ymap, ψ̄'; colormap=colormap)
        Colorbar(gl[1, 2], hm; labelvisible=false, ticklabelsize=ticksize_cbar)
        colsize!(gl, 2, Fixed(30))

        vlines!(axm, [xmap[i]], linestyle=:dash, color=(:white, 0.8))
        hlines!(axm, [ymap[j]], linestyle=:dash, color=(:white, 0.8))

        axm.xticks = (ticks, string.(Int.(ticks)))
        axm.yticks = (ticks, string.(Int.(ticks)))

        # =============================
        # ψ(s)
        # =============================
        axsψ = Axis(
            fig[1, 2],
            xlabel="Distance [pc]",
            ylabel=L"\psi\ [^\circ]",
            xlabelsize=labelsize,
            ylabelsize=labelsize,
            xticklabelsize=ticksize,
            yticklabelsize=ticksize,
            xgridvisible=false,
            ygridvisible=false
        )

        lines!(axsψ, dlos, ψ_sl; linewidth=2, color=(:black, 0.75))
        hlines!(axsψ, [0.0], linestyle=:dash, color=(:black, 0.35))
        axsψ.xticks = (ticks, string.(Int.(ticks)))
        xlims!(axsψ, 0, 50)

        # =============================
        # BLOS + ne (shared x-axis)
        # =============================
        axB = Axis(
            fig[2, 1:2],
            xlabel="Distance [pc]",
            ylabel=L"B_{\mathrm{LOS}}\ [\mu\mathrm{G}]",
            xlabelsize=labelsize,
            ylabelsize=labelsize,
            xticklabelsize=ticksize,
            yticklabelsize=ticksize,
            xgridvisible=false,
            ygridvisible=false
        )

        axne = Axis(
            fig[2, 1:2],
            yaxisposition=:right,
            ylabel=L"n_e\ [\mathrm{cm}^{-3}]",
            ylabelsize=labelsize,
            yticklabelsize=ticksize,
            xticksvisible=false,
            xticklabelsvisible=false,  # ✅ correct attribute
            xgridvisible=false,
            ygridvisible=false
        )

        linkxaxes!(axB, axne)

        lB  = lines!(axB,  dlos, BLOS_sl; linewidth=2, color=(:black, 0.70))
        lne = lines!(axne, dlos, ne_sl;   linewidth=2, color=(:red, 0.55), linestyle=:dash)

        axB.xticks = (ticks, string.(Int.(ticks)))
        xlims!(axB, 0, 50)
        xlims!(axne, 0, 50)  # safety

        Legend(
            fig[2, 1:2],
            [lB, lne],
            [L"B_{\mathrm{LOS}}", L"n_e"],
            labelsize=legendsize,
            framevisible=false,
            halign=:left,
            valign=:top
        )

        # ---- row sizing after rows exist
        rowsize!(fig.layout, 1, Relative(0.42))
        rowsize!(fig.layout, 2, Relative(0.58))

        rowgap!(fig.layout, 18)
        colgap!(fig.layout, 18)

        display(fig)
        return fig, (i=i, j=j)
    end
end

# ============================================================
# RUN
# ============================================================
S   = simu(7)
LOS = "x"

fig, ij = plot_meanpsi_sightline(S; los=LOS, seed=123)
