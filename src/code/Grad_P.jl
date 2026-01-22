# gradP_cube & gradRM_cube from Qnu/Unu FITS cubes (per-ν slice)
# + column density map N_H (from density.fits) + its gradient map
#
# FIGURE 1 (4 panels):
#   Top row:    |∇P| (left)   |∇RM| (right)
#   Bottom row: N_H [cm^-2] (left)   |∇N_H| [cm^-2 pc^-1] (right)
#   -> 4 colorbars, colormap=:magma
#   -> axes in pc (0..50, ticks every 10 pc)
#   -> big labels/ticks + big colorbar ticklabels
#   -> ν annotated at top (for the chosen k)
#
# FIGURE 2:
#   PDF of normalized |∇P|:
#     f = (|∇P| - <|∇P|>) / σ(|∇P|)
#
# NOTE (Makie compat):
#   Some Makie versions break if colorrange = nothing when a Colorbar is created.
#   We therefore ALWAYS pass a (vmin, vmax) tuple, using data min/max if needed.
# ============================================================

using FITSIO
using Printf
using CairoMakie
using LaTeXStrings
using Statistics
include("Desktop/Depolarization/src/io/fits_io.jl")

# ------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------
const LOS = "y"

const Q_path = LOS == "x" ?
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/x/Synchrotron/WithFaraday/Qnu.fits" :
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Qnu.fits"

const U_path = LOS == "x" ?
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/x/Synchrotron/WithFaraday/Unu.fits" :
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Unu.fits"

const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024"
const DENSITY_PATH = joinpath(SIMU_ROOT, "density.fits")

const Lbox_pc = 50.0
const NPIX = 256
const dx_pc = Lbox_pc / NPIX

const PC_TO_CM = 3.085677581e18
const dz_cm = dx_pc * PC_TO_CM

const P_eps = 1e-12

const ν_min_Hz = 120e6
const ν_max_Hz = 167e6
const dν_Hz = 98e3
const c_ms = 299_792_458.0

const OUTDIR = dirname(Q_path)
const SAVE_CUBES = true
const GRADP_OUT = joinpath(OUTDIR, "gradP_cube.fits")
const GRADRM_OUT = joinpath(OUTDIR, "gradRM_cube.fits")
const NH_OUT = joinpath(OUTDIR, "NH_map.fits")
const GRADNH_OUT = joinpath(OUTDIR, "grad_NH_map.fits")

# ------------------------------------------------------------
# PLOT SETTINGS
# ------------------------------------------------------------
const DO_PLOT_FOURPANELS = true
const k_plot = 10
const SAVE_PLOT = true
const PLOT_EXT = "pdf"
const PLOTPATH = joinpath(
    OUTDIR,
    @sprintf("gradP_gradRM_NH_gradNH_LOS_%s_k%04d.%s", LOS, k_plot, PLOT_EXT)
)

const GRADP_VMIN = 0.0
const GRADP_VMAX = 40.0
const GRADRM_VMIN = 0.0
const GRADRM_VMAX = 1.0

const NH_VMIN = 1e20
const NH_VMAX = 1e21
const GRADNH_VMIN = nothing
const GRADNH_VMAX = nothing

const LABELSIZE = 28
const TICKLABELSIZE = 22
const COLORBAR_LABELSIZE = 26
const COLORBAR_TICKLABELSIZE = 22
const NU_TITLESIZE = 26

const ticks_pc = (0:10:50, string.(0:10:50))

# ------------------------------------------------------------
# FIGURE 2: PDF of normalized |∇P|
# ------------------------------------------------------------
const DO_PDF_GRADP_NORM = true
const NBINS_PDF = 120
const XRANGE_PDF = (-6.0, 6.0)
const PDF_YLOG = false
const PDF_LABELSIZE = 28
const PDF_TICKLABELSIZE = 22
const PDF_LINEWIDTH = 4
const PDFPATH = joinpath(
    OUTDIR,
    @sprintf("PDF_gradPnorm_LOS_%s_k%04d.%s", LOS, k_plot, PLOT_EXT)
)

# ------------------------------------------------------------
# Build ν array (Hz) to match cube third axis exactly
# ------------------------------------------------------------
function build_nu_array(nν::Int; ν_min_Hz::Float64, dν_Hz::Float64)
    ν = ν_min_Hz .+ (0:nν-1) .* dν_Hz
    return ν
end

# ------------------------------------------------------------
# Finite differences (2D) on A[y,x]
# ------------------------------------------------------------
function d_dx(A::AbstractMatrix)
    ny, nx = size(A)
    out = similar(A, Float64)
    @inbounds begin
        for i in 1:ny
            for j in 2:nx-1
                out[i, j] = 0.5 * (A[i, j+1] - A[i, j-1])
            end
            out[i, 1]  = A[i, 2]  - A[i, 1]
            out[i, nx] = A[i, nx] - A[i, nx-1]
        end
    end
    out
end

function d_dy(A::AbstractMatrix)
    ny, nx = size(A)
    out = similar(A, Float64)
    @inbounds begin
        for j in 1:nx
            for i in 2:ny-1
                out[i, j] = 0.5 * (A[i+1, j] - A[i-1, j])
            end
            out[1, j]  = A[2, j]  - A[1, j]
            out[ny, j] = A[ny, j] - A[ny-1, j]
        end
    end
    out
end

# ------------------------------------------------------------
# Makie-safe colorrange builder
# ------------------------------------------------------------
function safe_crange(A; vmin=nothing, vmax=nothing)
    Amin = minimum(A)
    Amax = maximum(A)
    lo = vmin === nothing ? Amin : vmin
    hi = vmax === nothing ? Amax : vmax
    lo_f = Float32(lo)
    hi_f = Float32(hi)
    if !(isfinite(lo_f) && isfinite(hi_f)) || lo_f == hi_f
        eps = Float32(1f-6)
        return (lo_f - eps, hi_f + eps)
    end
    return (lo_f, hi_f)
end

# ------------------------------------------------------------
# Compute cubes per ν:
# |∇P| = sqrt( (dQ/dx)^2+(dQ/dy)^2+(dU/dx)^2+(dU/dy)^2 )
# |∇RM| = |∇P| / (2 λ(ν)^2 |P|)
# ------------------------------------------------------------
function compute_grad_cubes(Q::Array{<:Real,3}, U::Array{<:Real,3};
                            dx_pc::Float64, P_eps::Float64,
                            ν_min_Hz::Float64, dν_Hz::Float64, c_ms::Float64, ν_max_Hz::Float64)

    ny, nx, nν = size(Q)
    gradP  = Array{Float32}(undef, ny, nx, nν)
    gradRM = Array{Float32}(undef, ny, nx, nν)

    ν  = build_nu_array(nν; ν_min_Hz=ν_min_Hz, dν_Hz=dν_Hz)
    λ2 = (c_ms ./ ν) .^ 2

    @printf("[Info] Cube size: (%d, %d, %d)\n", ny, nx, nν)
    @printf("[Info] ν[1]=%.3f MHz, ν[end]=%.3f MHz, Δν=%.3f kHz\n",
            ν[1]/1e6, ν[end]/1e6, dν_Hz/1e3)
    @printf("[Info] Target ν_max (user) = %.3f MHz\n", ν_max_Hz/1e6)

    for k in 1:nν
        Qk = Float64.(Q[:, :, k])
        Uk = Float64.(U[:, :, k])

        dQdx = d_dx(Qk) ./ dx_pc
        dQdy = d_dy(Qk) ./ dx_pc
        dUdx = d_dx(Uk) ./ dx_pc
        dUdy = d_dy(Uk) ./ dx_pc

        gP = sqrt.(dQdx.^2 .+ dQdy.^2 .+ dUdx.^2 .+ dUdy.^2)
        P  = sqrt.(Qk.^2 .+ Uk.^2)

        gRM = gP ./ (2.0 * λ2[k] .* max.(P, P_eps))

        gradP[:, :, k]  .= Float32.(gP)
        gradRM[:, :, k] .= Float32.(gRM)
    end

    return gradP, gradRM, ν
end

# ------------------------------------------------------------
# N_H map (integrate density cube along LOS) + gradient magnitude
# density cube assumed axes ~ (x,y,z). We return a 2D map in [y,x] convention.
# N_H [cm^-2] = ∑ n_H [cm^-3] * dz [cm]
# |∇N_H| [cm^-2 pc^-1] from spatial derivatives in pc.
# ------------------------------------------------------------
function compute_NH_and_grad(density_cube::AbstractArray{<:Real,3};
                             LOS::String, dx_pc::Float64, dz_cm::Float64)

    axis = LOS == "x" ? 1 : (LOS == "y" ? 2 : 3)

    # Column (sum along LOS) and convert to cm^-2
    NH = dropdims(sum(Float64.(density_cube), dims=axis), dims=axis) .* dz_cm

    # Convert to [y,x] convention
    if axis == 1
        NH_map = NH                     # (y,z) -> treat z as x
    elseif axis == 2
        NH_map = permutedims(NH, (2,1)) # (x,z) -> (z,x)
    else
        NH_map = permutedims(NH, (2,1)) # (x,y) -> (y,x)
    end

    dNdx = d_dx(NH_map) ./ dx_pc
    dNdy = d_dy(NH_map) ./ dx_pc
    grad_NH = sqrt.(dNdx.^2 .+ dNdy.^2)

    return Float32.(NH_map), Float32.(grad_NH)
end

# ------------------------------------------------------------
# Plot 4 panels (2 rows x 2 cols) + ν annotation on top
# ------------------------------------------------------------
function plot_four_panels(gradP_cube::Array{<:Real,3}, gradRM_cube::Array{<:Real,3},
                          NH::AbstractMatrix{<:Real}, grad_NH::AbstractMatrix{<:Real};
                          k::Int, Lbox_pc::Float64, ν_Hz::AbstractVector{<:Real},
                          savepath::String)

    ny, nx, nν = size(gradP_cube)
    @assert size(gradRM_cube) == (ny, nx, nν)
    @assert size(NH) == (ny, nx)
    @assert size(grad_NH) == (ny, nx)
    @assert 1 <= k <= nν "k_plot=$k hors bornes (1..$nν)"

    x_pc = range(0, Lbox_pc; length=nx)
    y_pc = range(0, Lbox_pc; length=ny)

    crP   = safe_crange(gradP_cube[:, :, k]; vmin=GRADP_VMIN,   vmax=GRADP_VMAX)
    crRM  = safe_crange(gradRM_cube[:, :, k]; vmin=GRADRM_VMIN, vmax=GRADRM_VMAX)
    crNH  = safe_crange(NH;       vmin=NH_VMIN,       vmax=NH_VMAX)
    crgNH = safe_crange(grad_NH;  vmin=GRADNH_VMIN,   vmax=GRADNH_VMAX)

    ν_MHz = ν_Hz[k] / 1e6
    nu_label = latexstring(@sprintf("\\nu = %.2f~\\mathrm{MHz}", ν_MHz))

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1200, 1000))
        fig[0, 1:4] = Label(fig, nu_label; fontsize=NU_TITLESIZE, tellwidth=false)

        # --- Top row ---
        ax11 = Axis(fig[1, 1],
            xlabel="Distance [pc]", ylabel="Distance [pc]",
            xticks=ticks_pc, yticks=ticks_pc,
            limits=(0,50,0,50),
            xlabelsize=LABELSIZE, ylabelsize=LABELSIZE,
            xticklabelsize=TICKLABELSIZE, yticklabelsize=TICKLABELSIZE
        )
        ax13 = Axis(fig[1, 3],
            xlabel="Distance [pc]", ylabel="",
            xticks=ticks_pc,
            limits=(0,50,0,50),
            xlabelsize=LABELSIZE, xticklabelsize=TICKLABELSIZE,
            yticksvisible=false, yticklabelsvisible=false
        )

        hm11 = heatmap!(ax11, x_pc, y_pc, gradP_cube[:, :, k];  colormap=:magma, colorrange=crP)
        hm13 = heatmap!(ax13, x_pc, y_pc, gradRM_cube[:, :, k]; colormap=:magma, colorrange=crRM)

        cb12 = Colorbar(fig[1, 2], hm11; labelsize=COLORBAR_LABELSIZE, ticklabelsize=COLORBAR_TICKLABELSIZE)
        cb14 = Colorbar(fig[1, 4], hm13; labelsize=COLORBAR_LABELSIZE, ticklabelsize=COLORBAR_TICKLABELSIZE)
        cb12.label = L"|\nabla \mathbf{P}|"
        cb14.label = L"|\nabla \mathrm{RM}|"

        # --- Bottom row ---
        ax21 = Axis(fig[2, 1],
            xlabel="Distance [pc]", ylabel="Distance [pc]",
            xticks=ticks_pc, yticks=ticks_pc,
            limits=(0,50,0,50),
            xlabelsize=LABELSIZE, ylabelsize=LABELSIZE,
            xticklabelsize=TICKLABELSIZE, yticklabelsize=TICKLABELSIZE
        )
        ax23 = Axis(fig[2, 3],
            xlabel="Distance [pc]", ylabel="",
            xticks=ticks_pc,
            limits=(0,50,0,50),
            xlabelsize=LABELSIZE, xticklabelsize=TICKLABELSIZE,
            yticksvisible=false, yticklabelsvisible=false
        )

        hm21 = heatmap!(ax21, x_pc, y_pc, NH;      colormap=:magma, colorrange=crNH)
        hm23 = heatmap!(ax23, x_pc, y_pc, grad_NH; colormap=:magma, colorrange=crgNH)

        cb22 = Colorbar(fig[2, 2], hm21; labelsize=COLORBAR_LABELSIZE, ticklabelsize=COLORBAR_TICKLABELSIZE)
        cb24 = Colorbar(fig[2, 4], hm23; labelsize=COLORBAR_LABELSIZE, ticklabelsize=COLORBAR_TICKLABELSIZE)
        cb22.label = L"N_H\ [\mathrm{cm}^{-2}]"
        cb24.label = L"|\nabla N_H|\ [\mathrm{cm}^{-2}\,\mathrm{pc}^{-1}]"

        display(fig)
        if SAVE_PLOT
            save(savepath, fig)
            @printf("[OK] Saved 4-panel plot: %s\n", savepath)
        end
    end
end

# ------------------------------------------------------------
# FIGURE 2: PDF of normalized |∇P| (for one channel k)
# ------------------------------------------------------------
function gradP_normalized_slice(gradP_cube::AbstractArray{<:Real,3}, k::Int)
    A = Float64.(gradP_cube[:, :, k])
    μ = mean(A)
    σ = std(A)
    @assert σ > 0 "std(|∇P|)=0 (image is constant?)"
    return (A .- μ) ./ σ, μ, σ
end

function plot_pdf_gradP_norm(gradP_cube::AbstractArray{<:Real,3};
                             k::Int, ν_Hz::AbstractVector{<:Real},
                             savepath::String)

    fA, μ, σ = gradP_normalized_slice(gradP_cube, k)
    vecf = vec(fA)

    lo, hi = XRANGE_PDF
    mask = (vecf .>= lo) .& (vecf .<= hi)
    v = vecf[mask]

    edges = collect(range(lo, hi; length=NBINS_PDF + 1))
    counts = zeros(Float64, NBINS_PDF)
    # manual binning (no extra deps)
    binw = edges[2] - edges[1]
    @inbounds for x in v
        b = Int(floor((x - lo) / binw)) + 1
        if 1 <= b <= NBINS_PDF
            counts[b] += 1
        end
    end
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    pdf = counts ./ (sum(counts) * binw)

    ν_MHz = ν_Hz[k] / 1e6
    title_lbl = latexstring(@sprintf("\\nu = %.2f~\\mathrm{MHz}", ν_MHz))

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(900, 650))
        fig[0, 1] = Label(fig, title_lbl; fontsize=NU_TITLESIZE, tellwidth=false)

        ax = Axis(fig[1, 1],
            xlabel = L"f = (|\nabla \mathbf{P}| - \langle|\nabla \mathbf{P}|\rangle)/\sigma(|\nabla \mathbf{P}|)",
            ylabel = "PDF",
            xlabelsize = PDF_LABELSIZE,
            ylabelsize = PDF_LABELSIZE,
            xticklabelsize = PDF_TICKLABELSIZE,
            yticklabelsize = PDF_TICKLABELSIZE,
        )
        if PDF_YLOG
            ax.yscale = log10
        end

        lines!(ax, centers, pdf; linewidth=PDF_LINEWIDTH)

        display(fig)
        save(savepath, fig)
        @printf("[OK] Saved PDF(|∇P| normalized): %s\n", savepath)
        @printf("[Info] mean(|∇P|)=%.6g, std(|∇P|)=%.6g\n", μ, σ)
    end
end

# ------------------------------------------------------------
# RUN
# ------------------------------------------------------------
@printf("[Info] Reading Q: %s\n", Q_path)
@printf("[Info] Reading U: %s\n", U_path)

Q = read_FITS(Q_path)
U = read_FITS(U_path)

@assert ndims(Q) == 3 "Q doit être un cube (ny,nx,nν). ndims(Q)=$(ndims(Q))"
@assert size(Q) == size(U) "Q et U n'ont pas la même taille"

gradP_cube, gradRM_cube, ν_Hz = compute_grad_cubes(
    Q, U;
    dx_pc=dx_pc,
    P_eps=P_eps,
    ν_min_Hz=ν_min_Hz,
    dν_Hz=dν_Hz,
    c_ms=c_ms,
    ν_max_Hz=ν_max_Hz
)

if SAVE_CUBES
    write_FITS(GRADP_OUT,  gradP_cube)
    write_FITS(GRADRM_OUT, gradRM_cube)
    @printf("[OK] Saved cubes:\n  %s\n  %s\n", GRADP_OUT, GRADRM_OUT)
end

@printf("[Info] Reading density cube: %s\n", DENSITY_PATH)
density_cube = read_FITS(DENSITY_PATH)
@assert ndims(density_cube) == 3 "density.fits doit être un cube 3D. ndims=$(ndims(density_cube))"

NH_map, grad_NH_map = compute_NH_and_grad(density_cube; LOS=LOS, dx_pc=dx_pc, dz_cm=dz_cm)

if SAVE_CUBES
    write_FITS(NH_OUT, NH_map)
    write_FITS(GRADNH_OUT, grad_NH_map)
    @printf("[OK] Saved N_H maps:\n  %s\n  %s\n", NH_OUT, GRADNH_OUT)
end

if DO_PLOT_FOURPANELS
    plot_four_panels(gradP_cube, gradRM_cube, NH_map, grad_NH_map;
        k=k_plot, Lbox_pc=Lbox_pc, ν_Hz=ν_Hz, savepath=PLOTPATH
    )
end

if DO_PDF_GRADP_NORM
    plot_pdf_gradP_norm(gradP_cube; k=k_plot, ν_Hz=ν_Hz, savepath=PDFPATH)
end

println("\nDone.")
