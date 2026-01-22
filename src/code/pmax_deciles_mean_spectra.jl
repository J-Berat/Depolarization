using FITSIO
using CairoMakie
using Statistics
using Printf
using LaTeXStrings

include(joinpath(@__DIR__, "../constants.jl"))
include(joinpath(@__DIR__, "../io/fits_io.jl"))  # provides read_FITS, write_FITS

# ------------------------------------------------------------
# USER CHOICES
# ------------------------------------------------------------
const SIMU_ROOT = DepolarizationConstants.PmaxDecilesMeanSpectra.SIMU_ROOT
const SIMU_NAME = DepolarizationConstants.PmaxDecilesMeanSpectra.SIMU_NAME
const LOS = DepolarizationConstants.PmaxDecilesMeanSpectra.LOS

const LBOX_PC = DepolarizationConstants.PmaxDecilesMeanSpectra.LBOX_PC
const SAVE_FIG = DepolarizationConstants.PmaxDecilesMeanSpectra.SAVE_FIG
const OUT_DIR = DepolarizationConstants.PmaxDecilesMeanSpectra.OUT_DIR
const OUT_EXT = DepolarizationConstants.PmaxDecilesMeanSpectra.OUT_EXT

const USE_COLRANGE = DepolarizationConstants.PmaxDecilesMeanSpectra.USE_COLRANGE
const VMIN = DepolarizationConstants.PmaxDecilesMeanSpectra.VMIN
const VMAX = DepolarizationConstants.PmaxDecilesMeanSpectra.VMAX

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
pmax_path = joinpath(SIMU_ROOT, SIMU_NAME, LOS, "Synchrotron", "WithFaraday", "Pmax.fits")
fdf_path  = joinpath(SIMU_ROOT, SIMU_NAME, LOS, "Synchrotron", "WithFaraday", "FDF.fits")

isfile(pmax_path) || error("File not found:\n  $pmax_path")
isfile(fdf_path)  || error("File not found:\n  $fdf_path")

# ------------------------------------------------------------
# WCS helper (for phi axis)
# ------------------------------------------------------------
function fits_axis(header::Dict{String,Any}, n::Int, ax::Int)
    crval = get(header, "CRVAL$ax", 0.0)
    crpix = get(header, "CRPIX$ax", 1.0)
    cdelt = get(header, "CDELT$ax", nothing)
    cdelt === nothing && (cdelt = get(header, "CD$(ax)_$(ax)", 1.0))
    cdelt = float(cdelt)
    return [crval + (i - crpix) * cdelt for i in 1:n]
end

# ------------------------------------------------------------
# READ DATA (using repo helpers)
# ------------------------------------------------------------
Pmax = Float32.(Array(read_FITS(pmax_path)))
rawF = Float32.(Array(read_FITS(fdf_path)))

ndims(rawF) == 3 || error("FDF must be 3D, got ndims=$(ndims(rawF)) size=$(size(rawF))")

# Read header for FDF to build phi (we keep this minimal; only header access)
hdF = FITS(fdf_path, "r") do f
    h = read_header(f[1])
    hd = Dict{String,Any}()
    for k in keys(h)
        try
            hd[string(k)] = h[k]
        catch
        end
    end
    hd
end

# Reorder FDF to (nphi, nx, ny)
if haskey(hdF, "CRVAL3") || haskey(hdF, "CDELT3") || haskey(hdF, "CD3_3")
    nphi = size(rawF, 3)
    phi  = fits_axis(hdF, nphi, 3)
    FDF  = permutedims(rawF, (3, 1, 2))   # (nphi, nx, ny)
elseif haskey(hdF, "CRVAL1") || haskey(hdF, "CDELT1") || haskey(hdF, "CD1_1")
    nphi = size(rawF, 1)
    phi  = fits_axis(hdF, nphi, 1)
    FDF  = rawF                            # assume already (nphi, nx, ny)
else
    nphi = size(rawF, 1)
    phi  = collect(1:nphi)
    FDF  = rawF
end

nx, ny = size(Pmax)
(size(FDF, 2) == nx && size(FDF, 3) == ny) || error(
    "FDF spatial dims must match Pmax.\n" *
    "  size(Pmax) = ($(nx), $(ny))\n" *
    "  size(FDF)  = $(size(FDF)) (expected (nphi, $nx, $ny))"
)

# ------------------------------------------------------------
# DECILES
# ------------------------------------------------------------
p = vec(Float64.(Pmax))
p = p[.!isnan.(p)]

q10 = quantile(p, 0.10)
q90 = quantile(p, 0.90)

mask1  = (Pmax .<= q10) .& .!isnan.(Pmax)
mask10 = (Pmax .>= q90) .& .!isnan.(Pmax)

# ------------------------------------------------------------
# MEAN SPECTRA
# ------------------------------------------------------------
function mean_spectrum(FDF::AbstractArray{<:Real,3}, mask::AbstractArray{Bool,2})
    nphi = size(FDF, 1)
    out  = Vector{Float64}(undef, nphi)
    @inbounds for k in 1:nphi
        slab = @view FDF[k, :, :]
        vals = slab[mask]
        out[k] = isempty(vals) ? NaN : mean(Float64.(vals))
    end
    out
end

spec1  = mean_spectrum(FDF, mask1)
spec10 = mean_spectrum(FDF, mask10)

# ------------------------------------------------------------
# COORDS
# ------------------------------------------------------------
x = range(0, LBOX_PC; length=nx)
y = range(0, LBOX_PC; length=ny)

xt = 0:10:50
yt = 0:10:50

# ------------------------------------------------------------
# PLOT (heatmap height == 2 spectra height)
# ------------------------------------------------------------
with_theme(theme_latexfonts()) do
    fig = Figure(size=(1250, 650))

    outer = fig[1, 1] = GridLayout()
    colgap!(outer, 18)

    left  = outer[1, 1] = GridLayout()   # heatmap + colorbar
    right = outer[1, 2] = GridLayout()   # spectra stacked

    # LEFT: heatmap fills full height
    axH = Axis(left[1, 1];
        xlabel="Distance [pc]",
        ylabel="Distance [pc]",
        aspect=DataAspect(),
        xticks=(xt, string.(xt)),
        yticks=(yt, string.(yt)),
        title = @sprintf("LOS: %s   q10=%.3g   q90=%.3g", LOS, q10, q90),
    )

    hm = if USE_COLRANGE
        heatmap!(axH, x, y, Pmax'; colorrange=(VMIN, VMAX))
    else
        heatmap!(axH, x, y, Pmax')
    end

    Colorbar(left[1, 2], hm, label="[K]")

    contour!(axH, x, y, Float32.(mask1)';  levels=[0.5], color=:black, linewidth=2)
    contour!(axH, x, y, Float32.(mask10)'; levels=[0.5], color=:red,   linewidth=2)

    colsize!(left, 1, Relative(0.86))
    colsize!(left, 2, Relative(0.14))

    # RIGHT: two spectra stacked, same x, different y
    rowgap!(right, 10)

    axS1 = Axis(right[1, 1];
        xlabelvisible=false,
        xticklabelsvisible=false,
        xgridvisible=false,
        ygridvisible=false,
        ylabel=L"⟨F(\phi)⟩\;[K]",
        title="Mean spectrum — 1st decile",
    )

    axS2 = Axis(right[2, 1];
        xlabel=L"\phi\;[\mathrm{rad}\,\mathrm{m}^{-2}]",
        xgridvisible=false,
        ygridvisible=false,
        ylabel=L"⟨F(\phi)⟩\;[K]",
        title="Mean spectrum — 10th decile",
    )

    linkxaxes!(axS1, axS2)

    lines!(axS1, phi, spec1,  color=:black, linewidth=2)
    lines!(axS2, phi, spec10, color=:red,   linewidth=2)

    rowsize!(right, 1, Relative(0.5))
    rowsize!(right, 2, Relative(0.5))

    colsize!(outer, 1, Relative(0.62))
    colsize!(outer, 2, Relative(0.38))

    if SAVE_FIG
        mkpath(OUT_DIR)
        outpath = joinpath(OUT_DIR, @sprintf("Pmax_deciles_meanSpectra_%s_%s.%s", SIMU_NAME, LOS, OUT_EXT))
        save(outpath, fig)
        @info "Saved: $outpath"
    end

    fig
end
