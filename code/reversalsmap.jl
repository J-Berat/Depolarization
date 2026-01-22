using FITSIO
using CairoMakie
using Statistics
using Printf
using Base.Threads
using LaTeXStrings

# ------------------------------------------------------------
# USER CHOICES
# ------------------------------------------------------------

const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS       = "y"   # "x", "y", or "z"
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"

const NPIX    = 256
const LBOX_PC = 50.0
const B_SCALE = 1000.0      # -> μG

const SMOOTH_B   = true
const SMOOTH_WIN = 5
const SIGN_EPS   = 0.0

const TARGET_NREV = 24

const OUT_DIR   = "/Users/jb270005/Desktop/Depolarization_canals/sightline_plots"
const SAVE_PLOT = true
const SAVE_EXT  = "pdf"

const LABSIZE  = 25
const TICKSIZE = 25

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

read_fits_array(path::AbstractString) = FITS(path, "r") do f
    read(f[1])
end

function los_config(los::String)
    if los == "x"
        return ("Bx", (i, j, B) -> Array(@view B[:, i, j]))
    elseif los == "y"
        return ("By", (i, j, B) -> Array(@view B[i, :, j]))
    elseif los == "z"
        return ("Bz", (i, j, B) -> Array(@view B[i, j, :]))
    else
        error("LOS must be x/y/z")
    end
end

function smooth_ma(x::AbstractVector, w::Int)
    w = max(w, 1)
    w = isodd(w) ? w : (w + 1)
    n = length(x)
    if w == 1
        return collect(x)
    end
    y = similar(x, n)
    h = (w - 1) ÷ 2
    @inbounds for i in 1:n
        i1 = max(1, i - h)
        i2 = min(n, i + h)
        y[i] = mean(@view x[i1:i2])
    end
    return y
end

@inline function sgn(x::Real; eps::Real=0.0)
    if x > eps
        return 1
    elseif x < -eps
        return -1
    else
        return 0
    end
end

function reversal_indices(B::AbstractVector; eps::Real=0.0)
    idx = Int[]
    @inbounds for k in 1:(length(B)-1)
        s1 = sgn(B[k]; eps=eps)
        s2 = sgn(B[k+1]; eps=eps)
        if s1 != 0 && s2 != 0 && s1 != s2
            push!(idx, k)  # reversal between k and k+1
        end
    end
    return idx
end

# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

mkpath(OUT_DIR)

SIMU_DIR = joinpath(SIMU_ROOT, SIMU_NAME)

# load Pmax (optional but handy to show where the pixel is)
PMAX_FILE = joinpath(SIMU_DIR, LOS, "Synchrotron", "WithFaraday", "Pmax.fits")
@assert isfile(PMAX_FILE) "Missing: $PMAX_FILE"
Pmax = read_fits_array(PMAX_FILE)
@assert size(Pmax) == (NPIX, NPIX) "Unexpected Pmax size: $(size(Pmax))"

# load B cube
bname, profile_fun = los_config(LOS)
BLOS_FILE = joinpath(SIMU_DIR, bname * ".fits")
@assert isfile(BLOS_FILE) "Missing: $BLOS_FILE"
Bcube = read_fits_array(BLOS_FILE)
@assert size(Bcube) == (NPIX, NPIX, NPIX) "Unexpected $bname size: $(size(Bcube))"

# distance axis
dist = collect(range(0.0, LBOX_PC; length=NPIX))

# compute nrev map
nrev = Array{Int16}(undef, NPIX, NPIX)
@threads for j in 1:NPIX
    for i in 1:NPIX
        prof = B_SCALE .* profile_fun(i, j, Bcube)
        prof_use = SMOOTH_B ? smooth_ma(prof, SMOOTH_WIN) : prof
        nrev[i, j] = Int16(length(reversal_indices(prof_use; eps=SIGN_EPS)))
    end
end

# pick a pixel: exact TARGET_NREV if possible, else closest
idx_exact = findall(nrev .== TARGET_NREV)

i0::Int = 0
j0::Int = 0
found_exact = !isempty(idx_exact)

if found_exact
    (i0, j0) = Tuple(idx_exact[1])  # take the first (or change to rand if you want)
    println("✔ Found pixel with exactly $TARGET_NREV reversals at (i,j)=($i0,$j0)")
else
    # closest in absolute difference
    bestΔ = typemax(Int)
    bestI = 1
    bestJ = 1
    @inbounds for j in 1:NPIX, i in 1:NPIX
        d = abs(Int(nrev[i,j]) - TARGET_NREV)
        if d < bestΔ
            bestΔ = d
            bestI = i
            bestJ = j
            bestΔ == 0 && break
        end
    end
    i0, j0 = bestI, bestJ
    println("⚠ No pixel with exactly $TARGET_NREV reversals.")
    println("   Using closest pixel (i,j)=($i0,$j0) with nrev=$(nrev[i0,j0]).")
end

# extract the sightline
prof = B_SCALE .* profile_fun(i0, j0, Bcube)
prof_use = SMOOTH_B ? smooth_ma(prof, SMOOTH_WIN) : prof
ks = reversal_indices(prof_use; eps=SIGN_EPS)
n_here = length(ks)

# ------------------------------------------------------------
# PLOT (Pmax map + sightline)
# ------------------------------------------------------------

with_theme(theme_latexfonts()) do
    fig = Figure(size=(1400, 650))

    axm = Axis(fig[1, 1],
        title = L"P_{\max}\ \mathrm{(pixel\ selected)}",
        xlabel = L"i", ylabel = L"j",
        xlabelsize = LABSIZE, ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE,
        xgridvisible = false, ygridvisible = false,
    )
    hm = heatmap!(axm, Pmax; colormap=:viridis)
    scatter!(axm, [i0], [j0]; markersize=18, color=:white, strokecolor=:black, strokewidth=2)

    Colorbar(fig[1, 2], hm;
        label = L"P_{\max}\,[\mathrm{K}]",
        labelsize = LABSIZE,
        ticklabelsize = TICKSIZE,
    )

    ax = Axis(fig[1, 3],
        title = LaTeXString(@sprintf("Sightline in %s — (i,j)=(%d,%d) — Nrev=%d", bname, i0, j0, n_here)),
        xlabel = L"\mathrm{Distance\ [pc]}",
        ylabel = L"B_{\mathrm{LOS}}\,[\mu\mathrm{G}]",
        xlabelsize = LABSIZE, ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE,
        xgridvisible = false, ygridvisible = false,
    )

    lines!(ax, dist, prof_use; linewidth=3, color=:black)

    # mark reversals (between k and k+1 -> draw at dist[k])
    for k in ks
        vlines!(ax, dist[k]; linestyle=:dash, linewidth=1.6, color=:red, alpha=0.85)
    end

    text!(ax, 0.02, 0.95;
        text = LaTeXString(@sprintf("SMOOTH=%s (win=%d) \\quad SIGN\\_EPS=%.3f", string(SMOOTH_B), SMOOTH_WIN, SIGN_EPS)),
        space = :relative,
        align = (:left, :top),
        fontsize = 18,
        color = :black
    )

    xlims!(ax, 0.0, LBOX_PC)

    if SAVE_PLOT
        out = joinpath(OUT_DIR, @sprintf("Sightline_NrevTarget%d_LOS%s_%s_%s_i%d_j%d.%s",
                                         TARGET_NREV, LOS, SIMU_NAME, bname, i0, j0, SAVE_EXT))
        save(out, fig)
        println("Saved plot → $out")
    else
        fig
    end
end
