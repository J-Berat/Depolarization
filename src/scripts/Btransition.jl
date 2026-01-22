using FITSIO
using CairoMakie
using Printf
using Random
using Statistics
using LaTeXStrings
using Logging

include(joinpath(@__DIR__, "../src/utils/fits_utils.jl"))
include(joinpath(@__DIR__, "../src/utils/los_utils.jl"))

# ------------------------------------------------------------
# USER CHOICES
# ------------------------------------------------------------
const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS = "y"                                   # "x", "y", or "z"
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"

const OUT_PLOTS = "/Users/jb270005/Desktop/Depolarization_canals/sightline_plots"
const SAVE_PLOT = true
const SAVE_EXT  = "pdf"

const WRITE_LOG = true
const WRITE_CSV = true

const NPIX    = 256
const LBOX_PC = 50.0
const B_SCALE = 1000.0      # -> μG
const RNG_SEED = nothing    # set to 0 for reproducible

# Phi axis (given) : -10 -> 10 step 0.25
const PhiArray = collect(-10.0:0.25:10.0)   # rad/m^2
const NPHI = length(PhiArray)

# Typography
const LABSIZE  = 25
const TICKSIZE = 25
const ARROW_FONTSIZE = 16

# ------------------------------------------------------------
# Tight "one-reversal" boundaries
# ------------------------------------------------------------
const SMOOTH_B   = true
const SMOOTH_WIN = 5

const SIGN_EPS   = 0.0
const DERIV_TOL  = 0.05

const STOP_AT_NEXT_REVERSAL = true
const MAX_HALF_WIDTH_PIX = 18
const ISOLATE_FALLBACK = true

# Plot styling
const B_LINEWIDTH  = 3

const BOUNDS_STYLE = :dot
const BOUNDS_WIDTH = 1.2
const BOUNDS_ALPHA = 0.9

const DRAW_LABELS = true
const DRAW_ARROWS = true

# ne style
const NE_COLOR = (:black, 0.85)
const NE_LW    = 2.5
const NE_STYLE = :dashdot

# Scatter styling (inside windows only)
const SCAT_ALPHA  = 0.8
const SCAT_MARKER = :x
const SCAT_MSIZE  = 7

# Tick align (Makie): 1.0 -> inside, 0.0 -> outside (visually)
const TICKALIGN_IN = 1.0
const TICKALIGN_OUT = 0.0

# Heatmap sightline marker: full dashed lines (alpha -> 0.6)
const MAP_MARK_COLOR = (:white, 0.6)
const MAP_MARK_STYLE = :dash
const MAP_MARK_LW    = 2.5

# |B|=0 line style
const B0_COLOR = :blue
const B0_STYLE = :dash
const B0_LW    = 2.0

# FDF color
const FDF_COLOR = :black   # ou :green

count_reversals_in_window(ks::Vector{Int}, kL::Int, kR::Int) =
    count(k -> (k >= kL && k <= kR-1), ks)

function first_derivative(f::AbstractVector, Δ::Real)
    n = length(f)
    @assert n ≥ 3
    out = similar(f, n - 2)
    inv2Δ = 1 / (2Δ)
    @inbounds for i in 2:(n-1)
        out[i-1] = (f[i+1] - f[i-1]) * inv2Δ
    end
    return out
end

function tight_bounds_for_reversal(B::Vector{Float64}, Δ::Real, k0::Int, ks_all::Vector{Int};
                                   deriv_tol::Real=0.05,
                                   stop_at_next_reversal::Bool=true,
                                   max_half_width_pix::Int=0,
                                   isolate_fallback::Bool=true)

    n = length(B)
    d1 = first_derivative(B, Δ)
    i0 = clamp(k0 - 1, 1, length(d1))

    function includes_other_reversal(kL::Int, kR::Int)
        for k in ks_all
            if k != k0 && (k >= kL && k <= kR-1)
                return true
            end
        end
        return false
    end

    iL = i0
    while iL > 1
        if max_half_width_pix > 0 && (i0 - (iL - 1)) > max_half_width_pix
            break
        end
        if abs(d1[iL]) < deriv_tol
            break
        end
        if sign_eps(d1[iL-1]) != 0 && sign_eps(d1[iL]) != 0 &&
           sign_eps(d1[iL-1]) != sign_eps(d1[iL])
            break
        end
        if stop_at_next_reversal
            k_left_candidate = (iL - 1) + 1
            k_right_current  = (i0) + 1
            if includes_other_reversal(k_left_candidate, k_right_current)
                break
            end
        end
        iL -= 1
    end

    iR = i0
    while iR < length(d1)
        if max_half_width_pix > 0 && ((iR + 1) - i0) > max_half_width_pix
            break
        end
        if abs(d1[iR]) < deriv_tol
            break
        end
        if sign_eps(d1[iR]) != 0 && sign_eps(d1[iR+1]) != 0 &&
           sign_eps(d1[iR]) != sign_eps(d1[iR+1])
            break
        end
        if stop_at_next_reversal
            k_left_current    = (i0) + 1
            k_right_candidate = (iR + 1) + 1
            if includes_other_reversal(k_left_current, k_right_candidate)
                break
            end
        end
        iR += 1
    end

    kL = clamp(iL + 1, 1, n)
    kR = clamp(iR + 1, 1, n)

    if isolate_fallback && count_reversals_in_window(ks_all, kL, kR) > 1
        best = nothing
        for w in 0:(n-1)
            a = max(1, k0 - w)
            b = min(n, (k0 + 1) + w)
            if count_reversals_in_window(ks_all, a, b) == 1
                best = (a, b)
                break
            end
        end
        if best !== nothing
            kL, kR = best
        else
            kL, kR = max(1, k0), min(n, k0+1)
        end
    end

    return kL, kR
end

function merge_intervals(intervals::Vector{Tuple{Int,Int}})
    isempty(intervals) && return Tuple{Int,Int}[]
    ints = sort(intervals, by = first)
    merged = Tuple{Int,Int}[]
    a, b = ints[1]
    for (c, d) in ints[2:end]
        if c <= b + 1
            b = max(b, d)
        else
            push!(merged, (a, b))
            a, b = c, d
        end
    end
    push!(merged, (a, b))
    merged
end

function double_arrow!(ax, x1::Real, x2::Real, y::Real;
                       color=:black, linewidth=2,
                       head_dx::Real=0.8, head_dy::Real=1.0)
    lines!(ax, [x1, x2], [y, y]; color=color, linewidth=linewidth)
    poly!(ax, Point2f[(x1, y), (x1 + head_dx, y + head_dy), (x1 + head_dx, y - head_dy)]; color=color)
    poly!(ax, Point2f[(x2, y), (x2 - head_dx, y + head_dy), (x2 - head_dx, y - head_dy)]; color=color)
    return nothing
end

function write_log(path::String; kwargs...)
    open(path, "w") do io
        for (k, v) in kwargs
            println(io, "$(k) = $(v)")
        end
    end
end

function write_combined_csv(path::String;
                            simu::String,
                            los::String,
                            i1::Int,
                            j1::Int,
                            pmax1,
                            ks_all::Vector{Int},
                            windows::Vector{Tuple{Int,Int,Int}},
                            win_merged::Vector{Tuple{Int,Int}},
                            dist::AbstractVector)
    open(path, "w") do io
        println(io, "simu,los,i1,j1,pmax1,n_reversals_total,n_windows,n_windows_merged")
        println(io, "$(simu),$(los),$(i1),$(j1),$(pmax1),$(length(ks_all)),$(length(windows)),$(length(win_merged))")
        println(io)
        println(io, "simu,los,i,j,kmin,kmax,s_kmin_pc,s_kmax_pc,width_pc,npts")
        for (kmin, kmax) in win_merged
            sL = dist[kmin]
            sR = dist[kmax]
            width = sR - sL
            npts = kmax - kmin + 1
            println(io, "$(simu),$(los),$(i1),$(j1),$(kmin),$(kmax),$(sL),$(sR),$(width),$(npts))")
        end
    end
end

function extract_fdf_spectrum(FDFcube, i::Int, j::Int, nphi::Int)
    sz = size(FDFcube)
    dims_phi = findall(d -> sz[d] == nphi, 1:ndims(FDFcube))
    isempty(dims_phi) && error("No dimension in FDFcube matches NPHI=$nphi. size(FDFcube)=$sz")
    dφ = dims_phi[1]

    ndims(FDFcube) == 3 || error("Expected 3D FDF cube, got ndims=$(ndims(FDFcube)), size=$sz")

    if dφ == 3
        return vec(@view FDFcube[i, j, :])
    elseif dφ == 1
        return vec(@view FDFcube[:, i, j])
    elseif dφ == 2
        return vec(@view FDFcube[i, :, j])
    else
        error("Unexpected phi dimension index $dφ")
    end
end

# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
mkpath(OUT_PLOTS)
SIMU_DIR = joinpath(SIMU_ROOT, SIMU_NAME)

PMAX_FILE = joinpath(SIMU_DIR, LOS, "Synchrotron", "WithFaraday", "Pmax.fits")
@assert isfile(PMAX_FILE) "Missing: $PMAX_FILE"
Pmax = read_fits_array(PMAX_FILE)

q10 = quantile(vec(Pmax), 0.1)
q90 = quantile(vec(Pmax), 0.9)
idxs_1  = findall(Pmax .<= q10)
idxs_10 = findall(Pmax .>= q90)

rng = isnothing(RNG_SEED) ? Random.default_rng() : MersenneTwister(RNG_SEED)
(i1, j1)   = Tuple(rand(rng, idxs_1))
(i10, j10) = Tuple(rand(rng, idxs_10))

p1  = Pmax[i1,  j1]
p10 = Pmax[i10, j10]

bname, profile_fun = los_config(LOS)
BLOS_FILE = joinpath(SIMU_DIR, bname * ".fits")
@assert isfile(BLOS_FILE) "Missing: $BLOS_FILE"
Bcube = read_fits_array(BLOS_FILE)

prof1  = B_SCALE .* profile_fun(i1,  j1,  Bcube)
prof10 = B_SCALE .* profile_fun(i10, j10, Bcube)

prof1_use  = SMOOTH_B ? smooth_moving_average(prof1,  SMOOTH_WIN) : prof1
prof10_use = SMOOTH_B ? smooth_moving_average(prof10, SMOOTH_WIN) : prof10

# ne cube
NE_FILE = joinpath(SIMU_DIR, LOS, "Synchrotron", "ne.fits")
@assert isfile(NE_FILE) "Missing: $NE_FILE"
ne_cube = read_fits_array(NE_FILE)
@assert size(ne_cube) == (NPIX, NPIX, NPIX) "Unexpected ne size: $(size(ne_cube))"
ne1  = profile_fun(i1,  j1,  ne_cube)
ne10 = profile_fun(i10, j10, ne_cube)

dist = collect(range(0.0, LBOX_PC; length=NPIX))
Δ = dist[2] - dist[1]

FDF_FILE = joinpath(SIMU_DIR, LOS, "Synchrotron", "WithFaraday", "FDF.fits")
@assert isfile(FDF_FILE) "Missing: $FDF_FILE"
FDFcube = read_fits_array(FDF_FILE)

FDF1  = extract_fdf_spectrum(FDFcube, i1,  j1,  NPHI)
FDF10 = extract_fdf_spectrum(FDFcube, i10, j10, NPHI)

A1  = abs.(FDF1)
A10 = abs.(FDF10)

ks_all = reversal_indices(Vector{Float64}(prof1_use); eps=SIGN_EPS)

windows = Tuple{Int,Int,Int}[]
for k0 in ks_all
    kL, kR = tight_bounds_for_reversal(
        Vector{Float64}(prof1_use), Δ, k0, ks_all;
        deriv_tol=DERIV_TOL,
        stop_at_next_reversal=STOP_AT_NEXT_REVERSAL,
        max_half_width_pix=MAX_HALF_WIDTH_PIX,
        isolate_fallback=ISOLATE_FALLBACK
    )
    push!(windows, (k0, kL, kR))
end
sort!(windows, by = w -> w[2])

win_intervals = [(kL, kR) for (_, kL, kR) in windows]
win_merged = merge_intervals(win_intervals)

pdf_path = joinpath(OUT_PLOTS, @sprintf("sightline_LOS%s.%s", LOS, SAVE_EXT))
log_path = joinpath(OUT_PLOTS, @sprintf("sightline_LOS%s.log", LOS))
csv_path = joinpath(OUT_PLOTS, @sprintf("sightline_LOS%s.csv", LOS))

if WRITE_LOG
    write_log(log_path;
        SIMU_NAME=SIMU_NAME, LOS=LOS,
        NE_FILE=NE_FILE,
        PMAX_FILE=PMAX_FILE, BLOS_FILE=BLOS_FILE, FDF_FILE=FDF_FILE,
        RNG_SEED=RNG_SEED,
        decile1_pixel_i=i1, decile1_pixel_j=j1, decile1_Pmax=p1,
        decile10_pixel_i=i10, decile10_pixel_j=j10, decile10_Pmax=p10,
        n_reversals_total=length(ks_all),
        win_merged=win_merged
    )
end

if WRITE_CSV
    write_combined_csv(
        csv_path;
        simu=SIMU_NAME, los=LOS,
        i1=i1, j1=j1, pmax1=p1,
        ks_all=ks_all, windows=windows,
        win_merged=win_merged, dist=dist
    )
end

# ------------------------------------------------------------
# FIGURE
# ------------------------------------------------------------
fig = with_theme(theme_latexfonts()) do
    fig = Figure(size=(1250, 980))

    # ---- Heatmaps
    ax_map1 = Axis(fig[1, 1],
        title = LaTeXString(@sprintf("1^{\\mathrm{st}}\\ decile\\ (i=%d,\\ j=%d)", i1, j1)),
        titlesize = LABSIZE,
        xlabel = L"i", ylabel = L"j",
        xlabelsize = LABSIZE, ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE,
        xgridvisible = false, ygridvisible = false,
    )
    hm = heatmap!(ax_map1, Pmax; colormap=:magma)
    vlines!(ax_map1, [i1]; color=MAP_MARK_COLOR, linestyle=MAP_MARK_STYLE, linewidth=MAP_MARK_LW)
    hlines!(ax_map1, [j1]; color=MAP_MARK_COLOR, linestyle=MAP_MARK_STYLE, linewidth=MAP_MARK_LW)

    ax_map10 = Axis(fig[2, 1],
        title = LaTeXString(@sprintf("10^{\\mathrm{th}}\\ decile\\ (i=%d,\\ j=%d)", i10, j10)),
        titlesize = LABSIZE,
        xlabel = L"i", ylabel = L"j",
        xlabelsize = LABSIZE, ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE,
        xgridvisible = false, ygridvisible = false,
    )
    heatmap!(ax_map10, Pmax; colormap=:magma)
    vlines!(ax_map10, [i10]; color=MAP_MARK_COLOR, linestyle=MAP_MARK_STYLE, linewidth=MAP_MARK_LW)
    hlines!(ax_map10, [j10]; color=MAP_MARK_COLOR, linestyle=MAP_MARK_STYLE, linewidth=MAP_MARK_LW)

    Colorbar(fig[:, 2], hm; label=L"P_{\max}\,[\mathrm{K}]", labelsize=LABSIZE, ticklabelsize=TICKSIZE)

    # --------------------------------------------------------
    # Right panel TOP (decile 1): MUST create row 2 BEFORE rowsize!
    # --------------------------------------------------------
    gl_top = GridLayout()
    fig[1, 3] = gl_top

    ax1B = Axis(gl_top[1, 1],
        xlabel = L"\mathrm{Distance\ [pc]}",
        ylabel = L"B_{\parallel}\ [\mu\mathrm{G}]",
        xaxisposition = :top,
        xlabelsize = LABSIZE, ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE,
        xgridvisible = false, ygridvisible = false,
        xtickalign = TICKALIGN_OUT,
        ytickalign = TICKALIGN_OUT,
    )

    ax1F = Axis(gl_top[2, 1],
        xlabel = L"\phi\ [\mathrm{rad}\ \mathrm{m}^{-2}]",
        ylabel = L"|F(\phi)|\ [\mathrm{K}]",
        xaxisposition = :bottom,
        xlabelsize = LABSIZE, ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE,
        xgridvisible = false, ygridvisible = false,
        xtickalign = TICKALIGN_IN,
        ytickalign = TICKALIGN_IN,
    )

    rowgap!(gl_top, 0)
    rowsize!(gl_top, 1, Relative(0.80))
    rowsize!(gl_top, 2, Relative(0.20))

    # Twin axis for ne on the TOP B axis
    ax1ne = Axis(gl_top[1, 1],
        yaxisposition = :right,
        ylabel = L"n_e\ [\mathrm{cm}^{-3}]",
        xaxisposition = :top,
        ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE,
        yticklabelsize = TICKSIZE,
        backgroundcolor = :transparent,
        xgridvisible = false, ygridvisible = false,
        xtickalign = TICKALIGN_OUT,
        ytickalign = TICKALIGN_OUT,
    )
    linkxaxes!(ax1B, ax1ne)
    hidexdecorations!(ax1ne; ticks=true, ticklabels=true, grid=true)

    # ---- XTICKS discovery (TOP)
    xlims!(ax1B, 0.0, 50.0)
    ax1B.xticks = (0:10:50, string.(0:10:50))

    pB0 = hlines!(ax1B, 0.0; color=B0_COLOR, linestyle=B0_STYLE, linewidth=B0_LW, label=L"|B_{\parallel}|=0")
    pB  = lines!(ax1B, dist, prof1_use; linewidth=B_LINEWIDTH, label=L"B_{\parallel}")
    pNe = lines!(ax1ne, dist, ne1; linewidth=NE_LW, color=NE_COLOR, linestyle=NE_STYLE, label=L"n_e")

    # Windows + arrows + labels on ax1B
    for (kL, kR) in win_merged
        scatter!(ax1B, dist[kL:kR], prof1[kL:kR];
            color = (:black, SCAT_ALPHA), marker = SCAT_MARKER, markersize = SCAT_MSIZE,
            label = L"B_{\parallel}"
        )
    end

    if !isempty(win_merged)
        yminB, ymaxB = extrema(prof1_use)
        yrB = ymaxB - yminB
        head_dx = 0.6
        head_dy = 0.018 * yrB
        y_arrow = ymaxB - 0.03 * yrB
        y_text  = ymaxB - 0.01 * yrB

        first = true
        for (kL, kR) in win_merged
            sL, sR = dist[kL], dist[kR]
            lbl = first ? L"\mathrm{reversal\ windows}" : nothing
            first = false

            vlines!(ax1B, [sL, sR]; linewidth=BOUNDS_WIDTH, linestyle=BOUNDS_STYLE,
                    color=:black, alpha=BOUNDS_ALPHA, label=lbl)

            if DRAW_ARROWS
                double_arrow!(ax1B, sL, sR, y_arrow; color=:black, linewidth=2, head_dx=head_dx, head_dy=head_dy)
            end
            if DRAW_LABELS
                Δs = sR - sL
                text!(ax1B, 0.5*(sL+sR), y_text;
                    text = LaTeXString(@sprintf("%.2f\\,\\mathrm{pc}", Δs)),
                    align = (:center, :bottom), fontsize = ARROW_FONTSIZE, color = :black)
            end
        end
    end

    axislegend(ax1B, [pB, pNe, pB0], [L"B_{\parallel}", L"n_e", L"|B_{\parallel}|=0"];
        position=:lb, merge=true, unique=true, framevisible=false, backgroundcolor=:white, padding=(8,8,8,8)
    )

    # ---- XTICKS discovery (BOTTOM): -10..10 step 0.25 (ticks), labels every 2.5
    xlims!(ax1F, -10.0, 10.0)
    xtφ = -10.0:0.25:10.0
    ax1F.xticks = (xtφ, [mod(i, 10) == 1 ? @sprintf("%.1f", x) : "" for (i, x) in enumerate(xtφ)])

    lines!(ax1F, PhiArray, A1; linewidth=2, label=L"|F(\phi)|", color=FDF_COLOR)
    vlines!(ax1F, PhiArray[argmax(A1)]; linestyle=:dash, linewidth=2, label=L"\phi_{\mathrm{peak}}", color=FDF_COLOR)
    axislegend(ax1F; position=:lb, framevisible=false, backgroundcolor=:white, padding=(8,8,8,8), merge=true, unique=true)

    # --------------------------------------------------------
    # Right panel BOTTOM (decile 10): same structure
    # --------------------------------------------------------
    gl_bot = GridLayout()
    fig[2, 3] = gl_bot

    ax10B = Axis(gl_bot[1, 1],
        xlabel = L"\mathrm{Distance\ [pc]}",
        ylabel = L"B_{\parallel}\ [\mu\mathrm{G}]",
        xaxisposition = :top,
        xlabelsize = LABSIZE, ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE,
        xgridvisible = false, ygridvisible = false,
        xtickalign = TICKALIGN_OUT,
        ytickalign = TICKALIGN_OUT,
    )

    ax10F = Axis(gl_bot[2, 1],
        xlabel = L"\phi\ [\mathrm{rad}\ \mathrm{m}^{-2}]",
        ylabel = L"|F(\phi)|\ [\mathrm{K}]",
        xaxisposition = :bottom,
        xlabelsize = LABSIZE, ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE,
        xgridvisible = false, ygridvisible = false,
        xtickalign = TICKALIGN_IN,
        ytickalign = TICKALIGN_IN,
    )

    rowgap!(gl_bot, 0)
    rowsize!(gl_bot, 1, Relative(0.80))
    rowsize!(gl_bot, 2, Relative(0.20))

    ax10ne = Axis(gl_bot[1, 1],
        yaxisposition = :right,
        ylabel = L"n_e\ [\mathrm{cm}^{-3}]",
        xaxisposition = :top,
        ylabelsize = LABSIZE,
        xticklabelsize = TICKSIZE,
        yticklabelsize = TICKSIZE,
        backgroundcolor = :transparent,
        xgridvisible = false, ygridvisible = false,
        xtickalign = TICKALIGN_OUT,
        ytickalign = TICKALIGN_OUT,
    )
    linkxaxes!(ax10B, ax10ne)
    hidexdecorations!(ax10ne; ticks=true, ticklabels=true, grid=true)

    # ---- XTICKS discovery (TOP)
    xlims!(ax10B, 0.0, 50.0)
    ax10B.xticks = (0:10:50, string.(0:10:50))

    p10B0 = hlines!(ax10B, 0.0; color=B0_COLOR, linestyle=B0_STYLE, linewidth=B0_LW, label=L"|B_{\parallel}|=0")
    p10B  = lines!(ax10B, dist, prof10_use; linewidth=B_LINEWIDTH, label=L"B_{\parallel}")
    p10Ne = lines!(ax10ne, dist, ne10; linewidth=NE_LW, color=NE_COLOR, linestyle=NE_STYLE, label=L"n_e")

    axislegend(ax10B, [p10B, p10Ne, p10B0], [L"B_{\parallel}", L"n_e", L"|B_{\parallel}|=0"];
        position=:lb, merge=true, unique=true, framevisible=false, backgroundcolor=:white, padding=(8,8,8,8)
    )

    # ---- XTICKS discovery (BOTTOM): -10..10 step 0.25 (ticks), labels every 2.5
    xlims!(ax10F, -10.0, 10.0)
    xtφ2 = -10.0:0.25:10.0
    ax10F.xticks = (xtφ2, [mod(i, 10) == 1 ? @sprintf("%.1f", x) : "" for (i, x) in enumerate(xtφ2)])

    lines!(ax10F, PhiArray, A10; linewidth=2, label=L"|F(\phi)|", color=FDF_COLOR)
    vlines!(ax10F, PhiArray[argmax(A10)]; linestyle=:dash, linewidth=2, label=L"\phi_{\mathrm{peak}}", color=FDF_COLOR)
    axislegend(ax10F; position=:lb, framevisible=false, backgroundcolor=:white, padding=(8,8,8,8), merge=true, unique=true)

    fig
end

if SAVE_PLOT
    save(pdf_path, fig)
    @info "Saved plot" path=pdf_path
end
