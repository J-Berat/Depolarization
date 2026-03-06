using FFTW
using FITSIO
using CairoMakie
using LaTeXStrings
using Statistics
using Printf
using Interpolations
include("src/code/fits_io.jl")

# ============================================================
# PATHS
# ============================================================
Q_in =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Qnu.fits"

U_in =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Unu.fits"

Pmax_nofilter_path =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Pmax.fits"

Q_in_phi =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/realFDF.fits"

U_in_phi =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/imagFDF.fits"

base_out =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/WithInstrument"

mkpath(base_out)

# --- B / density cubes (for LIC figure) ---
Bx_in   = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/Bx.fits"
By_in   = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/By.fits"
Bz_in   = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/Bz.fits"
dens_in = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/density.fits"

# ============================================================
# GRID / PHYSICAL SCALE
# ============================================================
n, m = 256, 256
Lbox_pc = 50.0
Δx = Lbox_pc / n
Δy = Lbox_pc / m

# fixed small-scale removal: remove [0, Lcut_small] pc
Lcut_small = 1.0

# large-scale removal variants: remove [Llarge, 50]
Llarge_list = [100.0, 50.0, 40.0, 20.0]

# Nyquist in cycles/pc
fNy_x = 1 / (2Δx)
fNy_y = 1 / (2Δy)
fNy   = min(fNy_x, fNy_y)

@info "Grid" n m Δx Δy fNy

# ============================================================
# RM-SYNTHESIS AXES
# ============================================================
νmin_MHz = 120.0
νmax_MHz = 167.0
Δν_MHz   = 0.098

nuArray  = (νmin_MHz:Δν_MHz:νmax_MHz) .* 1e6
PhiArray = collect(-10.0:0.25:10.0)

# ============================================================
# UTILITIES
# ============================================================
function apply_instrument_2d(img::AbstractMatrix, H::AbstractMatrix)
    real.(ifft(fft(img) .* H))
end

function apply_to_array_xy(data, H; n::Int=256, m::Int=256)
    nd = ndims(data)

    if nd == 2
        return apply_instrument_2d(data, H)

    elseif nd == 3
        sz = size(data)
        out = similar(float.(data))

        if sz[1] == n && sz[2] == m
            @views for k in 1:sz[3]
                out[:, :, k] = apply_instrument_2d(data[:, :, k], H)
            end
            return out

        elseif sz[2] == n && sz[3] == m
            @views for k in 1:sz[1]
                tmp = apply_instrument_2d(reshape(data[k, :, :], n, m), H)
                out[k, :, :] .= reshape(tmp, 1, n, m)
            end
            return out

        else
            error("Unsupported 3D shape $(sz)")
        end
    else
        error("Unsupported ndims(data)=$nd")
    end
end

# ============================================================
# LIC utilities
# ============================================================
function interpol2d(map, xmap, ymap)
    nx, ny = size(map)
    x = range(1, nx, length=nx)
    y = range(1, ny, length=ny)
    itp = linear_interpolation((x,y), map, extrapolation_bc=Line())
    mapout = zeros(nx, ny)
    for p in 1:length(map)
        mapout[p] = itp(xmap[p], ymap[p])
    end
    return mapout
end

function xymap(nx, ny)
    xmap = collect(1:nx) .* ones(Int64, ny)'
    ymap = ones(Int64, nx) .* collect(1:ny)'
    return xmap, ymap
end

function lic(vx, vy; niter=1, len=8, normalize=false, amplitude=nothing, level=0.1, scalar=nothing)
    nx, ny = size(vx)

    uu = @. sqrt(vx^2 + vy^2)
    ii = findall(uu .== 0.0)
    if !isempty(ii)
        uu[ii] .= 1.0
    end

    if normalize
        ux = vx ./ uu
        uy = vy ./ uu
    else
        ux = vx ./ maximum(uu)
        uy = vy ./ maximum(uu)
    end

    vl = rand(nx, ny)

    for it in 1:niter
        texture = copy(vl)

        vv = zeros(nx, ny)
        pi, pj = xymap(nx, ny)
        mi = copy(pi)
        mj = copy(pj)

        ppi = 1.0*pi
        ppj = 1.0*pj
        mmi = 1.0*mi
        mmj = 1.0*mj

        for l in 0:len
            dpi = interpol2d(ux, ppi, ppj)
            dpj = interpol2d(uy, ppi, ppj)
            dmi = interpol2d(ux, mmi, mmj)
            dmj = interpol2d(uy, mmi, mmj)

            ppi = @. ppi + 0.25 * dpi
            ppj = @. ppj + 0.25 * dpj
            mmi = @. mmi - 0.25 * dmi
            mmj = @. mmj - 0.25 * dmj

            pi = @. mod(round(ppi) + nx-1, nx) + 1
            pj = @. mod(round(ppj) + ny-1, ny) + 1
            mi = @. mod(round(mmi) + nx-1, nx) + 1
            mj = @. mod(round(mmj) + ny-1, ny) + 1

            ppi = @. pi + (ppi - round(ppi))
            ppj = @. pj + (ppj - round(ppj))
            mmi = @. mi + (mmi - round(mmi))
            mmj = @. mj + (mmj - round(mmj))

            tp = interpol2d(texture, ppi, ppj)
            tm = interpol2d(texture, mmi, mmj)
            vv = @. vv + tp + tm
        end

        vl = 0.25 * vv / len
    end

    if amplitude !== nothing
        if scalar !== nothing && length(scalar) > 1
            uu .= scalar
        end

        if length(level) != 1
            level = 0.1
        end

        if !isempty(ii)
            uu[ii] .= 0.0
        end

        level = max([0.0, min([level, 1.0])])
        uu .= (1.0 - level) * uu ./ maximum(uu) .+ level
        vl .= vl .* uu
    end

    if !isempty(ii)
        vl[ii] .= 0.0
    end

    return vl
end

Pphi_max_map(absF::AbstractArray) =
    dropdims(maximum(absF; dims=ndims(absF)); dims=ndims(absF))

# ============================================================
# PSD utilities (isotropic 1D: <|FFT|^2>(k))
# ============================================================
function psd2d(img::AbstractMatrix)
    A = float.(img)
    A .-= mean(A)
    F = fft(A)
    P = abs2.(F) ./ (length(A)^2)
    return fftshift(P)
end

function psd1d_isotropic(img::AbstractMatrix, kx::AbstractVector, ky::AbstractVector;
                         nbins::Int=80, kmin::Real=0.0, kmax::Real=Inf)

    n, m = size(img)
    @assert length(kx) == n
    @assert length(ky) == m

    P2 = psd2d(img)

    KX = repeat(kx, 1, m)
    KY = repeat(ky', n, 1)
    K  = sqrt.(KX.^2 .+ KY.^2)

    kvec = vec(K)
    pvec = vec(P2)

    good = isfinite.(kvec) .& isfinite.(pvec) .& (kvec .>= kmin) .& (kvec .<= kmax)
    kvec = kvec[good]
    pvec = pvec[good]

    kmax_eff = isempty(kvec) ? 0.0 : maximum(kvec)
    edges = range(max(kmin, 0.0), min(kmax, kmax_eff), length=nbins+1)
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])

    Pk = fill(NaN, nbins)
    @inbounds for i in 1:nbins
        lo, hi = edges[i], edges[i+1]
        sel = (kvec .>= lo) .& (kvec .< hi)
        if any(sel)
            Pk[i] = mean(pvec[sel])
        end
    end

    return centers, Pk
end

# ============================================================
# HARD BAND-PASS MASK (cycles/pc) for Fourier-space filtering
# ============================================================
function instrument_bandpass_L(n::Int, m::Int;
                               Δx::Real, Δy::Real=Δx,
                               Lcut_small::Real,
                               Llarge::Real,
                               fNy::Real)

    fx = fftfreq(n, Δx)
    fy = fftfreq(m, Δy)
    FX = repeat(fx, 1, m)
    FY = repeat(fy', n, 1)
    F  = sqrt.(FX.^2 .+ FY.^2)

    flo = 1 / Llarge
    fhi_raw = 1 / Lcut_small
    fhi = min(fhi_raw, fNy)

    @info "Filter" Lcut_small Llarge flo fhi fhi_raw

    H = Float32.((F .>= flo) .& (F .<= fhi))
    return H, fftshift(H)
end

# ============================================================
# EXTRA: 1D spectrum along x averaged over y (Fourier-space)
# ============================================================
function psd1d_x_mean_over_y(img::AbstractMatrix, kx::AbstractVector; remove_mean::Bool=true)
    n, m = size(img)
    @assert length(kx) == n

    A = float.(img)
    if remove_mean
        A .-= mean(A)
    end

    Fx  = fft(A, 1)
    Px  = abs2.(Fx) ./ (n^2)
    Pkx = vec(mean(Px; dims=2))
    return kx, fftshift(Pkx)
end

# ============================================================
# helper: extract channel slice (n,m) from cube (n,m,nν) or (nν,n,m)
# ============================================================
function get_chan_xy(cube, ichan::Int, n::Int, m::Int)
    @assert ndims(cube) == 3 "cube must be 3D. Got ndims=$(ndims(cube)) size=$(size(cube))"
    sz = size(cube)

    if sz[1] == n && sz[2] == m
        @assert 1 ≤ ichan ≤ sz[3] "ichan=$(ichan) out of bounds for third dim=$(sz[3])"
        return float.(@view cube[:, :, ichan])
    elseif sz[2] == n && sz[3] == m
        @assert 1 ≤ ichan ≤ sz[1] "ichan=$(ichan) out of bounds for first dim=$(sz[1])"
        return float.(reshape(@view(cube[ichan, :, :]), n, m))
    else
        error("Unsupported cube shape $(sz). Expected (n,m,nν) or (nν,n,m) with n=$n m=$m.")
    end
end

# ============================================================
# LaTeX tick labels for linear axes
# ============================================================
function set_latex_linear_ticks_pc!(ax::Axis; step_pc::Real=10.0, Lbox_pc::Real=50.0)
    vals = collect(0.0:step_pc:Lbox_pc)
    ax.xticks = vals
    ax.yticks = vals
    ax.xtickformat = (vs)->[LaTeXString("\\mathrm{$(round(Int, v))}") for v in vs]
    ax.ytickformat = (vs)->[LaTeXString("\\mathrm{$(round(Int, v))}") for v in vs]
    return nothing
end

# ============================================================
# POW10 ticks for log axes
# ============================================================
pow10_label_from_val(v::Real) = LaTeXString("10^{$(round(Int, log10(v)))}")

function set_pow10_ticks!(ax::Axis; which::Symbol = :both)
    r = ax.finallimits[]
    xmin = r.origin[1]
    ymin = r.origin[2]
    xmax = r.origin[1] + r.widths[1]
    ymax = r.origin[2] + r.widths[2]

    if (which == :x || which == :both) && xmin > 0 && xmax > 0
        emin = floor(Int, log10(xmin))
        emax = ceil(Int,  log10(xmax))
        vals = 10.0 .^ collect(emin:emax)
        ax.xticks = vals
        ax.xtickformat = (vs)->[pow10_label_from_val(v) for v in vs]
    end

    if (which == :y || which == :both) && ymin > 0 && ymax > 0
        emin = floor(Int, log10(ymin))
        emax = ceil(Int,  log10(ymax))
        vals = 10.0 .^ collect(emin:emax)
        ax.yticks = vals
        ax.ytickformat = (vs)->[pow10_label_from_val(v) for v in vs]
    end
    return nothing
end

# ============================================================
# Safe log scaling
# ============================================================
function set_log_safe!(ax::Axis; xdata::AbstractVector, ydata::AbstractVector, pad_frac::Real=0.08)
    xok = isfinite.(xdata) .& (xdata .> 0)
    yok = isfinite.(ydata) .& (ydata .> 0)
    @assert any(xok) "No positive finite x data for log axis"
    @assert any(yok) "No positive finite y data for log axis"

    xmin = minimum(xdata[xok]); xmax = maximum(xdata[xok])
    ymin = minimum(ydata[yok]); ymax = maximum(ydata[yok])

    xmin = xmin * (1 - pad_frac)
    xmax = xmax * (1 + pad_frac)
    ymin = ymin * (1 - pad_frac)
    ymax = ymax * (1 + pad_frac)

    xmin = max(xmin, eps(Float64))
    ymin = max(ymin, eps(Float64))

    ax.xscale = log10
    ax.yscale = log10
    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)
    return nothing
end

# ============================================================
# vertical dashed lines
# ============================================================
function add_filter_verticals_kx!(ax; Llist_pc::AbstractVector, linestyle=:dash)
    ks = (2π) ./ Float64.(Llist_pc)
    vlines!(ax, ks; linestyle=linestyle, color=:black)
    return nothing
end

# ============================================================
# Strongly distinctive colors
# ============================================================
function curve_colors(n::Int)
    base = [
        :black,
        :red3,
        :dodgerblue3,
        :seagreen3,
        :darkorange2,
        :purple3,
        :magenta3,
        :goldenrod2,
        :brown3,
    ]
    @assert n ≤ length(base) "Not enough predefined distinct colors for n=$n"
    return base[1:n]
end

# ============================================================
# Peak helpers
# ============================================================
function kpeak_in_window(k::AbstractVector, y::AbstractVector; kmin::Real=0.0, kmax::Real=Inf)
    ok = isfinite.(k) .& isfinite.(y) .& (k .> 0) .& (y .> 0) .& (k .>= kmin) .& (k .<= kmax)
    if !any(ok)
        return NaN
    end
    kk = k[ok]
    yy = y[ok]
    return kk[argmax(yy)]
end

function add_peak_vline!(ax, k::AbstractVector, y::AbstractVector;
                         kmin::Real=0.0, kmax::Real=Inf,
                         color=:black, linestyle=:dot, linewidth=4.5)
    kp = kpeak_in_window(k, y; kmin=kmin, kmax=kmax)
    if isfinite(kp)
        vlines!(ax, [kp]; color=color, linestyle=linestyle, linewidth=linewidth)
    end
    return kp
end

# ============================================================
# helper: Llabel where you asked round(50/256*Llarge)
# ============================================================
Llabel_pc(Llarge) = round(50/256 * Llarge)

# ============================================================
# AXES
# ============================================================
kx = (2π) .* fftshift(fftfreq(n, Δx))
ky = (2π) .* fftshift(fftfreq(m, Δy))
x_pc = (0:n-1) .* Δx
y_pc = (0:m-1) .* Δy

# ============================================================
# Colorbar ticklabels in LaTeX (GLOBAL)
# ============================================================
cb_latex_ticks(vs) = [LaTeXString("\\mathrm{$(Printf.@sprintf("%.2g", v))}") for v in vs]

# ============================================================
# THEMES
# ============================================================
function set_theme_heatmaps!()
    set_theme!(Theme(
        fontsize = 44,
        Axis = (
            titlesize      = 40,
            xlabelsize     = 56,
            ylabelsize     = 56,
            xticklabelsize = 46,
            yticklabelsize = 46,
            spinewidth     = 2,
            xgridvisible   = false,
            ygridvisible   = false,
            xminorticksvisible = false,
            yminorticksvisible = false,
        ),
        Colorbar = (
            labelsize     = 68,
            ticklabelsize = 44,
            width         = 34,
            tickformat    = cb_latex_ticks,
        ),
    ))
end

function set_theme_spectra!()
    set_theme!(Theme(
        fontsize = 44,
        Axis = (
            titlesize      = 40,
            xlabelsize     = 56,
            ylabelsize     = 56,
            xticklabelsize = 46,
            yticklabelsize = 46,
            spinewidth     = 2,

            xgridvisible   = false,
            ygridvisible   = false,

            xminorticksvisible = true,
            yminorticksvisible = true,
            xminorticks = IntervalsBetween(9),
            yminorticks = IntervalsBetween(9),

            xticksize      = 14,
            yticksize      = 14,
            xtickwidth     = 2,
            ytickwidth     = 2,
            xminorticksize = 9,
            yminorticksize = 9,
            xminortickwidth = 2,
            yminortickwidth = 2,
        ),
        Legend = (labelsize = 44,),
        Colorbar = (
            labelsize     = 62,
            ticklabelsize = 44,
            width         = 34,
            tickformat    = cb_latex_ticks,
        ),
    ))
end

# ============================================================
# READ RAW Q/U
# ============================================================
Qdata = read_FITS(Q_in)
Udata = read_FITS(U_in)

# ============================================================
# OPTION A (I/O OPTI): cache only the channel slice ichan in RAM
# ============================================================
ichan = 50

Qslice_filt = Dict{Float64, Matrix{Float32}}()
Uslice_filt = Dict{Float64, Matrix{Float32}}()

# ============================================================
# FILTER Q/U → RM-SYNTHESIS → Pmax
# ============================================================
for Llarge in Llarge_list
    tag = @sprintf("HardBandPass_remove_L0_to_1pc_and_%dto50pc", Int(round(Llarge)))
    subdir = joinpath(base_out, tag)
    mkpath(subdir)

    H, _ = instrument_bandpass_L(n, m;
                                 Δx=Δx, Δy=Δy,
                                 Lcut_small=Lcut_small,
                                 Llarge=Llarge,
                                 fNy=fNy)

    Qf = apply_to_array_xy(Qdata, H; n=n, m=m)
    Uf = apply_to_array_xy(Udata, H; n=n, m=m)

    Qslice_filt[Llarge] = Float32.(get_chan_xy(Qf, ichan, n, m))
    Uslice_filt[Llarge] = Float32.(get_chan_xy(Uf, ichan, n, m))

    write_FITS(joinpath(subdir, "Qnu_filtered.fits"), Qf)
    write_FITS(joinpath(subdir, "Unu_filtered.fits"), Uf)

    absF, _, _ = RMSynthesis(Qf, Uf, nuArray, PhiArray; log_progress=true)

    Pmaxf = Pphi_max_map(absF)
    mkpath(joinpath(subdir, "RMSynthesis"))
    write_FITS(joinpath(subdir, "RMSynthesis", "Pphi_max.fits"), Pmaxf)
end

# ============================================================
# LOAD Pmax RESULTS
# ============================================================
Pmax0 = read_FITS(Pmax_nofilter_path)

Pmax_filt = Dict{Float64, Matrix}()
L_ok      = Float64[]

for Llarge in Llarge_list
    tag = @sprintf("HardBandPass_remove_L0_to_1pc_and_%dto50pc", Int(round(Llarge)))
    p = joinpath(base_out, tag, "RMSynthesis", "Pphi_max.fits")
    isfile(p) || continue
    Pmax_filt[Llarge] = read_FITS(p)
    push!(L_ok, Llarge)
end

@assert !isempty(L_ok) "No filtered Pphi_max.fits found under: $base_out"

# ============================================================
# COLOR RANGES (robust percentiles)
# ============================================================
vals = vec(float.(Pmax0))
for Llarge in L_ok
    append!(vals, vec(float.(Pmax_filt[Llarge])))
end
vals = vals[isfinite.(vals)]
sort!(vals)
lo = vals[clamp(round(Int, 0.01length(vals)), 1, length(vals))]
hi = vals[clamp(round(Int, 0.99length(vals)), 1, length(vals))]
pmax_colorrange = (lo, hi)

# ============================================================
# FIGURE 1: Pmax maps
# ============================================================
with_theme(theme_latexfonts()) do
    set_theme_heatmaps!()

    fig = Figure(size=(3000, 1050), figure_padding=30)

    ax0 = Axis(fig[1,1],
        xlabel = LaTeXString("x\\,[\\mathrm{pc}]"),
        ylabel = LaTeXString("y\\,[\\mathrm{pc}]"),
        aspect = DataAspect()
    )
    hm0 = heatmap!(ax0, x_pc, y_pc, Pmax0; colorrange=pmax_colorrange, colormap=:magma)
    set_latex_linear_ticks_pc!(ax0; step_pc=10.0, Lbox_pc=Lbox_pc)

    text!(ax0, 0.97, 0.97;
        text  = LaTeXString("\\mathrm{no\\ filter}"),
        space = :relative,
        align = (:right, :top),
        color = :white,
        fontsize = 50
    )

    for (j, Llarge) in enumerate(L_ok)
        ax = Axis(fig[1, j+1],
            xlabel = LaTeXString("x\\,[\\mathrm{pc}]"),
            aspect = DataAspect(),
            ylabelvisible = false,
            yticklabelsvisible = false,
        )
        heatmap!(ax, x_pc, y_pc, Pmax_filt[Llarge]; colorrange=pmax_colorrange, colormap=:magma)
        set_latex_linear_ticks_pc!(ax; step_pc=10.0, Lbox_pc=Lbox_pc)

        text!(ax, 0.97, 0.97;
            text  = LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"),
            space = :relative,
            align = (:right, :top),
            color = :white,
            fontsize = 50
        )
    end

    Colorbar(fig[1, length(L_ok)+2], hm0;
        label      = LaTeXString("P_{\\max}\\,[\\mathrm{K}]"),
        tellheight = true,
        tickformat = cb_latex_ticks
    )

    display(fig)
end

# ============================================================
# FIGURE 2 : PSDs ONLY
# ============================================================
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figPSD = Figure(size=(3000, 950), figure_padding=30)

    nbins_psd = 90
    kmin_plot = 2π * (1 / Lbox_pc)
    kmax_plot = maximum(sqrt.(repeat(kx,1,m).^2 .+ repeat(ky',n,1).^2))

    kcen0, Pk0 = psd1d_isotropic(Pmax0, kx, ky; nbins=nbins_psd, kmin=kmin_plot, kmax=kmax_plot)
    ok0 = isfinite.(Pk0) .& (kcen0 .> 0) .& (Pk0 .> 0)

    axpsd0 = Axis(figPSD[1,1],
        xlabel = LaTeXString("k\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\langle |\\tilde{P}|^2 \\rangle(k)"),
    )
    lines!(axpsd0, kcen0[ok0], Pk0[ok0])
    set_log_safe!(axpsd0; xdata=kcen0[ok0], ydata=Pk0[ok0])
    autolimits!(axpsd0)
    set_pow10_ticks!(axpsd0; which=:both)

    for (j, Llarge) in enumerate(L_ok)
        kcen, Pk = psd1d_isotropic(Pmax_filt[Llarge], kx, ky; nbins=nbins_psd, kmin=kmin_plot, kmax=kmax_plot)
        ok = isfinite.(Pk) .& (kcen .> 0) .& (Pk .> 0)

        axpsd = Axis(figPSD[1, j+1],
            xlabel = LaTeXString("k\\,[\\mathrm{rad\\,pc^{-1}}]"),
            ylabelvisible = false,
            yticklabelsvisible = false,
        )
        lines!(axpsd, kcen[ok], Pk[ok])
        set_log_safe!(axpsd; xdata=kcen[ok], ydata=Pk[ok])
        autolimits!(axpsd)
        set_pow10_ticks!(axpsd; which=:both)
    end

    display(figPSD)
end

# ============================================================
# FIGURE 3 : Pmax spectrum along x, mean over y
# ============================================================
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figS = Figure(size=(1500, 1200), figure_padding=30)
    axS = Axis(figS[1, 1],
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\langle |\\tilde{P}(k_x, y)|^2 \\rangle_y"),
    )

    handles = Any[]
    labels  = Any[]

    Lsorted = sort(L_ok)
    cols = curve_colors(1 + length(Lsorted))

    kxv0, P0 = psd1d_x_mean_over_y(Pmax0, kx; remove_mean=true)
    ok0 = isfinite.(P0) .& (kxv0 .> 0) .& (P0 .> 0)
    l0 = lines!(axS, kxv0[ok0], P0[ok0], color=cols[1], linewidth=3)
    push!(handles, l0)
    push!(labels, LaTeXString("\\mathrm{no\\ filter}"))

    for (i, Llarge) in enumerate(Lsorted)
        kxv, Pk = psd1d_x_mean_over_y(Pmax_filt[Llarge], kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)
        li = lines!(axS, kxv[ok], Pk[ok], color=cols[i+1], linewidth=3)
        push!(handles, li)
        push!(labels, LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"))
    end

    add_filter_verticals_kx!(axS; Llist_pc=Llarge_list, linestyle=:dash)
    axislegend(axS, handles, labels; position=:rt, framevisible=true)

    allx = Float64[]
    ally = Float64[]
    append!(allx, kxv0[ok0]); append!(ally, P0[ok0])
    for Llarge in Lsorted
        kxv, Pk = psd1d_x_mean_over_y(Pmax_filt[Llarge], kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)
        append!(allx, kxv[ok]); append!(ally, Pk[ok])
    end

    set_log_safe!(axS; xdata=allx, ydata=ally)
    autolimits!(axS)
    set_pow10_ticks!(axS; which=:both)

    # ---------------- INSET: k_peak vs Llarge (Pmax space) ----------------
    let
        kmin_win = 0.12
        kmax_win = Inf
        @info "Inset Pmax: k-window" kmin_win kmax_win

        Leff  = Float64[]
        kpeak = Float64[]
        cpeak = Symbol[]

        for (i, Llarge) in enumerate(Lsorted)
            kxv, Pk = psd1d_x_mean_over_y(Pmax_filt[Llarge], kx; remove_mean=true)
            kp = kpeak_in_window(kxv, Pk; kmin=kmin_win, kmax=kmax_win)

            push!(Leff, (50/256) * Llarge)
            push!(kpeak, kp)
            push!(cpeak, cols[i+1])

            @info "Inset Pmax peak" Llarge kp
        end

        ok = isfinite.(Leff) .& isfinite.(kpeak) .& (kpeak .> 0)

        axInset = Axis(figS[1,1],
            width  = Relative(0.38),
            height = Relative(0.34),
            halign = :left,
            valign = :bottom,

            xlabel = LaTeXString("L_{\\mathrm{large}}"),
            ylabel = LaTeXString("k_{\\mathrm{x,peak}}"),

            xminorticksvisible = false,
            yminorticksvisible = true,
            yminorticks = IntervalsBetween(9),
            yminorticksize = 5,
            yminortickwidth = 1.5,

            xlabelsize = 25,
            ylabelsize = 25,
            xticksize = 8,
            yticksize = 8,
            xtickwidth = 2,
            ytickwidth = 2,
            xticklabelsize = 22,
            yticklabelsize = 22,

            xtickalign = 1,
            ytickalign = 1,

            xaxisposition = :top,
            yaxisposition = :right,
        )

        lines!(axInset, Leff[ok], kpeak[ok], color=:black, linewidth=2)
        scatter!(axInset, Leff[ok], kpeak[ok]; color=cpeak[ok], markersize=14)

        xlims!(axInset, 0, 40)

        ymin = max(minimum(kpeak[ok]) * 0.9, eps(Float64))
        ymax = maximum(kpeak[ok]) * 1.1
        axInset.yscale = linear
        ylims!(axInset, ymin, ymax)

        axInset.xtickformat = (vs)->[
            LaTeXString("\\mathrm{$(round(v; digits=1))}") for v in vs
        ]
        set_pow10_ticks!(axInset; which=:y)

        axInset.backgroundcolor = (:white, 0.85)
    end

    display(figS)
end

# ============================================================
# FIGURES 4–7 : Q, U, P, Q^2 at specific channel
# ============================================================

# ---- FIGURE 4: Qν channel ichan ----
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figQ = Figure(size=(1500, 1200), figure_padding=30)
    axQ = Axis(figQ[1, 1],
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\langle |\\tilde{Q}(k_x, y;\\nu_{50})|^2 \\rangle_y"),
    )

    handles = Any[]
    labels  = Any[]

    cols = curve_colors(1 + length(sort(L_ok)))

    Qslice0 = get_chan_xy(Qdata, ichan, n, m)
    kxv0, P0 = psd1d_x_mean_over_y(Qslice0, kx; remove_mean=true)
    ok0 = isfinite.(P0) .& (kxv0 .> 0) .& (P0 .> 0)
    l0 = lines!(axQ, kxv0[ok0], P0[ok0], color=cols[1], linewidth=3)
    push!(handles, l0)
    push!(labels, LaTeXString("\\mathrm{no\\ filter}"))

    allx = Float64[]; ally = Float64[]
    append!(allx, kxv0[ok0]); append!(ally, P0[ok0])

    for (i, Llarge) in enumerate(sort(L_ok))
        haskey(Qslice_filt, Llarge) || continue
        Qslice = Qslice_filt[Llarge]

        kxv, Pk = psd1d_x_mean_over_y(Qslice, kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)

        li = lines!(axQ, kxv[ok], Pk[ok], color=cols[i+1], linewidth=3)
        push!(handles, li)
        push!(labels, LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"))

        append!(allx, kxv[ok]); append!(ally, Pk[ok])
    end

    add_filter_verticals_kx!(axQ; Llist_pc=Llarge_list, linestyle=:dash)
    axislegend(axQ, handles, labels; position=:rt, framevisible=true)

    set_log_safe!(axQ; xdata=allx, ydata=ally)
    autolimits!(axQ)
    set_pow10_ticks!(axQ; which=:both)

    display(figQ)
end

# ---- FIGURE 5: Uν channel ichan ----
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figU = Figure(size=(1500, 1200), figure_padding=30)
    axU = Axis(figU[1, 1],
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\langle |\\tilde{U}(k_x, y;\\nu_{50})|^2 \\rangle_y"),
    )

    handles = Any[]
    labels  = Any[]

    cols = curve_colors(1 + length(sort(L_ok)))

    Uslice0 = get_chan_xy(Udata, ichan, n, m)
    kxv0, P0 = psd1d_x_mean_over_y(Uslice0, kx; remove_mean=true)
    ok0 = isfinite.(P0) .& (kxv0 .> 0) .& (P0 .> 0)
    l0 = lines!(axU, kxv0[ok0], P0[ok0], color=cols[1], linewidth=3)
    push!(handles, l0)
    push!(labels, LaTeXString("\\mathrm{no\\ filter}"))

    allx = Float64[]; ally = Float64[]
    append!(allx, kxv0[ok0]); append!(ally, P0[ok0])

    for (i, Llarge) in enumerate(sort(L_ok))
        haskey(Uslice_filt, Llarge) || continue
        Uslice = Uslice_filt[Llarge]

        kxv, Pk = psd1d_x_mean_over_y(Uslice, kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)

        li = lines!(axU, kxv[ok], Pk[ok], color=cols[i+1], linewidth=3)
        push!(handles, li)
        push!(labels, LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"))

        append!(allx, kxv[ok]); append!(ally, Pk[ok])
    end

    add_filter_verticals_kx!(axU; Llist_pc=Llarge_list, linestyle=:dash)
    axislegend(axU, handles, labels; position=:rt, framevisible=true)

    set_log_safe!(axU; xdata=allx, ydata=ally)
    autolimits!(axU)
    set_pow10_ticks!(axU; which=:both)

    display(figU)
end

# ---- FIGURE 6: Pν channel ichan where P = sqrt(Q^2 + U^2) ----
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figP = Figure(size=(1500, 1200), figure_padding=30)
    axP = Axis(figP[1, 1],
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\langle |\\tilde{P}(k_x, y;\\nu_{50})|^2 \\rangle_y"),
    )

    handles = Any[]
    labels  = Any[]

    kmin_win = 0.12
    kmax_win = Inf

    Lsorted = sort(L_ok)
    cols = curve_colors(1 + length(Lsorted))

    specs = Vector{NamedTuple}()

    # no filter
    Q0 = get_chan_xy(Qdata, ichan, n, m)
    U0 = get_chan_xy(Udata, ichan, n, m)
    P0img = sqrt.(Q0.^2 .+ U0.^2)

    kxv0, P0spec = psd1d_x_mean_over_y(P0img, kx; remove_mean=true)
    ok0 = isfinite.(P0spec) .& (kxv0 .> 0) .& (P0spec .> 0)
    push!(specs, (label=LaTeXString("\\mathrm{no\\ filter}"), color=cols[1], kxv=kxv0, Pk=P0spec, ok=ok0))

    # filtered
    for (i, Llarge) in enumerate(Lsorted)
        (haskey(Qslice_filt, Llarge) && haskey(Uslice_filt, Llarge)) || continue
        Qs = Float64.(Qslice_filt[Llarge])
        Us = Float64.(Uslice_filt[Llarge])
        Pimg = sqrt.(Qs.^2 .+ Us.^2)

        kxv, Pk = psd1d_x_mean_over_y(Pimg, kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)

        push!(specs, (
            label=LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"),
            color=cols[i+1],
            kxv=kxv,
            Pk=Pk,
            ok=ok
        ))
    end

    allx = Float64[]
    ally = Float64[]
    for s in specs
        append!(allx, s.kxv[s.ok])
        append!(ally, s.Pk[s.ok])
    end

    kmax_band = maximum(allx)
    vspan!(axP, kmin_win, kmax_band; color = (:gray, 0.10))

    for s in specs
        li = lines!(axP, s.kxv[s.ok], s.Pk[s.ok], color=s.color, linewidth=3)
        push!(handles, li)
        push!(labels, s.label)

        kp = add_peak_vline!(axP, s.kxv, s.Pk;
                             kmin=kmin_win, kmax=kmax_win,
                             color=s.color, linestyle=:dot, linewidth=4.5)
        @info "P(nu) peak" s.label kp
    end

    axislegend(axP, handles, labels; position=:rt, framevisible=true)

    set_log_safe!(axP; xdata=allx, ydata=ally)
    autolimits!(axP)
    set_pow10_ticks!(axP; which=:both)

    display(figP)
end

# ---- FIGURE 7 : Q^2 spectrum along x, mean over y ----
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figQ2 = Figure(size=(1500, 1200), figure_padding=30)
    axQ2 = Axis(figQ2[1, 1],
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\left\\langle \\left|\\tilde{Q^2}(k_x,y;\\nu_{50})\\right|^2 \\right\\rangle_y"),
    )

    handles = Any[]
    labels  = Any[]

    cols = curve_colors(1 + length(sort(L_ok)))

    Qslice0 = get_chan_xy(Qdata, ichan, n, m)
    Q2img0  = Qslice0 .^ 2
    kxv0, P0 = psd1d_x_mean_over_y(Q2img0, kx; remove_mean=true)
    ok0 = isfinite.(P0) .& (kxv0 .> 0) .& (P0 .> 0)
    l0 = lines!(axQ2, kxv0[ok0], P0[ok0], color=cols[1], linewidth=3)
    push!(handles, l0)
    push!(labels, LaTeXString("\\mathrm{no\\ filter}"))

    allx = Float64[]; ally = Float64[]
    append!(allx, kxv0[ok0]); append!(ally, P0[ok0])

    for (i, Llarge) in enumerate(sort(L_ok))
        haskey(Qslice_filt, Llarge) || continue
        Q2img = Float64.(Qslice_filt[Llarge]).^2

        kxv, Pk = psd1d_x_mean_over_y(Q2img, kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)

        li = lines!(axQ2, kxv[ok], Pk[ok], color=cols[i+1], linewidth=3)
        push!(handles, li)
        push!(labels, LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"))

        append!(allx, kxv[ok]); append!(ally, Pk[ok])
    end

    add_filter_verticals_kx!(axQ2; Llist_pc=Llarge_list, linestyle=:dash)
    axislegend(axQ2, handles, labels; position=:rt, framevisible=true)

    set_log_safe!(axQ2; xdata=allx, ydata=ally)
    autolimits!(axQ2)
    set_pow10_ticks!(axQ2; which=:both)

    display(figQ2)
end

# ============================================================
# FIGURES 8–10 : Qphi, Uphi, Pphi at specific phi-channel
# ============================================================
Qphi_cube = read_FITS(Q_in_phi)
Uphi_cube = read_FITS(U_in_phi)

iphi = 30
ϕ30  = PhiArray[iphi]
@info "Phi channel" iphi ϕ30

Qφ0 = get_chan_xy(Qphi_cube, iphi, n, m)
Uφ0 = get_chan_xy(Uphi_cube, iphi, n, m)
Pφ0 = sqrt.(Qφ0.^2 .+ Uφ0.^2)

# ---- FIGURE 8: Qphi ----
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figQφ = Figure(size=(1500, 1200), figure_padding=30)
    axQφ = Axis(figQφ[1, 1],
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\langle |\\tilde{Q}_{\\phi}(k_x,y)|^2 \\rangle_y"),
    )

    handles = Any[]
    labels  = Any[]

    cols = curve_colors(1 + length(sort(L_ok)))

    kxv0, P0 = psd1d_x_mean_over_y(Qφ0, kx; remove_mean=true)
    ok0 = isfinite.(P0) .& (kxv0 .> 0) .& (P0 .> 0)
    l0 = lines!(axQφ, kxv0[ok0], P0[ok0], color=cols[1], linewidth=3)
    push!(handles, l0)
    push!(labels, LaTeXString("\\mathrm{no\\ filter}"))

    allx = Float64[]; ally = Float64[]
    append!(allx, kxv0[ok0]); append!(ally, P0[ok0])

    for (i, Llarge) in enumerate(sort(L_ok))
        H, _ = instrument_bandpass_L(n, m; Δx=Δx, Δy=Δy, Lcut_small=Lcut_small, Llarge=Llarge, fNy=fNy)
        Qφf = apply_instrument_2d(Qφ0, H)

        kxv, Pk = psd1d_x_mean_over_y(Qφf, kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)

        li = lines!(axQφ, kxv[ok], Pk[ok], color=cols[i+1], linewidth=3)
        push!(handles, li)
        push!(labels, LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"))

        append!(allx, kxv[ok]); append!(ally, Pk[ok])
    end

    add_filter_verticals_kx!(axQφ; Llist_pc=Llarge_list, linestyle=:dash)
    axislegend(axQφ, handles, labels; position=:rt, framevisible=true)

    set_log_safe!(axQφ; xdata=allx, ydata=ally)
    autolimits!(axQφ)
    set_pow10_ticks!(axQφ; which=:both)

    display(figQφ)
end

# ---- FIGURE 9: Uphi ----
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figUφ = Figure(size=(1500, 1200), figure_padding=30)
    axUφ = Axis(figUφ[1, 1],
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\langle |\\tilde{U}_{\\phi}(k_x,y)|^2 \\rangle_y"),
    )

    handles = Any[]
    labels  = Any[]

    cols = curve_colors(1 + length(sort(L_ok)))

    kxv0, P0 = psd1d_x_mean_over_y(Uφ0, kx; remove_mean=true)
    ok0 = isfinite.(P0) .& (kxv0 .> 0) .& (P0 .> 0)
    l0 = lines!(axUφ, kxv0[ok0], P0[ok0], color=cols[1], linewidth=3)
    push!(handles, l0)
    push!(labels, LaTeXString("\\mathrm{no\\ filter}"))

    allx = Float64[]; ally = Float64[]
    append!(allx, kxv0[ok0]); append!(ally, P0[ok0])

    for (i, Llarge) in enumerate(sort(L_ok))
        H, _ = instrument_bandpass_L(n, m; Δx=Δx, Δy=Δy, Lcut_small=Lcut_small, Llarge=Llarge, fNy=fNy)
        Uφf = apply_instrument_2d(Uφ0, H)

        kxv, Pk = psd1d_x_mean_over_y(Uφf, kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)

        li = lines!(axUφ, kxv[ok], Pk[ok], color=cols[i+1], linewidth=3)
        push!(handles, li)
        push!(labels, LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"))

        append!(allx, kxv[ok]); append!(ally, Pk[ok])
    end

    add_filter_verticals_kx!(axUφ; Llist_pc=Llarge_list, linestyle=:dash)
    axislegend(axUφ, handles, labels; position=:rt, framevisible=true)

    set_log_safe!(axUφ; xdata=allx, ydata=ally)
    autolimits!(axUφ)
    set_pow10_ticks!(axUφ; which=:both)

    display(figUφ)
end

# ---- FIGURE 10: Pphi ----
with_theme(theme_latexfonts()) do
    set_theme_spectra!()

    figPφ = Figure(size=(1500, 1200), figure_padding=30)
    axPφ = Axis(figPφ[1, 1],
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("\\langle |\\tilde{P}_{\\phi}(k_x,y)|^2 \\rangle_y"),
    )

    handles = Any[]
    labels  = Any[]

    kmin_win = 0.12
    kmax_win = Inf

    Lsorted = sort(L_ok)
    cols = curve_colors(1 + length(Lsorted))

    specs = Vector{NamedTuple}()

    # no filter
    kxv0, P0 = psd1d_x_mean_over_y(Pφ0, kx; remove_mean=true)
    ok0 = isfinite.(P0) .& (kxv0 .> 0) .& (P0 .> 0)
    push!(specs, (label=LaTeXString("\\mathrm{no\\ filter}"), color=cols[1], kxv=kxv0, Pk=P0, ok=ok0))

    # filtered
    for (i, Llarge) in enumerate(Lsorted)
        H, _ = instrument_bandpass_L(n, m; Δx=Δx, Δy=Δy, Lcut_small=Lcut_small, Llarge=Llarge, fNy=fNy)
        Qφf = apply_instrument_2d(Qφ0, H)
        Uφf = apply_instrument_2d(Uφ0, H)
        Pφf = sqrt.(Qφf.^2 .+ Uφf.^2)

        kxv, Pk = psd1d_x_mean_over_y(Pφf, kx; remove_mean=true)
        ok = isfinite.(Pk) .& (kxv .> 0) .& (Pk .> 0)

        push!(specs, (
            label=LaTeXString("L_{\\mathrm{large}}=$(Llabel_pc(Llarge))\\,\\mathrm{pc}"),
            color=cols[i+1],
            kxv=kxv,
            Pk=Pk,
            ok=ok
        ))
    end

    allx = Float64[]
    ally = Float64[]
    for s in specs
        append!(allx, s.kxv[s.ok])
        append!(ally, s.Pk[s.ok])
    end

    kmax_band = maximum(allx)
    vspan!(axPφ, kmin_win, kmax_band; color = (:gray, 0.10))

    for s in specs
        li = lines!(axPφ, s.kxv[s.ok], s.Pk[s.ok], color=s.color, linewidth=3)
        push!(handles, li)
        push!(labels, s.label)

        kp = add_peak_vline!(axPφ, s.kxv, s.Pk;
                             kmin=kmin_win, kmax=kmax_win,
                             color=s.color, linestyle=:dot, linewidth=4.5)
        @info "P(phi) peak" s.label kp
    end

    axislegend(axPφ, handles, labels; position=:rt, framevisible=true)

    set_log_safe!(axPφ; xdata=allx, ydata=ally)
    autolimits!(axPφ)
    set_pow10_ticks!(axPφ; which=:both)

    # ---------------- INSET: k_peak vs Llarge (phi space) ----------------
    let
        @info "Inset Pphi: k-window" kmin_win kmax_win

        Leff  = Float64[]
        kpeak = Float64[]
        cpeak = Symbol[]

        for (i, Llarge) in enumerate(Lsorted)
            H, _ = instrument_bandpass_L(n, m; Δx=Δx, Δy=Δy, Lcut_small=Lcut_small, Llarge=Llarge, fNy=fNy)

            Qφf = apply_instrument_2d(Qφ0, H)
            Uφf = apply_instrument_2d(Uφ0, H)
            Pφf = sqrt.(Qφf.^2 .+ Uφf.^2)

            kxv, Pk = psd1d_x_mean_over_y(Pφf, kx; remove_mean=true)
            kp = kpeak_in_window(kxv, Pk; kmin=kmin_win, kmax=kmax_win)

            push!(Leff, (50/256) * Llarge)
            push!(kpeak, kp)
            push!(cpeak, cols[i+1])

            @info "Inset Pphi peak" Llarge kp
        end

        ok = isfinite.(Leff) .& isfinite.(kpeak) .& (kpeak .> 0)

        axInset = Axis(figPφ[1,1],
            width  = Relative(0.38),
            height = Relative(0.34),
            halign = :left,
            valign = :bottom,

            xlabel = LaTeXString("L_{\\mathrm{large}}"),
            ylabel = LaTeXString("k_{\\mathrm{x,peak}}"),

            xminorticksvisible = false,
            yminorticksvisible = true,
            yminorticks = IntervalsBetween(9),
            yminorticksize = 5,
            yminortickwidth = 1.5,

            xlabelsize = 25,
            ylabelsize = 25,
            xticksize = 8,
            yticksize = 8,
            xtickwidth = 2,
            ytickwidth = 2,
            xticklabelsize = 22,
            yticklabelsize = 22,

            xtickalign = 1,
            ytickalign = 1,

            xaxisposition = :top,
            yaxisposition = :right,
        )

        lines!(axInset, Leff[ok], kpeak[ok], color=:black, linewidth=2)
        scatter!(axInset, Leff[ok], kpeak[ok]; color=cpeak[ok], markersize=14)

        xlims!(axInset, 0, 40)

        ymin = max(minimum(kpeak[ok]) * 0.9, eps(Float64))
        ymax = maximum(kpeak[ok]) * 1.1
        axInset.yscale = linear
        ylims!(axInset, ymin, ymax)

        axInset.xtickformat = (vs)->[
            LaTeXString("\\mathrm{$(round(v; digits=1))}") for v in vs
        ]
        set_pow10_ticks!(axInset; which=:y)

        axInset.backgroundcolor = (:white, 0.85)
    end

    display(figPφ)
end

# ============================================================
# FIGURE 11 : LIC
# ============================================================
# Bx   = read_FITS(Bx_in)
# By   = read_FITS(By_in)
# Bz   = read_FITS(Bz_in)
# dens = read_FITS(dens_in)

# Bx   = permutedims(Bx,   (1,3,2))
# By   = permutedims(By,   (1,3,2))
# Bz   = permutedims(Bz,   (1,3,2))
# dens = permutedims(dens, (1,3,2))

# den = sum(dens, dims=3)
# den = den .+ (den .== 0)

# meanBx = dropdims(sum(Bx .* dens, dims=3) ./ den, dims=3)
# meanBy = dropdims(sum(By .* dens, dims=3) ./ den, dims=3)
# meanBz = dropdims(sum(Bz .* dens, dims=3) ./ den, dims=3)

# thetay = atan.(meanBx, meanBz)
# xy = @. cos(thetay)
# yy = @. sin(thetay)
# ly = lic(xy, yy, niter=5, len=18, normalize=false)

# with_theme(theme_latexfonts()) do
#     set_theme_heatmaps!()

#     figLIC = Figure(size=(1700, 1350), figure_padding=30)
#     ax = Axis(figLIC[1,1],
#         xlabel = LaTeXString("x\\,[\\mathrm{pc}]"),
#         ylabel = LaTeXString("y\\,[\\mathrm{pc}]"),
#         aspect = DataAspect()
#     )
#     set_latex_linear_ticks_pc!(ax; step_pc=10.0, Lbox_pc=Lbox_pc)

#     hm = heatmap!(ax, x_pc, y_pc, ly; colormap=:grays)
#     Colorbar(figLIC[1,2], hm;
#         label      = LaTeXString("\\mathrm{LIC}\\,[\\mathrm{arb.}]"),
#         tellheight = true,
#         tickformat = cb_latex_ticks
#     )

#     display(figLIC)
# end