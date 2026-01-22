using FFTW;
using FITSIO;
using CairoMakie;
using LaTeXStrings;
using Statistics;
using Printf;
include(joinpath(@__DIR__, "../io/fits_io.jl"));

# ============================================================
# PATHS
# ============================================================
Q_in =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Qnu.fits";

U_in =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Unu.fits";

Pmax_nofilter_path =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Pmax.fits";

base_out =
"/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/WithInstrument";

mkpath(base_out);

# ============================================================
# GRID / PHYSICAL SCALE
# ============================================================
n, m = 256, 256;
Lbox_pc = 50.0;
Δx = Lbox_pc / n;
Δy = Lbox_pc / m;

# fixed small-scale removal: remove [0, 1] pc
Lcut_small = 1.0;

# large-scale removal variants: remove [Llarge, 50] pc
Llarge_list = [49, 45, 40];

# Nyquist in cycles/pc
fNy_x = 1 / (2Δx);
fNy_y = 1 / (2Δy);
fNy   = min(fNy_x, fNy_y);

@info "Grid" n m Δx Δy fNy

# ============================================================
# RM-SYNTHESIS AXES
# ============================================================
νmin_MHz = 120.0;
νmax_MHz = 167.0;
Δν_MHz   = 0.098;

nuArray  = (νmin_MHz:Δν_MHz:νmax_MHz) .* 1e6;
PhiArray = collect(-10.0:0.25:10.0);

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

Pphi_max_map(absF::AbstractArray) =
    dropdims(maximum(absF; dims=ndims(absF)); dims=ndims(absF))

# ============================================================
# HARD BAND-PASS MASK (cycles/pc) WITH NYQUIST CLAMP
# keep flo <= f <= min(1/Lcut_small, fNy)
# flo = 1/Llarge  (remove large scales L >= Llarge)
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

    H = Float64.((F .>= flo) .& (F .<= fhi))
    return H, fftshift(H)
end

# ============================================================
# FILTER Q/U → RM-SYNTHESIS → Pmax
# ============================================================
Qdata = read_FITS(Q_in);
Udata = read_FITS(U_in);

for Llarge in Llarge_list
    tag = @sprintf("HardBandPass_remove_L0_to_1pc_and_%dto50pc", Int(round(Llarge)));
    subdir = joinpath(base_out, tag);
    mkpath(subdir);

    H, _ = instrument_bandpass_L(n, m;
                                 Δx=Δx, Δy=Δy,
                                 Lcut_small=Lcut_small,
                                 Llarge=Llarge,
                                 fNy=fNy);

    Qf = apply_to_array_xy(Qdata, H; n=n, m=m);
    Uf = apply_to_array_xy(Udata, H; n=n, m=m);

    write_FITS(joinpath(subdir, "Qnu_filtered.fits"), Qf);
    write_FITS(joinpath(subdir, "Unu_filtered.fits"), Uf);

    absF, _, _ = RMSynthesis(Qf, Uf, nuArray, PhiArray; log_progress=true);

    Pmaxf = Pphi_max_map(absF);
    mkpath(joinpath(subdir, "RMSynthesis"));
    write_FITS(joinpath(subdir, "RMSynthesis", "Pphi_max.fits"), Pmaxf);
end

# ============================================================
# LOAD RESULTS FOR PLOT
# ============================================================
Pmax0 = Float64.(read_FITS(Pmax_nofilter_path))

Pmax_filt = Dict{Float64, Matrix{Float64}}()
Hshift    = Dict{Float64, Matrix{Float64}}()
L_ok      = Float64[]

for Llarge in Llarge_list
    tag = @sprintf("HardBandPass_remove_L0_to_1pc_and_%dto50pc", Int(round(Llarge)))
    p = joinpath(base_out, tag, "RMSynthesis", "Pphi_max.fits")
    isfile(p) || continue

    Pmax_filt[Llarge] = Float64.(read_FITS(p))

    _, Hs = instrument_bandpass_L(n, m;
                                  Δx=Δx, Δy=Δy,
                                  Lcut_small=Lcut_small,
                                  Llarge=Llarge,
                                  fNy=fNy)
    Hshift[Llarge] = Hs
    push!(L_ok, Llarge)
end

@assert !isempty(L_ok) "No filtered Pphi_max.fits found under: $base_out"

# ============================================================
# COLOR RANGES
# ============================================================
vals = Float64[]
append!(vals, vec(Pmax0))
for Llarge in L_ok
    append!(vals, vec(Pmax_filt[Llarge]))
end
vals = vals[isfinite.(vals)]
sort!(vals)
lo = vals[clamp(round(Int, 0.01length(vals)), 1, length(vals))]
hi = vals[clamp(round(Int, 0.99length(vals)), 1, length(vals))]
pmax_colorrange = (lo, hi)

inst_colorrange = (0.0, 1.0)

# ============================================================
# AXES
# ============================================================
kx = (2π) .* fftshift(fftfreq(n, Δx))
ky = (2π) .* fftshift(fftfreq(m, Δy))
x_pc = (0:n-1) .* Δx
y_pc = (0:m-1) .* Δy

# ============================================================
# PLOT (big labels + colorbars)
# ============================================================
with_theme(theme_latexfonts()) do
    set_theme!(Theme(
        fontsize = 36,
        Axis = (
            titlesize      = 30,
            xlabelsize     = 40,
            ylabelsize     = 40,
            xticklabelsize = 32,
            yticklabelsize = 32,
            spinewidth     = 2,
        ),
        Colorbar = (
            labelsize     = 40,
            ticklabelsize = 30,
            width         = 60,
        ),
    ))

    fig = Figure(size=(3000, 1500), figure_padding=30)

    # --- Instruments
    axA = Axis(fig[1,1],
        title  = LaTeXString("Instrument\\ (all\\ pass)"),
        xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
        ylabel = LaTeXString("k_y\\,[\\mathrm{rad\\,pc^{-1}}]"),
        aspect = DataAspect()
    )
    hmA = heatmap!(axA, kx, ky, ones(n,m); colorrange=inst_colorrange)

    inst_ref = hmA

    for (j, Llarge) in enumerate(L_ok)
        ax = Axis(fig[1, j+1],
            title  = LaTeXString("remove\\ L\\in[0,1]\\cup[$(Int(round(Llarge))),50]\\,pc"),
            xlabel = LaTeXString("k_x\\,[\\mathrm{rad\\,pc^{-1}}]"),
            aspect = DataAspect(),
            ylabelvisible = false,
            yticklabelsvisible = false
        )
        inst_ref = heatmap!(ax, kx, ky, Hshift[Llarge]; colorrange=inst_colorrange)
    end

    Colorbar(fig[1, length(L_ok)+2], inst_ref,
             label=LaTeXString("|H(\\mathbf{k})|"))

    # --- Pmax maps
    ax0 = Axis(fig[2,1],
        title  = LaTeXString("P_{\\max}"),
        xlabel = LaTeXString("x\\,[pc]"),
        ylabel = LaTeXString("y\\,[pc]"),
        aspect = DataAspect()
    )
    hm0 = heatmap!(ax0, x_pc, y_pc, Pmax0; colorrange=pmax_colorrange)

    for (j, Llarge) in enumerate(L_ok)
        ax = Axis(fig[2, j+1],
            title  = LaTeXString("P_{\\max}^{f}"),
            xlabel = LaTeXString("x\\,[pc]"),
            aspect = DataAspect(),
            ylabelvisible = false,
            yticklabelsvisible = false
        )
        heatmap!(ax, x_pc, y_pc, Pmax_filt[Llarge]; colorrange=pmax_colorrange)
    end

    Colorbar(fig[2, length(L_ok)+2], hm0,
             label=LaTeXString("P_{\\max}"))

    fig
end
