using CairoMakie
using Makie
using Statistics
using FITSIO

# ===============================
# I/O
# ===============================
readfits(path::AbstractString) = Float32.(FITS(path, "r") do f
    read(f[1])
end)

# ===============================
# Finite differences (periodic, same size)
# ===============================
function dd(A::AbstractArray{<:Real,3}, dim::Int, dx::Real)
    A1 = circshift(A, ntuple(d -> d == dim ? -1 : 0, 3))
    A2 = circshift(A, ntuple(d -> d == dim ? +1 : 0, 3))
    return (A1 .- A2) ./ (2f0 * Float32(dx))
end

function div3(vx::AbstractArray{<:Real,3},
              vy::AbstractArray{<:Real,3},
              vz::AbstractArray{<:Real,3},
              dx::Real)
    @assert size(vx) == size(vy) == size(vz) "V components must have identical cube size"
    return dd(vx, 1, dx) .+ dd(vy, 2, dx) .+ dd(vz, 3, dx)
end

function curlmag3(vx::AbstractArray{<:Real,3},
                  vy::AbstractArray{<:Real,3},
                  vz::AbstractArray{<:Real,3},
                  dx::Real)
    @assert size(vx) == size(vy) == size(vz) "V components must have identical cube size"
    dvy_dz = dd(vy, 3, dx); dvz_dy = dd(vz, 2, dx)
    dvz_dx = dd(vz, 1, dx); dvx_dz = dd(vx, 3, dx)
    dvx_dy = dd(vx, 2, dx); dvy_dx = dd(vy, 1, dx)

    cx = dvz_dy .- dvy_dz
    cy = dvx_dz .- dvz_dx
    cz = dvy_dx .- dvx_dy
    return sqrt.(cx .* cx .+ cy .* cy .+ cz .* cz)
end

# ===============================
# LOS reductions
# ===============================
function mean_los(A::AbstractArray{<:Real,3}, los::AbstractString)
    dim = los == "x" ? 1 : los == "y" ? 2 : los == "z" ? 3 : error("los must be \"x\", \"y\" or \"z\"")
    return dropdims(mean(A; dims=dim), dims=dim)
end

function max_los(A::AbstractArray{<:Real,3}, los::AbstractString)
    dim = los == "x" ? 1 : los == "y" ? 2 : los == "z" ? 3 : error("los must be \"x\", \"y\" or \"z\"")
    return dropdims(maximum(A; dims=dim), dims=dim) 
end

# ===============================
# Ticks in pc
# ===============================
function ticks_pc(N::Int; Lbox_pc::Real=50.0, step_pc::Real=10.0)
    pcs = 0:step_pc:Lbox_pc
    pos = round.(Int, 1 .+ (N - 1) .* (pcs ./ Lbox_pc))
    lab = string.(Int.(pcs))
    return pos, lab
end

# ===============================
# Optional colorrange helper
# ===============================
function maybe_colorrange(vmin, vmax)
    (vmin === nothing || vmax === nothing) && return nothing
    return (Float32(vmin), Float32(vmax))
end

# ===============================
# Plot styling
# ===============================
function style!(ax, xt, xl, yt, yl)
    ax.aspect = DataAspect()
    ax.xticks = (xt, xl)
    ax.yticks = (yt, yl)

    ax.xlabel = "Distance [pc]"
    ax.ylabel = "Distance [pc]"

    ax.xlabelsize = 28
    ax.ylabelsize = 28
    ax.xticklabelsize = 24
    ax.yticklabelsize = 24

    ax.titlesize = 28
    ax.titlefont = :regular  # not bold

    ax.xgridvisible = false
    ax.ygridvisible = false
    return ax
end

# ===============================
# Main
# ===============================
function main(;
    los::AbstractString = "y",
    Lbox_pc::Real = 50.0,

    pmax_file::AbstractString,
    vx_file::AbstractString,
    vy_file::AbstractString,
    vz_file::AbstractString,

    canal_q::Real = 0.10,

    # default: automatic (use nothing)
    p_vmin = nothing, p_vmax = nothing,
    w_vmin = nothing, w_vmax = nothing,
    c_vmin = nothing, c_vmax = nothing,

    outprefix::AbstractString = "sim"
)
    # ---- read data
    P = readfits(pmax_file)
    vx = readfits(vx_file)
    vy = readfits(vy_file)
    vz = readfits(vz_file)

    @assert ndims(P) == 2 "Pmax must be a 2D map"
    @assert ndims(vx) == 3 && ndims(vy) == 3 && ndims(vz) == 3 "V must be 3D cubes"
    @assert size(vx) == size(vy) == size(vz) "V cubes must have identical size"

    # ---- geometry
    N = size(vx, 1)
    dx_pc = Lbox_pc / N  # pc per cell (consistent with your 256 pixels over 50 pc convention)

    # ---- derived cubes
    W = curlmag3(vx, vy, vz, dx_pc)          # |∇×v|  [ (km/s)/pc ] if v in km/s
    C = -div3(vx, vy, vz, dx_pc)             # -∇·v   [ (km/s)/pc ]

    # ---- LOS maps
    Wμ = mean_los(W, los);  Wx = max_los(W, los)
    Cμ = mean_los(C, los);  Cx = max_los(C, los)

    # ---- canal mask from Pmax (same 2D plane as the plots)
    Pthr = quantile(vec(P), canal_q)
    canal = Float32.(P .< Pthr)  # 0/1 map for contour

    # ---- ticks
    xt, xl = ticks_pc(size(P,1); Lbox_pc=Lbox_pc, step_pc=10.0)
    yt, yl = ticks_pc(size(P,2); Lbox_pc=Lbox_pc, step_pc=10.0)

    # ---- labels (MathTeXEngine-safe: use \mathrm)
    LBLP = L"P_{\mathrm{max}}\,[\mathrm{K}]"
    LBLW = L"|\nabla\times v|\,[\mathrm{km\,s^{-1}\,pc^{-1}}]"
    LBLC = L"-\nabla\cdot v\,[\mathrm{km\,s^{-1}\,pc^{-1}}]"

    # ---- colorranges (only pass if both bounds are set)
    crP = maybe_colorrange(p_vmin, p_vmax)
    crW = maybe_colorrange(w_vmin, w_vmax)
    crC = maybe_colorrange(c_vmin, c_vmax)

    # helper to call heatmap! with optional colorrange
    function hm!(ax, A; cmap, cr)
        if cr === nothing
            return heatmap!(ax, A; colormap=cmap)
        else
            return heatmap!(ax, A; colormap=cmap, colorrange=cr)
        end
    end

    with_theme(theme_latexfonts()) do
        # =========================
        # FIG 1 : MEAN
        # =========================
        f1 = Figure(size=(1200, 800))

        ax11 = style!(Axis(f1[1,1]), xt,xl, yt,yl)
        ax12 = style!(Axis(f1[1,3]), xt,xl, yt,yl)
        ax22 = style!(Axis(f1[2,3]), xt,xl, yt,yl)

        hP = hm!(ax11, P;  cmap=:magma,   cr=crP)
        hW = hm!(ax12, Wμ; cmap=:viridis, cr=crW)
        hC = hm!(ax22, Cμ; cmap=:viridis, cr=crC)

        ax11.title = L"P_{\mathrm{max}}"
        ax12.title = L"\mathrm{mean}(|\nabla\times v|)"
        ax22.title = L"\mathrm{mean}(−\nabla\cdot v)"

        # canals contours over W and C
        contour!(ax12, canal; levels=[0.5], color=:white, linewidth=2)
        contour!(ax22, canal; levels=[0.5], color=:white, linewidth=2)

        Colorbar(f1[1,2], hP, label=LBLP, labelsize=26, ticklabelsize=22)
        Colorbar(f1[1,4], hW, label=LBLW, labelsize=26, ticklabelsize=22)
        Colorbar(f1[2,4], hC, label=LBLC, labelsize=26, ticklabelsize=22)

        save(outprefix * "_mean.pdf", f1)
        display(f1)

        # =========================
        # FIG 2 : MAX
        # =========================
        f2 = Figure(size=(1100, 850))

        bx11 = style!(Axis(f2[1,1]), xt,xl, yt,yl)
        bx12 = style!(Axis(f2[1,3]), xt,xl, yt,yl)
        bx22 = style!(Axis(f2[2,3]), xt,xl, yt,yl)

        gP = hm!(bx11, P;  cmap=:magma,   cr=crP)
        gW = hm!(bx12, Wx; cmap=:viridis, cr=crW)
        gC = hm!(bx22, Cx; cmap=:viridis, cr=crC)

        bx11.title = L"P_{\mathrm{max}}"
        bx12.title = L"\mathrm{max}(|\nabla\times v|)"
        bx22.title = L"\mathrm{max}(−\nabla\cdot v)"

        contour!(bx12, canal; levels=[0.5], color=:white, linewidth=2)
        contour!(bx22, canal; levels=[0.5], color=:white, linewidth=2)

        Colorbar(f2[1,2], gP, label=LBLP, labelsize=26, ticklabelsize=22)
        Colorbar(f2[1,4], gW, label=LBLW, labelsize=26, ticklabelsize=22)
        Colorbar(f2[2,4], gC, label=LBLC, labelsize=26, ticklabelsize=22)

        save(outprefix * "_max.pdf", f2)
        display(f2)

        return f1, f2
    end
end

# ===============================
# Example call
# ===============================
main(;
    los="y",
    Lbox_pc=50.0,

    pmax_file="/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Pmax.fits",

    # NOTE: no "y/physical_param" here
    vx_file="/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/Vx.fits",
    vy_file="/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/Vy.fits",
    vz_file="/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/Vz.fits",

    canal_q=0.10,

    # default: automatic colorranges
    p_vmin=nothing, p_vmax=nothing,
    w_vmin=nothing, w_vmax=nothing,
    c_vmin=nothing, c_vmax=nothing,

    outprefix="sim7"
)
