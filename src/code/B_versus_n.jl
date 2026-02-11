using FITSIO
using CairoMakie
using Statistics
using LaTeXStrings

function plot_B_vs_n(; 
    data_dir::AbstractString = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024",
    file_bx::AbstractString = "Bx.fits",
    file_by::AbstractString = "By.fits",
    file_bz::AbstractString = "Bz.fits",
    file_n::AbstractString  = "density.fits",
    nbins::Int = 30,
    min_per_bin::Int = 50,
    out_fig::AbstractString = "B_vs_n.pdf"
)

    # -------------------------
    # LOAD DATA
    # -------------------------
    Bx = read(FITS(joinpath(data_dir, file_bx))[1])
    By = read(FITS(joinpath(data_dir, file_by))[1])
    Bz = read(FITS(joinpath(data_dir, file_bz))[1])
    n  = read(FITS(joinpath(data_dir, file_n ))[1])

    # -------------------------
    # |B| AND UNIT CONVERSION
    # -------------------------
    B = sqrt.(Bx.^2 .+ By.^2 .+ Bz.^2)
    B .= 1000 .* B  # µG

    nv = vec(n)
    Bv = vec(B)

    mask = isfinite.(nv) .& isfinite.(Bv) .& (nv .> 0) .& (Bv .> 0)
    nv = nv[mask]
    Bv = Bv[mask]

    # -------------------------
    # BINNING
    # -------------------------
    function binned_stats(nv, Bv; nbins::Int=30, min_per_bin::Int=50)
        edges = range(minimum(nv), maximum(nv), length=nbins+1)

        x_center = Float64[]
        B_med    = Float64[]
        B_wmean  = Float64[]
        B_rms    = Float64[]

        for i in 1:nbins
            sel = (nv .>= edges[i]) .& (nv .< edges[i+1])
            count(sel) < min_per_bin && continue

            nbin = nv[sel]
            Bbin = Bv[sel]

            xc = 0.5 * (edges[i] + edges[i+1])

            push!(x_center, xc)
            push!(B_med, median(Bbin))
            push!(B_wmean, sum(Bbin .* nbin) / sum(nbin))
            push!(B_rms, sqrt(mean(Bbin .^ 2)))
        end

        x_center, B_med, B_wmean, B_rms
    end

    x_all, B_med, B_wmean, B_rms = binned_stats(nv, Bv; nbins=nbins, min_per_bin=min_per_bin)

    # -------------------------
    # LABELS
    # -------------------------
    xlab = L"n\;[\mathrm{cm^{-3}}]"
    ylab = L"\text{Intensity}\;[\mu\mathrm{G}]"

    # -------------------------
    # PLOT
    # -------------------------
    fig = nothing
    with_theme(theme_latexfonts()) do
        fig = Figure(size=(820, 560))

        ax = Axis(
            fig[1, 1],
            xlabel = xlab,
            ylabel = ylab,
            limits = (0, 200, nothing, nothing),
            xticks = 0:50:200,
            xgridvisible = false,
            ygridvisible = false,
            xlabelsize = 22,
            ylabelsize = 22,
            xticklabelsize = 18,
            yticklabelsize = 18
        )

        lines!(ax, x_all, B_med;
            linewidth = 2.5,
            color = :royalblue,
            label = L"\mathrm{median}(B)"
        )

        lines!(ax, x_all, B_wmean;
            linewidth = 2.5,
            color = :darkorange,
            label = L"\langle B \rangle_{n}"
        )

        lines!(ax, x_all, B_rms;
            linewidth = 2.5,
            color = :seagreen,
            label = L"B_{\mathrm{rms}}"
        )

        vlines!(ax, [100.0];
            linestyle = :dash,
            linewidth = 2,
            color = :black
        )

        axislegend(ax; position=:rb, framevisible=false, labelsize=16)

        save(joinpath(data_dir, out_fig), fig)
    end

    fig
end

plot_B_vs_n()
