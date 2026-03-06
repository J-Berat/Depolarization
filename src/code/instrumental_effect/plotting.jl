cb_latex_ticks(vs) = [LaTeXString("\\mathrm{" * Printf.format(Printf.Format("%.2g"), v) * "}") for v in vs]
pow10_label_from_val(v::Real) = LaTeXString("10^{$(round(Int, log10(v)))}")

function set_latex_linear_ticks_pc!(ax::Axis; step_pc::Real=10.0, Lbox_pc::Real=50.0)
    vals = collect(0.0:step_pc:Lbox_pc)
    ax.xticks = vals
    ax.yticks = vals
    ax.xtickformat = (vs) -> [LaTeXString("\\mathrm{$(round(Int, v))}") for v in vs]
    ax.ytickformat = (vs) -> [LaTeXString("\\mathrm{$(round(Int, v))}") for v in vs]
    return nothing
end

function set_pow10_ticks!(ax::Axis; which::Symbol=:both)
    r = ax.finallimits[]
    xmin = r.origin[1]
    ymin = r.origin[2]
    xmax = r.origin[1] + r.widths[1]
    ymax = r.origin[2] + r.widths[2]

    if (which == :x || which == :both) && xmin > 0 && xmax > 0
        emin = floor(Int, log10(xmin))
        emax = ceil(Int, log10(xmax))
        vals = 10.0 .^ collect(emin:emax)
        ax.xticks = vals
        ax.xtickformat = (vs) -> [pow10_label_from_val(v) for v in vs]
    end

    if (which == :y || which == :both) && ymin > 0 && ymax > 0
        emin = floor(Int, log10(ymin))
        emax = ceil(Int, log10(ymax))
        vals = 10.0 .^ collect(emin:emax)
        ax.yticks = vals
        ax.ytickformat = (vs) -> [pow10_label_from_val(v) for v in vs]
    end
    return nothing
end

function set_log_safe!(ax::Axis; xdata::AbstractVector, ydata::AbstractVector, pad_frac::Real=0.08)
    xok = isfinite.(xdata) .& (xdata .> 0)
    yok = isfinite.(ydata) .& (ydata .> 0)
    @assert any(xok) "No positive finite x data for log axis"
    @assert any(yok) "No positive finite y data for log axis"

    xmin = minimum(xdata[xok])
    xmax = maximum(xdata[xok])
    ymin = minimum(ydata[yok])
    ymax = maximum(ydata[yok])

    xmin = max(xmin * (1 - pad_frac), eps(Float64))
    xmax = xmax * (1 + pad_frac)
    ymin = max(ymin * (1 - pad_frac), eps(Float64))
    ymax = ymax * (1 + pad_frac)

    ax.xscale = log10
    ax.yscale = log10
    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)
    return nothing
end

function add_filter_verticals_kx!(ax; Llist_pc::AbstractVector, linestyle=:dash)
    ks = (2π) ./ Float64.(Llist_pc)
    vlines!(ax, ks; linestyle=linestyle, color=:black)
    return nothing
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
    @assert n <= length(base) "Not enough predefined distinct colors for n=$n"
    return base[1:n]
end

function set_theme_heatmaps!()
    set_theme!(Theme(
        fontsize=44,
        Axis=(
            titlesize=40,
            xlabelsize=56,
            ylabelsize=56,
            xticklabelsize=46,
            yticklabelsize=46,
            spinewidth=2,
            xgridvisible=false,
            ygridvisible=false,
            xminorticksvisible=false,
            yminorticksvisible=false,
        ),
        Colorbar=(
            labelsize=68,
            ticklabelsize=44,
            width=34,
            tickformat=cb_latex_ticks,
        ),
    ))
end

function set_theme_spectra!()
    set_theme!(Theme(
        fontsize=44,
        Axis=(
            titlesize=40,
            xlabelsize=56,
            ylabelsize=56,
            xticklabelsize=46,
            yticklabelsize=46,
            spinewidth=2,
            xgridvisible=false,
            ygridvisible=false,
            xminorticksvisible=true,
            yminorticksvisible=true,
            xminorticks=IntervalsBetween(9),
            yminorticks=IntervalsBetween(9),
            xticksize=14,
            yticksize=14,
            xtickwidth=2,
            ytickwidth=2,
            xminorticksize=9,
            yminorticksize=9,
            xminortickwidth=2,
            yminortickwidth=2,
        ),
        Legend=(labelsize=44,),
        Colorbar=(
            labelsize=62,
            ticklabelsize=44,
            width=34,
            tickformat=cb_latex_ticks,
        ),
    ))
end

"""Build and plot multiple 1D spectra sharing the same layout and styling."""
function plot_multi_spectrum!(ax::Axis, series;
                              add_verticals::Bool=false,
                              Llist_pc::AbstractVector=Float64[],
                              peak_window::Union{Nothing, Tuple{Real, Real}}=nothing)
    handles = Any[]
    labels = Any[]
    allx = Float64[]
    ally = Float64[]

    for s in series
        li = lines!(ax, s.kx[s.ok], s.Pk[s.ok], color=s.color, linewidth=3)
        push!(handles, li)
        push!(labels, s.label)
        append!(allx, s.kx[s.ok])
        append!(ally, s.Pk[s.ok])

        if peak_window !== nothing
            kmin, kmax = peak_window
            add_peak_vline!(ax, s.kx, s.Pk; kmin=kmin, kmax=kmax, color=s.color, linestyle=:dot, linewidth=4.5)
        end
    end

    if add_verticals
        add_filter_verticals_kx!(ax; Llist_pc=Llist_pc, linestyle=:dash)
    end

    axislegend(ax, handles, labels; position=:rt, framevisible=true)
    set_log_safe!(ax; xdata=allx, ydata=ally)
    autolimits!(ax)
    set_pow10_ticks!(ax; which=:both)

    return (allx=allx, ally=ally)
end
