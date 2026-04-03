"""
    cb_latex_ticks(...)

    Formats colorbar ticks as LaTeX strings.
"""
cb_latex_ticks(vs) = [LaTeXString("\\mathrm{" * Printf.format(Printf.Format("%.2g"), v) * "}") for v in vs]
const UNIFORM_FIG_TICKLABEL_SIZE = 34
const UNIFORM_FIG_COLORBAR_LABEL_SIZE = 34
const UNIFORM_FIG_COLORBAR_TICKLABEL_SIZE = 34
"""
    pow10_label_from_val(...)

    Formats values as `10^n` labels.
"""
pow10_label_from_val(v::Real) = LaTeXString("10^{$(round(Int, log10(v)))}")

"""
    latex_linear_tickformat(...)

    Formats linear ticks as LaTeX strings.
"""
latex_linear_tickformat(vs; digits::Int=2) = [LaTeXString("\\mathrm{$(round(v; digits=digits))}") for v in vs]

"""
    set_latex_linear_ticks_pc!(...)

    Configures linear parsec ticks in LaTeX format.
"""
function set_latex_linear_ticks_pc!(ax::Axis; step_pc::Real=10.0, Lbox_pc::Real=50.0)
    vals = collect(0.0:step_pc:Lbox_pc)
    ax.xticks = vals
    ax.yticks = vals
    ax.xtickformat = (vs) -> [LaTeXString("\\mathrm{$(round(Int, v))}") for v in vs]
    ax.ytickformat = (vs) -> [LaTeXString("\\mathrm{$(round(Int, v))}") for v in vs]
    return nothing
end

"""
    set_pow10_ticks!(...)

    Configures power-of-10 ticks on log axes.
"""
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
        ax.xminorticksvisible = true
        ax.xminorticks = IntervalsBetween(9)
    end

    if (which == :y || which == :both) && ymin > 0 && ymax > 0
        emin = floor(Int, log10(ymin))
        emax = ceil(Int, log10(ymax))
        vals = 10.0 .^ collect(emin:emax)
        ax.yticks = vals
        ax.ytickformat = (vs) -> [pow10_label_from_val(v) for v in vs]
        ax.yminorticksvisible = true
        ax.yminorticks = IntervalsBetween(9)
    end
    return nothing
end

"""
    set_log_safe!(...)

    Safely enables log scaling with positive finite limits.
"""
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

"""
    set_linear_safe!(...)

    Safely enables linear scaling with finite limits.
"""
function set_linear_safe!(ax::Axis; xdata::AbstractVector, ydata::AbstractVector, pad_frac::Real=0.08)
    xok = isfinite.(xdata)
    yok = isfinite.(ydata)
    @assert any(xok) "No finite x data for linear axis"
    @assert any(yok) "No finite y data for linear axis"

    xmin = minimum(xdata[xok])
    xmax = maximum(xdata[xok])
    ymin = minimum(ydata[yok])
    ymax = maximum(ydata[yok])

    dx = max(xmax - xmin, eps(Float64))
    dy = max(ymax - ymin, eps(Float64))

    ax.xscale = identity
    ax.yscale = identity
    xlims!(ax, xmin - pad_frac * dx, xmax + pad_frac * dx)
    ylims!(ax, ymin - pad_frac * dy, ymax + pad_frac * dy)
    return nothing
end

"""
    configure_axis_style!(...)

    Unified axis styling for spectra and insets.
"""
function configure_axis_style!(ax::Axis; xdata::AbstractVector, ydata::AbstractVector,
                               xscale::Symbol=:log10, yscale::Symbol=:log10,
                               xminor_subdiv::Int=9, yminor_subdiv::Int=9,
                               xaxisposition::Symbol=:bottom, yaxisposition::Symbol=:left,
                               x_linear_digits::Int=2, y_linear_digits::Int=2)
    ax.xaxisposition = xaxisposition
    ax.yaxisposition = yaxisposition

    ax.xminorticksvisible = xminor_subdiv > 0
    ax.yminorticksvisible = yminor_subdiv > 0
    if xminor_subdiv > 0
        ax.xminorticks = IntervalsBetween(xminor_subdiv)
    end
    if yminor_subdiv > 0
        ax.yminorticks = IntervalsBetween(yminor_subdiv)
    end

    if xscale == :log10 || yscale == :log10
        set_log_safe!(ax; xdata=xdata, ydata=ydata)
        autolimits!(ax)
        if xscale == :log10 && yscale == :log10
            set_pow10_ticks!(ax; which=:both)
        elseif xscale == :log10
            set_pow10_ticks!(ax; which=:x)
            ax.ytickformat = (vs) -> latex_linear_tickformat(vs; digits=y_linear_digits)
        else
            set_pow10_ticks!(ax; which=:y)
            ax.xtickformat = (vs) -> latex_linear_tickformat(vs; digits=x_linear_digits)
        end
    else
        set_linear_safe!(ax; xdata=xdata, ydata=ydata)
        autolimits!(ax)
        ax.xtickformat = (vs) -> latex_linear_tickformat(vs; digits=x_linear_digits)
        ax.ytickformat = (vs) -> latex_linear_tickformat(vs; digits=y_linear_digits)
    end
    return nothing
end

"""
    add_filter_verticals_kx!(...)

    Adds vertical lines for filter scales.
"""
function add_filter_verticals_kx!(ax; Llist_pc::AbstractVector, linestyle=:dash)
    ks = 1.0 ./ Float64.(Llist_pc)
    vlines!(ax, ks; linestyle=linestyle, color=:black)
    return nothing
end

"""
    add_peak_vline!(...)

    Adds vertical line at detected peak.
"""
function add_peak_vline!(ax, k::AbstractVector, y::AbstractVector;
                         kmin::Real=0.0, kmax::Real=Inf,
                         color=:black, linestyle=:dot, linewidth=4.5,
                         smooth_radius::Int=2)
    kp = kpeak_in_window(k, y; kmin=kmin, kmax=kmax, smooth_radius=smooth_radius)
    if isfinite(kp)
        vlines!(ax, [kp]; color=color, linestyle=linestyle, linewidth=linewidth)
    end
    return kp
end

"""
    curve_colors(...)

    Returns discrete plotting color palette.
"""
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

"""
    set_theme_heatmaps!(...)

    Applies heatmap plotting theme.
"""
function set_theme_heatmaps!()
    set_theme!(Theme(
        fontsize=44,
        Axis=(
            titlesize=40,
            xlabelsize=56,
            ylabelsize=56,
            xticklabelsize=UNIFORM_FIG_TICKLABEL_SIZE,
            yticklabelsize=UNIFORM_FIG_TICKLABEL_SIZE,
            spinewidth=2,
            xgridvisible=false,
            ygridvisible=false,
            xminorticksvisible=false,
            yminorticksvisible=false,
        ),
        Colorbar=(
            labelsize=UNIFORM_FIG_COLORBAR_LABEL_SIZE,
            ticklabelsize=UNIFORM_FIG_COLORBAR_TICKLABEL_SIZE,
            width=34,
            tickformat=cb_latex_ticks,
        ),
    ))
end

"""
    set_theme_spectra!(...)

    Applies spectra plotting theme.
"""
function set_theme_spectra!()
    set_theme!(Theme(
        fontsize=44,
        Axis=(
            titlesize=40,
            xlabelsize=56,
            ylabelsize=56,
            xticklabelsize=UNIFORM_FIG_TICKLABEL_SIZE,
            yticklabelsize=UNIFORM_FIG_TICKLABEL_SIZE,
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
        Legend=(labelsize=34,),
        Colorbar=(
            labelsize=UNIFORM_FIG_COLORBAR_LABEL_SIZE,
            ticklabelsize=UNIFORM_FIG_COLORBAR_TICKLABEL_SIZE,
            width=34,
            tickformat=cb_latex_ticks,
        ),
    ))
end

"""
    plot_multi_spectrum!(...)

    Plots multiple spectra with legend and optional peak/vertical annotations.
"""
function _smooth_positive_series(y::AbstractVector; radius::Int=2)
    n = length(y)
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        acc = 0.0
        w = 0
        i1 = max(1, i - radius)
        i2 = min(n, i + radius)
        for j in i1:i2
            yj = y[j]
            if isfinite(yj) && yj > 0
                acc += log(Float64(yj))
                w += 1
            end
        end
        out[i] = w > 0 ? exp(acc / w) : NaN
    end
    return out
end

function plot_multi_spectrum!(ax::Axis, series;
                              add_verticals::Bool=false,
                              Llist_pc::AbstractVector=Float64[],
                              peak_window::Union{Nothing, Tuple{Real, Real}}=nothing)
    handles = Any[]
    labels = Any[]
    allx = Float64[]
    ally = Float64[]

    for s in series
        append!(allx, s.kx[s.ok])
        append!(ally, s.Pk[s.ok])
    end

    if peak_window !== nothing && !isempty(allx)
        kmin, kmax = peak_window
        xhi = isfinite(kmax) ? float(kmax) : maximum(allx)
        xlo = float(kmin)
        if isfinite(xlo) && isfinite(xhi) && xhi > xlo
            # Highlight peak-search window.
            vspan!(ax, xlo, xhi; color=(:lightgray, 0.30))
        end
    end

    for s in series
        Pplot = _smooth_positive_series(s.Pk; radius=2)
        if hasproperty(s, :Pstd)
            Pstd_plot = _smooth_positive_series(max.(s.Pstd, eps(Float64)); radius=2)
            ylo = max.(Pplot .- Pstd_plot, 0.25 .* Pplot, eps(Float64))
            yhi = Pplot .+ Pstd_plot
            okband = s.ok .& isfinite.(Pplot) .& isfinite.(ylo) .& isfinite.(yhi) .& (Pplot .> 0) .& (ylo .> 0) .& (yhi .> 0)
            band!(ax, s.kx[okband], ylo[okband], yhi[okband]; color=(s.color, 0.10))
        end
        okline = s.ok .& isfinite.(Pplot) .& (Pplot .> 0)
        li = lines!(ax, s.kx[okline], Pplot[okline], color=s.color, linewidth=3)
        push!(handles, li)
        push!(labels, s.label)

        if peak_window !== nothing
            kmin, kmax = peak_window
            add_peak_vline!(ax, s.kx, s.Pk; kmin=kmin, kmax=kmax, color=s.color, linestyle=:dot, linewidth=4.5)
        end
    end

    if add_verticals
        add_filter_verticals_kx!(ax; Llist_pc=Llist_pc, linestyle=:dash)
    end

    axislegend(ax, handles, labels; position=:rt, framevisible=true)
    configure_axis_style!(ax;
        xdata=allx,
        ydata=ally,
        xscale=:log10,
        yscale=:log10,
        xminor_subdiv=9,
        yminor_subdiv=9,
        xaxisposition=:bottom,
        yaxisposition=:left,
    )

    return (allx=allx, ally=ally)
end
