using Statistics
using CairoMakie

include(joinpath(@__DIR__, "..", "lib", "DepolLib.jl"))
using .DepolLib

"""
    rms_label(...)

    Extracts RMS token from simulation name.
"""
function rms_label(s::String)
    m = match(r"rms(.*?)nograv1024", s)
    return m === nothing ? s : m.captures[1]
end

"""
    to_cgs(...)

    Converts velocity to CGS.
"""
to_cgs(V) = V .* 1e5
"""
    cs(...)

    Computes sound speed.
"""
cs(T) = sqrt.((5 / 3) * 1.380649e-16 .* T ./ (1.27 * 1.6735575e-24))
"""
    rho(...)

    Computes mass density.
"""
rho(n) = (1.27 * 1.6735575e-24) .* n

"""
    compute_ms_ma(...)

    Computes mean sonic and Alfven Mach numbers.
"""
function compute_ms_ma(root::String, sim::String)
    n = read_fits_f32(simulation_field_path(root, sim, "density"))
    T = read_fits_f32(simulation_field_path(root, sim, "temperature"))

    Bx = read_fits_f32(simulation_field_path(root, sim, "Bx")) .* 1f-3
    By = read_fits_f32(simulation_field_path(root, sim, "By")) .* 1f-3
    Bz = read_fits_f32(simulation_field_path(root, sim, "Bz")) .* 1f-3

    Vx = to_cgs(read_fits_f32(simulation_field_path(root, sim, "Vx")))
    Vy = to_cgs(read_fits_f32(simulation_field_path(root, sim, "Vy")))
    Vz = to_cgs(read_fits_f32(simulation_field_path(root, sim, "Vz")))

    vmag = sqrt.(Vx.^2 .+ Vy.^2 .+ Vz.^2)
    Ms = mean(vmag ./ cs(T))

    Bmag = sqrt.(Bx.^2 .+ By.^2 .+ Bz.^2)
    vA = Bmag ./ sqrt.(4f0 * Float32(pi) .* rho(n))
    MA = mean(vmag ./ vA)

    return Ms, MA
end

"""
    run_mach_suite(...)

    Runs full Mach-suite workflow and writes plots + summary CSV.
"""
function run_mach_suite(cfg)::Dict{String,Any}
    root = resolve_simulations_root(cfg)
    los_y = require_los(string(cfg_get(cfg, ["tasks", "mach_suite", "los_y"]; default="y")))
    los_x = require_los(string(cfg_get(cfg, ["tasks", "mach_suite", "los_x"]; default="x")))

    sims = default_simulation_list()
    required_files = String[]
    for sim in sims
        push!(required_files, withfaraday_path(root, sim, los_y, "Pmax.fits"))
        push!(required_files, withfaraday_path(root, sim, los_x, "Pmax.fits"))
        for field in ("density", "temperature", "Bx", "By", "Bz", "Vx", "Vy", "Vz")
            push!(required_files, simulation_field_path(root, sim, field))
        end
    end
    require_existing_files(required_files; context="mach_suite")

    pth_y_vals = Float32[]
    pth_x_vals = Float32[]
    pmax_finite = Dict{String,Tuple{Vector{Float32},Vector{Float32}}}()
    for sim in sims
        vy = finite_values(read_fits_f32(withfaraday_path(root, sim, los_y, "Pmax.fits")))
        vx = finite_values(read_fits_f32(withfaraday_path(root, sim, los_x, "Pmax.fits")))
        pmax_finite[sim] = (vy, vx)
        append!(pth_y_vals, vy)
        append!(pth_x_vals, vx)
    end
    pth_y = quantile(pth_y_vals, 0.10)
    pth_x = quantile(pth_x_vals, 0.10)

    ms = Float64[]
    ma = Float64[]
    p10_y = Float64[]
    p10_x = Float64[]
    fdep_y = Float64[]
    fdep_x = Float64[]
    rms_tags = String[]

    for sim in sims
        Ms, MA = compute_ms_ma(root, sim)
        vy, vx = pmax_finite[sim]

        push!(ms, Ms)
        push!(ma, MA)
        push!(p10_y, quantile(vy, 0.10))
        push!(p10_x, quantile(vx, 0.10))
        push!(fdep_y, count(<=(pth_y), vy) / length(vy))
        push!(fdep_x, count(<=(pth_x), vx) / length(vx))
        push!(rms_tags, rms_label(sim))
    end

    colors = Dict(
        "09000" => :royalblue,
        "18000" => :darkorange,
        "36000" => :seagreen,
        "72000" => :firebrick,
    )

    p10_plot = standard_output_path(cfg, "mach_suite", "p10_scatter", "pdf"; simu="multi", los=los_y)
    frac_plot = standard_output_path(cfg, "mach_suite", "fraction_scatter", "pdf"; simu="multi", los=los_y)
    summary_csv = standard_output_path(cfg, "mach_suite", "summary", "csv"; simu="multi", los=los_y)

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1100, 850))
        ax1 = Axis(fig[1, 1], xlabel="<Ms>", ylabel="P10 LOS y", xgridvisible=false, ygridvisible=false)
        ax2 = Axis(fig[1, 2], xlabel="<MA>", ylabel="P10 LOS x", xgridvisible=false, ygridvisible=false)

        for (tag, col) in colors
            idx = findall(rms_tags .== tag)
            isempty(idx) && continue
            scatter!(ax1, ms[idx], p10_y[idx], color=col)
            scatter!(ax2, ma[idx], p10_x[idx], color=col, label=tag)
        end
        axislegend(ax2; title="rms")
        save(p10_plot, fig)
    end

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1100, 850))
        ax1 = Axis(fig[1, 1], xlabel="<Ms>", ylabel="f(P<=Pth) LOS y", xgridvisible=false, ygridvisible=false)
        ax2 = Axis(fig[1, 2], xlabel="<MA>", ylabel="f(P<=Pth) LOS x", xgridvisible=false, ygridvisible=false)

        for (tag, col) in colors
            idx = findall(rms_tags .== tag)
            isempty(idx) && continue
            scatter!(ax1, ms[idx], fdep_y[idx], color=col)
            scatter!(ax2, ma[idx], fdep_x[idx], color=col, label=tag)
        end
        axislegend(ax2; title="rms")
        save(frac_plot, fig)
    end

    open(summary_csv, "w") do io
        println(io, "simu,rms,Ms,MA,P10_y,P10_x,fdep_y,fdep_x")
        for i in eachindex(sims)
            println(io, string(
                sims[i], ",", rms_tags[i], ",", ms[i], ",", ma[i], ",",
                p10_y[i], ",", p10_x[i], ",", fdep_y[i], ",", fdep_x[i]
            ))
        end
    end

    return Dict(
        "task" => "mach_suite",
        "n_simulations" => length(sims),
        "thresholds" => Dict("y" => pth_y, "x" => pth_x),
        "outputs" => Dict(
            "p10_plot" => p10_plot,
            "fraction_plot" => frac_plot,
            "summary_csv" => summary_csv,
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_job_entrypoint("mach_suite", run_mach_suite)
end
