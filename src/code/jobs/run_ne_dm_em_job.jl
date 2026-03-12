using CairoMakie
using LaTeXStrings

include(joinpath(@__DIR__, "..", "lib", "DepolLib.jl"))
using .DepolLib

"""
    _safe_colorrange(...)

Returns a valid colorrange even for uniform maps.
"""
function _safe_colorrange(m::AbstractMatrix)
    vals = DepolLib.finite_values(m)
    isempty(vals) && return nothing
    lo, hi = extrema(vals)
    if lo == hi
        δ = lo == 0 ? 1.0 : abs(lo) * 1e-3
        return (lo - δ, hi + δ)
    end
    return (lo, hi)
end

"""
    _pc_ticks(...)

Builds parsec ticks from 0 to Lbox with a fixed step.
"""
function _pc_ticks(lbox_pc::Real; step_pc::Real=10.0)
    lbox = Float64(lbox_pc)
    step = Float64(step_pc)
    step > 0 || error("tick step must be > 0, got $step")
    ticks = collect(0.0:step:lbox)
    if isempty(ticks) || ticks[end] < lbox
        push!(ticks, lbox)
    end
    return ticks
end

"""
    _plot_dm_em_maps(...)

Writes a side-by-side DM/EM map figure.
"""
function _plot_dm_em_maps(dm::AbstractMatrix, em::AbstractMatrix, out_plot::AbstractString;
                          los::AbstractString, lbox_pc::Real=50.0,
                          title_size::Int=44, label_size::Int=36, tick_size::Int=30,
                          cbar_label_size::Int=44, cbar_tick_size::Int=36, cbar_width::Int=60,
                          tick_step_pc::Real=10.0)
    size(dm) == size(em) || error("DM/EM map size mismatch: dm=$(size(dm)) em=$(size(em))")
    xlab, ylab = DepolLib.sky_plane_labels(los)

    x = DepolLib.axis_pc(size(dm, 2); lbox_pc=lbox_pc)
    y = DepolLib.axis_pc(size(dm, 1); lbox_pc=lbox_pc)
    ticks_pc = _pc_ticks(lbox_pc; step_pc=tick_step_pc)

    cr_dm = _safe_colorrange(dm)
    cr_em = _safe_colorrange(em)

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1900, 820))

        ax_dm = Axis(fig[1, 1],
            title=L"\mathrm{Dispersion\ Measure}\ (DM)",
            xlabel=string(xlab, " [pc]"),
            ylabel=string(ylab, " [pc]"),
            titlesize=title_size,
            xlabelsize=label_size,
            ylabelsize=label_size,
            xticklabelsize=tick_size,
            yticklabelsize=tick_size,
            xticks=ticks_pc,
            yticks=ticks_pc,
            xgridvisible=false, ygridvisible=false)
        hm_dm = heatmap!(ax_dm, x, y, dm; colormap=:viridis, colorrange=cr_dm)
        Colorbar(fig[1, 2], hm_dm;
            label=L"DM\ [\mathrm{pc}\ \mathrm{cm}^{-3}]",
            labelsize=cbar_label_size,
            ticklabelsize=cbar_tick_size,
            width=cbar_width)

        ax_em = Axis(fig[1, 3],
            title=L"\mathrm{Emission\ Measure}\ (EM)",
            xlabel=string(xlab, " [pc]"),
            ylabel=string(ylab, " [pc]"),
            titlesize=title_size,
            xlabelsize=label_size,
            ylabelsize=label_size,
            xticklabelsize=tick_size,
            yticklabelsize=tick_size,
            xticks=ticks_pc,
            yticks=ticks_pc,
            xgridvisible=false, ygridvisible=false)
        hm_em = heatmap!(ax_em, x, y, em; colormap=:magma, colorrange=cr_em)
        Colorbar(fig[1, 4], hm_em;
            label=L"EM\ [\mathrm{pc}\ \mathrm{cm}^{-6}]",
            labelsize=cbar_label_size,
            ticklabelsize=cbar_tick_size,
            width=cbar_width)

        save(out_plot, fig)
    end

    return out_plot
end

"""
    _plot_measure_heatmap(...)

Writes one DM or EM heatmap figure.
"""
function _plot_measure_heatmap(m::AbstractMatrix, out_plot::AbstractString;
                               los::AbstractString, lbox_pc::Real,
                               title::LaTeXString, cbar_label::LaTeXString,
                               colormap=:viridis,
                               title_size::Int=44, label_size::Int=36, tick_size::Int=30,
                               cbar_label_size::Int=44, cbar_tick_size::Int=36, cbar_width::Int=60,
                               tick_step_pc::Real=10.0)
    xlab, ylab = DepolLib.sky_plane_labels(los)
    x = DepolLib.axis_pc(size(m, 2); lbox_pc=lbox_pc)
    y = DepolLib.axis_pc(size(m, 1); lbox_pc=lbox_pc)
    ticks_pc = _pc_ticks(lbox_pc; step_pc=tick_step_pc)

    cr = _safe_colorrange(m)

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1000, 820))
        ax = Axis(fig[1, 1],
            title=title,
            xlabel=string(xlab, " [pc]"),
            ylabel=string(ylab, " [pc]"),
            titlesize=title_size,
            xlabelsize=label_size,
            ylabelsize=label_size,
            xticklabelsize=tick_size,
            yticklabelsize=tick_size,
            xticks=ticks_pc,
            yticks=ticks_pc,
            xgridvisible=false, ygridvisible=false)
        hm = heatmap!(ax, x, y, m; colormap=colormap, colorrange=cr)
        Colorbar(fig[1, 2], hm;
            label=cbar_label,
            labelsize=cbar_label_size,
            ticklabelsize=cbar_tick_size,
            width=cbar_width)
        save(out_plot, fig)
    end

    return out_plot
end

"""
    ask_user(...)

Prompts for a numeric value and returns parsed `Float64`, or `default` on empty input.
"""
function ask_user(label::AbstractString, default::Float64)
    while true
        print(" - ", label, " [default=", default, "]: ")
        flush(stdout)

        input = try
            readline(stdin)
        catch
            return default
        end
        s = strip(input)
        isempty(s) && return default

        try
            return parse(Float64, s)
        catch
            println("Invalid value '", s, "'. Please enter a number.")
        end
    end
end

"""
    resolve_wolfire_constants(...)

Gets Wolfire constants from config defaults, optionally prompting user in TTY mode.
"""
function resolve_wolfire_constants(cfg::AbstractDict)
    zeta_default = Float64(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "zeta"]; default=2.5e-16))
    geff_default = Float64(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "Geff"]; default=1.0))
    phiPAH_default = Float64(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "phiPAH"]; default=0.5))
    xc_default = Float64(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "XC"]; default=1.4e-4))

    prompt_constants = Bool(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "prompt_constants"]; default=true))
    if prompt_constants
        if Base.stdin isa Base.TTY
            println("Please enter the values for the constants:")
            zeta = ask_user("zeta (ionization rate by Cosmic Rays)", zeta_default)
            geff = ask_user("Geff (effective radiation field)", geff_default)
            phiPAH = ask_user("phiPAH (collision rate parameter for PAH)", phiPAH_default)
            xc = ask_user("XC (Conversion factor of H into C)", xc_default)
            return zeta, geff, phiPAH, xc
        end
        @warn "tasks.ne_dm_em.prompt_constants=true but stdin is not interactive; using config defaults."
    end

    return zeta_default, geff_default, phiPAH_default, xc_default
end

"""
    run_ne_dm_em_job(...)

Computes `ne` from density/temperature with Wolfire approximation, then
integrates along LOS to produce DM and EM maps.
"""
function run_ne_dm_em_job(cfg)::Dict{String,Any}
    enabled = DepolLib.task_enabled(cfg, ["tasks", "ne_dm_em", "enabled"]; default=true)
    if !enabled
        return DepolLib.skipped_job_result("ne_dm_em", "disabled by tasks.ne_dm_em.enabled")
    end

    sim_root = DepolLib.resolve_simulations_root(cfg)
    sim = string(DepolLib.cfg_require(cfg, ["simulation", "name"]))
    los = DepolLib.require_los(string(DepolLib.cfg_require(cfg, ["simulation", "los"])))

    zeta, Geff, phiPAH, XC = resolve_wolfire_constants(cfg)
    lbox_pc = Float64(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "lbox_pc"];
                              default=DepolLib.cfg_get(cfg, ["tasks", "reversals_map", "lbox_pc"]; default=50.0)))
    pc_to_cm = Float64(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "pc_to_cm"]; default=3.085677581491367e18))
    overwrite = Bool(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "overwrite"]; default=true))
    save_plot = Bool(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "save_plot"]; default=true))
    title_size = Int(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "title_size"]; default=44))
    label_size = Int(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "label_size"]; default=36))
    tick_size = Int(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "tick_size"]; default=30))
    tick_step_pc = Float64(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "tick_step_pc"]; default=10.0))
    cbar_label_size = Int(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "cbar_label_size"]; default=44))
    cbar_tick_size = Int(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "cbar_tick_size"]; default=36))
    cbar_width = Int(DepolLib.cfg_get(cfg, ["tasks", "ne_dm_em", "cbar_width"]; default=60))

    dens_file = DepolLib.simulation_field_path(sim_root, sim, "density")
    temp_file = DepolLib.simulation_field_path(sim_root, sim, "temperature")
    DepolLib.require_existing_files([dens_file, temp_file]; context="ne_dm_em_job")

    n = DepolLib.read_fits_f32(dens_file)
    T = DepolLib.read_fits_f32(temp_file)
    DepolLib.require_ndims(n, 3, "density")
    DepolLib.require_ndims(T, 3, "temperature")
    DepolLib.require_same_size([n, T], ["density", "temperature"])

    ne = Float32.(DepolLib.Wolfire_ne(zeta, Geff, phiPAH, XC, T, n))
    dm, em, dl_pc = DepolLib.dm_em_maps(ne, los; lbox_pc=lbox_pc)
    dm = Float32.(dm)
    em = Float32.(em)
    ne_col = Float32.(dm .* pc_to_cm)  # N_e = ∫ n_e dl in cm^-2

    syn_dir = DepolLib.synchrotron_dir(sim_root, sim, los)
    mkpath(syn_dir)

    ne_out = DepolLib.synchrotron_path(sim_root, sim, los, "ne.fits")
    dm_out = DepolLib.synchrotron_path(sim_root, sim, los, "DM.fits")
    em_out = DepolLib.synchrotron_path(sim_root, sim, los, "EM.fits")
    ne_col_out = DepolLib.synchrotron_path(sim_root, sim, los, "NeColumn.fits")

    DepolLib.write_FITS(ne_out, ne; overwrite=overwrite)
    DepolLib.write_FITS(dm_out, dm; overwrite=overwrite)
    DepolLib.write_FITS(em_out, em; overwrite=overwrite)
    DepolLib.write_FITS(ne_col_out, ne_col; overwrite=overwrite)

    dm_em_plot_out = nothing
    dm_plot_out = nothing
    em_plot_out = nothing
    ne_col_plot_out = nothing
    if save_plot
        dm_em_plot_out = DepolLib.standard_output_path(cfg, "ne_dm_em", "dm_em_maps", "pdf"; simu=sim, los=los)
        dm_plot_out = DepolLib.standard_output_path(cfg, "ne_dm_em", "dm_heatmap", "pdf"; simu=sim, los=los)
        em_plot_out = DepolLib.standard_output_path(cfg, "ne_dm_em", "em_heatmap", "pdf"; simu=sim, los=los)
        ne_col_plot_out = DepolLib.standard_output_path(cfg, "ne_dm_em", "ne_column_density_heatmap", "pdf"; simu=sim, los=los)

        _plot_dm_em_maps(dm, em, dm_em_plot_out;
            los=los, lbox_pc=lbox_pc,
            title_size=title_size, label_size=label_size, tick_size=tick_size,
            cbar_label_size=cbar_label_size, cbar_tick_size=cbar_tick_size,
            cbar_width=cbar_width, tick_step_pc=tick_step_pc)
        _plot_measure_heatmap(dm, dm_plot_out;
            los=los, lbox_pc=lbox_pc,
            title=L"\mathrm{Dispersion\ Measure}\ (DM)",
            cbar_label=L"DM\ [\mathrm{pc}\ \mathrm{cm}^{-3}]",
            colormap=:viridis,
            title_size=title_size, label_size=label_size, tick_size=tick_size,
            cbar_label_size=cbar_label_size, cbar_tick_size=cbar_tick_size,
            cbar_width=cbar_width, tick_step_pc=tick_step_pc)
        _plot_measure_heatmap(em, em_plot_out;
            los=los, lbox_pc=lbox_pc,
            title=L"\mathrm{Emission\ Measure}\ (EM)",
            cbar_label=L"EM\ [\mathrm{pc}\ \mathrm{cm}^{-6}]",
            colormap=:magma,
            title_size=title_size, label_size=label_size, tick_size=tick_size,
            cbar_label_size=cbar_label_size, cbar_tick_size=cbar_tick_size,
            cbar_width=cbar_width, tick_step_pc=tick_step_pc)
        _plot_measure_heatmap(ne_col, ne_col_plot_out;
            los=los, lbox_pc=lbox_pc,
            title=L"\mathrm{Electron\ Column\ Density}\ (N_e)",
            cbar_label=L"N_e\ [\mathrm{cm}^{-2}]",
            colormap=:plasma,
            title_size=title_size, label_size=label_size, tick_size=tick_size,
            cbar_label_size=cbar_label_size, cbar_tick_size=cbar_tick_size,
            cbar_width=cbar_width, tick_step_pc=tick_step_pc)
    end

    return Dict(
        "task" => "ne_dm_em",
        "simulation" => sim,
        "los" => los,
        "constants" => Dict(
            "zeta" => zeta,
            "Geff" => Geff,
            "phiPAH" => phiPAH,
            "XC" => XC,
        ),
        "dl_pc" => dl_pc,
        "outputs" => Dict(
            "ne_fits" => ne_out,
            "dm_fits" => dm_out,
            "em_fits" => em_out,
            "ne_column_fits" => ne_col_out,
            "dm_em_plot" => dm_em_plot_out,
            "dm_heatmap_plot" => dm_plot_out,
            "em_heatmap_plot" => em_plot_out,
            "ne_column_heatmap_plot" => ne_col_plot_out,
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    DepolLib.run_job_entrypoint("ne_dm_em", run_ne_dm_em_job)
end
