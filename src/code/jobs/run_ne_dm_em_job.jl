using CairoMakie
using LaTeXStrings

include(joinpath(@__DIR__, "..", "lib", "DepolLib.jl"))
using .DepolLib

"""
    _sky_plane_labels(...)

Returns axis labels for the sky plane associated with LOS.
"""
function _sky_plane_labels(los::AbstractString)
    los = require_los(los)
    if los == "x"
        return "y", "z"
    elseif los == "y"
        return "x", "z"
    else
        return "x", "y"
    end
end

"""
    _plot_dm_em_maps(...)

Writes a side-by-side DM/EM map figure.
"""
function _plot_dm_em_maps(dm::AbstractMatrix, em::AbstractMatrix, out_plot::AbstractString;
                          los::AbstractString, lbox_pc::Real=50.0)
    size(dm) == size(em) || error("DM/EM map size mismatch: dm=$(size(dm)) em=$(size(em))")
    xlab, ylab = _sky_plane_labels(los)

    x = axis_pc(size(dm, 2); lbox_pc=lbox_pc)
    y = axis_pc(size(dm, 1); lbox_pc=lbox_pc)

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(1900, 820))

        ax_dm = Axis(fig[1, 1],
            title=L"\mathrm{Dispersion\ Measure}\ (DM)",
            xlabel=string(xlab, " [pc]"),
            ylabel=string(ylab, " [pc]"),
            xgridvisible=false, ygridvisible=false)
        hm_dm = heatmap!(ax_dm, x, y, dm; colormap=:viridis)
        Colorbar(fig[1, 2], hm_dm; label=L"DM\ [\mathrm{pc}\ \mathrm{cm}^{-3}]")

        ax_em = Axis(fig[1, 3],
            title=L"\mathrm{Emission\ Measure}\ (EM)",
            xlabel=string(xlab, " [pc]"),
            ylabel=string(ylab, " [pc]"),
            xgridvisible=false, ygridvisible=false)
        hm_em = heatmap!(ax_em, x, y, em; colormap=:magma)
        Colorbar(fig[1, 4], hm_em; label=L"EM\ [\mathrm{pc}\ \mathrm{cm}^{-6}]")

        save(out_plot, fig)
    end

    return out_plot
end

"""
    run_ne_dm_em_job(...)

Computes `ne` from density/temperature with Wolfire approximation, then
integrates along LOS to produce DM and EM maps.
"""
function run_ne_dm_em_job(cfg)::Dict{String,Any}
    enabled = task_enabled(cfg, ["tasks", "ne_dm_em", "enabled"]; default=true)
    if !enabled
        return skipped_job_result("ne_dm_em", "disabled by tasks.ne_dm_em.enabled")
    end

    sim_root = resolve_simulations_root(cfg)
    sim = string(cfg_require(cfg, ["simulation", "name"]))
    los = require_los(string(cfg_require(cfg, ["simulation", "los"])))

    zeta = Float64(cfg_get(cfg, ["tasks", "ne_dm_em", "zeta"]; default=2.5e-16))
    Geff = Float64(cfg_get(cfg, ["tasks", "ne_dm_em", "Geff"]; default=1.0))
    phiPAH = Float64(cfg_get(cfg, ["tasks", "ne_dm_em", "phiPAH"]; default=0.5))
    XC = Float64(cfg_get(cfg, ["tasks", "ne_dm_em", "XC"]; default=1.4e-4))
    lbox_pc = Float64(cfg_get(cfg, ["tasks", "ne_dm_em", "lbox_pc"];
                              default=cfg_get(cfg, ["tasks", "reversals_map", "lbox_pc"]; default=50.0)))
    overwrite = Bool(cfg_get(cfg, ["tasks", "ne_dm_em", "overwrite"]; default=true))
    save_plot = Bool(cfg_get(cfg, ["tasks", "ne_dm_em", "save_plot"]; default=true))

    dens_file = simulation_field_path(sim_root, sim, "density")
    temp_file = simulation_field_path(sim_root, sim, "temperature")
    require_existing_files([dens_file, temp_file]; context="ne_dm_em_job")

    n = read_fits_f32(dens_file)
    T = read_fits_f32(temp_file)
    require_ndims(n, 3, "density")
    require_ndims(T, 3, "temperature")
    require_same_size([n, T], ["density", "temperature"])

    ne = Float32.(Wolfire_ne(zeta, Geff, phiPAH, XC, T, n))
    dm, em, dl_pc = dm_em_maps(ne, los; lbox_pc=lbox_pc)
    dm = Float32.(dm)
    em = Float32.(em)

    syn_dir = synchrotron_dir(sim_root, sim, los)
    mkpath(syn_dir)

    ne_out = synchrotron_path(sim_root, sim, los, "ne.fits")
    dm_out = synchrotron_path(sim_root, sim, los, "DM.fits")
    em_out = synchrotron_path(sim_root, sim, los, "EM.fits")

    write_FITS(ne_out, ne; overwrite=overwrite)
    write_FITS(dm_out, dm; overwrite=overwrite)
    write_FITS(em_out, em; overwrite=overwrite)

    plot_out = nothing
    if save_plot
        plot_out = standard_output_path(cfg, "ne_dm_em", "dm_em_maps", "pdf"; simu=sim, los=los)
        _plot_dm_em_maps(dm, em, plot_out; los=los, lbox_pc=lbox_pc)
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
            "dm_em_plot" => plot_out,
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_job_entrypoint("ne_dm_em", run_ne_dm_em_job)
end
