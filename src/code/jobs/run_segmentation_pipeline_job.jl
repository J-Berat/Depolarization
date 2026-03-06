using CairoMakie
using Printf

include(joinpath(@__DIR__, "..", "lib", "DepolLib.jl"))
using .DepolLib

"""
    split_along_dim(...)

    Splits data into regular chunks along one axis and writes outputs.
"""
function split_along_dim(data::AbstractArray, outdir::String; base::String, dim::Int, chunk_size::Int, task_name::String, sim::String, los::String)
    n = size(data, dim)
    n % chunk_size == 0 || error("Dimension $dim size=$n not divisible by chunk_size=$chunk_size")

    n_chunks = div(n, chunk_size)
    mkpath(outdir)

    for i in 1:n_chunks
        i1 = (i - 1) * chunk_size + 1
        i2 = i * chunk_size
        ranges = ntuple(d -> d == dim ? (i1:i2) : Colon(), ndims(data))
        slice = Array(@view data[ranges...])

        chunk_dir = joinpath(outdir, @sprintf("chunk%02d", i))
        mkpath(chunk_dir)
        outname = string(normalize_task_name(task_name), "_", sim, "_LOS", los, "_", base, ".fits")
        write_FITS(joinpath(chunk_dir, outname), slice)
    end
end

"""
    cut_cube(...)

    Slices a 3D cube between index bounds along a chosen axis.
"""
function cut_cube(A::Array{T,3}, ax::Int, kmin::Int, kmax::Int) where {T}
    if ax == 1
        return Array(@view A[kmin:kmax, :, :])
    elseif ax == 2
        return Array(@view A[:, kmin:kmax, :])
    else
        return Array(@view A[:, :, kmin:kmax])
    end
end

"""
    run_segmentation_pipeline_job(...)

    Runs the merged chunking + transition cuts + Pmax comparison pipeline.
"""
function run_segmentation_pipeline_job(cfg)::Dict{String,Any}
    enabled = task_enabled(cfg, ["tasks", "segmentation_pipeline_job", "enabled"]; default=true)
    if !enabled
        return skipped_job_result("segmentation_pipeline_job", "disabled by tasks.segmentation_pipeline_job.enabled")
    end

    sim_root = resolve_simulations_root(cfg)
    sim = string(cfg_require(cfg, ["simulation", "name"]))
    los = require_los(string(cfg_require(cfg, ["simulation", "los"])))

    npix = Int(cfg_get(cfg, ["tasks", "cut_transition", "npix"]; default=256))
    lbox_pc = Float64(cfg_get(cfg, ["tasks", "cut_transition", "lbox_pc"]; default=50.0))
    merge_gap = Int(cfg_get(cfg, ["tasks", "cut_transition", "merge_gap_pix"]; default=1))
    overwrite = Bool(cfg_get(cfg, ["tasks", "cut_transition", "overwrite"]; default=true))
    chunk_size = Int(cfg_get(cfg, ["tasks", "make_chunk", "chunk_size"]; default=32))

    transition_csv = joinpath(
        string(cfg_get(cfg, ["paths", "transitions_csv_root"]; default="./data/transitions")),
        string(cfg_get(cfg, ["tasks", "cut_transition", "transitions_csv"]; default="BLOS_decile1_LOSy_transitions.csv"))
    )

    fields = ["Bx", "By", "Bz", "Vx", "Vy", "Vz", "density", "temperature"]
    cut_vars = Dict(
        "Bx" => "Bx.fits",
        "By" => "By.fits",
        "Bz" => "Bz.fits",
        "density" => "density.fits",
        "temperature" => "temperature.fits",
        "Vx" => "Vx.fits",
        "Vy" => "Vy.fits",
        "Vz" => "Vz.fits",
    )

    required_files = [transition_csv]
    append!(required_files, [simulation_field_path(sim_root, sim, field) for field in fields])
    require_existing_files(required_files; context="segmentation_pipeline_job")

    # 1) Chunking
    chunk_out = standard_output_dir(cfg, "segmentation_pipeline_job_chunks"; simu=sim, los=los)
    size_tag = @sprintf("%dpix", chunk_size)
    out_x = joinpath(chunk_out, size_tag, "split_x")
    out_y = joinpath(chunk_out, size_tag, "split_y")
    out_z = joinpath(chunk_out, size_tag, "split_z")

    for field in fields
        infile = simulation_field_path(sim_root, sim, field)
        if !isfile(infile)
            continue
        end
        data = read_FITS(infile)
        require_ndims(data, 3, field)
        size(data) == (npix, npix, npix) || error("Unexpected cube size for $field: $(size(data))")

        split_along_dim(data, out_x; base=field, dim=1, chunk_size=chunk_size, task_name="segmentation_pipeline_job", sim=sim, los=los)
        split_along_dim(data, out_y; base=field, dim=2, chunk_size=chunk_size, task_name="segmentation_pipeline_job", sim=sim, los=los)
        split_along_dim(data, out_z; base=field, dim=3, chunk_size=chunk_size, task_name="segmentation_pipeline_job", sim=sim, los=los)
    end

    # 2) Cut from transition intervals
    intervals_raw = read_intervals_from_csv(transition_csv)
    intervals = merge_intervals(intervals_raw; gap_pix=merge_gap)
    ax = los_axis(los)

    cut_out = standard_output_dir(cfg, "segmentation_pipeline_job_cuts"; simu=sim, los=los)
    dist = collect(range(0.0, lbox_pc; length=npix))

    intervals_csv = standard_output_path(cfg, "segmentation_pipeline_job", "intervals", "csv"; simu=sim, los=los)
    open(intervals_csv, "w") do io
        println(io, "segment,kmin,kmax,smin_pc,smax_pc,width_pc,npix")
        for (seg, (kmin, kmax)) in enumerate(intervals)
            println(io, @sprintf("%d,%d,%d,%.6f,%.6f,%.6f,%d", seg, kmin, kmax, dist[kmin], dist[kmax], dist[kmax]-dist[kmin], kmax-kmin+1))
        end
    end

    for (seg, (kmin, kmax)) in enumerate(intervals)
        segdir = joinpath(cut_out, @sprintf("seg_%02d_k%04d_%04d", seg, kmin, kmax))
        mkpath(segdir)

        for (var, fname) in cut_vars
            inpath = simulation_field_path(sim_root, sim, first(splitext(fname)))
            if !isfile(inpath)
                continue
            end
            A = read_FITS(inpath)
            require_ndims(A, 3, var)
            Acut = cut_cube(A, ax, kmin, kmax)
            outname = string(normalize_task_name("segmentation_pipeline_job"), "_", sim, "_LOS", los, "_seg", lpad(seg, 2, '0'), "_", replace(var, " " => "_"), ".fits")
            outpath = joinpath(segdir, outname)
            write_FITS(outpath, Acut; overwrite=overwrite)
        end
    end

    # 3) Pmax comparison figure (first 3 cuts + full)
    pmax_maps = Array{Float32,2}[]
    titles = String[]
    for (seg, (kmin, kmax)) in enumerate(intervals[1:min(3, length(intervals))])
        segdir = joinpath(cut_out, @sprintf("seg_%02d_k%04d_%04d", seg, kmin, kmax))
        f = joinpath(segdir, string(normalize_task_name("segmentation_pipeline_job"), "_", sim, "_LOS", los, "_seg", lpad(seg, 2, '0'), "_Pmax.fits"))
        if isfile(f)
            push!(pmax_maps, read_fits_f32(f))
            push!(titles, @sprintf("segment %d", seg))
        end
    end

    full_pmax_file = withfaraday_path(sim_root, sim, los, "Pmax.fits")
    if isfile(full_pmax_file)
        push!(pmax_maps, read_fits_f32(full_pmax_file))
        push!(titles, "full")
    end

    fig_out = standard_output_path(cfg, "segmentation_pipeline_job", "pmax_panels", "pdf"; simu=sim, los=los)
    if !isempty(pmax_maps)
        with_theme(theme_latexfonts()) do
            n = length(pmax_maps)
            cols = min(n, 2)
            rows = cld(n, cols)
            fig = Figure(size=(700 * cols, 600 * rows))

            for i in 1:n
                r = cld(i, cols)
                c = ((i - 1) % cols) + 1
                axp = Axis(fig[r, c], title=titles[i], xlabel="Distance [pc]", ylabel="Distance [pc]", xgridvisible=false, ygridvisible=false)
                m = pmax_maps[i]
                x = axis_pc(size(m, 2); lbox_pc=lbox_pc)
                y = axis_pc(size(m, 1); lbox_pc=lbox_pc)
                hm = heatmap!(axp, x, y, m; colormap=:magma)
                Colorbar(fig[r, c + cols], hm)
            end
            save(fig_out, fig)
        end
    end

    return Dict(
        "task" => "segmentation_pipeline_job",
        "chunk_dir" => chunk_out,
        "cut_dir" => cut_out,
        "intervals_csv" => intervals_csv,
        "pmax_plot" => fig_out,
        "n_intervals" => length(intervals),
        "status" => "ok",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_job_entrypoint("segmentation_pipeline_job", run_segmentation_pipeline_job)
end
