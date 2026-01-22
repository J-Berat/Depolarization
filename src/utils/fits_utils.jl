using FITSIO

read_fits_array(path::AbstractString) = FITS(path, "r") do f
    read(f[1])
end

function write_fits_array(path::AbstractString, data; overwrite::Bool=true)
    if overwrite && isfile(path)
        rm(path; force=true)
    end
    FITS(path, "w") do f
        write(f, data)
    end
end
