using FITSIO

function read_FITS(path::AbstractString)
    FITS(path, "r") do f
        read(f[1])
    end
end

function write_FITS(path::AbstractString, data; overwrite::Bool=true)
    if overwrite && isfile(path)
        rm(path; force=true)
    end
    FITS(path, "w") do f
        write(f, data)
    end
end
