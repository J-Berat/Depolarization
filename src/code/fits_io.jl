using FITSIO

"""
    read_FITS(...)

    Reads the primary HDU from a FITS file.
"""
function read_FITS(path::AbstractString)
    FITS(path, "r") do f
        read(f[1])
    end
end

"""
    write_FITS(...)

    Writes array data to FITS.
"""
function write_FITS(path::AbstractString, data; overwrite::Bool=true)
    if overwrite && isfile(path)
        rm(path; force=true)
    end
    FITS(path, "w") do f
        write(f, data)
    end
end
