"""
    ReadSimulation(simu::String, LOS::String, conversionn::Number, conversionT::Number, conversionV::Number, conversionB::Number)
        -> Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, Union{AbstractArray, Nothing}, Union{AbstractArray, Nothing}}

Reads and processes simulation data from FITS files for a specified line of sight (LOS).

# Arguments
- `simu::String`: The directory containing the simulation FITS files.
- `LOS::String`: The line of sight direction, either "x", "y", or "z".
- `conversionn::Number`: Conversion factor for density.
- `conversionT::Number`: Conversion factor for temperature.
- `conversionV::Number`: Conversion factor for velocity.
- `conversionB::Number`: Conversion factor for magnetic field.

# Returns
- `Tuple`: A tuple containing the following elements:
  - `B1::AbstractArray`: Magnetic field component perpendicular to the LOS.
  - `B2::AbstractArray`: Another magnetic field component perpendicular to the LOS.
  - `BLOS::AbstractArray`: Magnetic field component along the LOS.
  - `V1::AbstractArray`: Velocity component perpendicular to the LOS.
  - `V2::AbstractArray`: Another velocity component perpendicular to the LOS.
  - `VLOS::AbstractArray`: Velocity component along the LOS.
  - `T::AbstractArray`: Temperature array.
  - `n::AbstractArray`: Density array.
  - `nH2::Union{AbstractArray, Nothing}`: Molecular hydrogen density array (or `nothing` if the file is not present).
  - `nHp::Union{AbstractArray, Nothing}`: Ionized hydrogen density array (or `nothing` if the file is not present).
"""
function ReadSimulation(simu, LOS, conversionn, conversionT, conversionV, conversionB)
    paths = simulation_paths(simu)

    T = read_file(paths.temperature, conversionT)
    n = read_file(paths.density, conversionn)
    nH2 = read_optional_file(paths.densityH2, conversionn, LOS)
    nHp = read_optional_file(paths.densityHp, conversionn, LOS)

    B1, B2, BLOS = read_magnetic_components(paths, LOS, conversionB)
    V1, V2, VLOS = read_velocity_components(paths, LOS, conversionV)

    B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
    V1, V2, VLOS = permute_dims(V1, LOS), permute_dims(V2, LOS), permute_dims(VLOS, LOS)
    T, n = permute_dims(T, LOS), permute_dims(n, LOS)

    return (B1, B2, BLOS, V1, V2, VLOS, T, n, nH2, nHp)
end

function ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    paths = simulation_paths(simu)

    T = read_file(paths.temperature, conversionT)
    n = read_file(paths.density, conversionn)
    nH2 = read_optional_file(paths.densityH2, conversionn, LOS)
    nHp = read_optional_file(paths.densityHp, conversionn, LOS)

    B1, B2, BLOS = read_magnetic_components(paths, LOS, conversionB)

    B1, B2, BLOS = permute_dims(B1, LOS), permute_dims(B2, LOS), permute_dims(BLOS, LOS)
    T, n = permute_dims(T, LOS), permute_dims(n, LOS)

    return (B1, B2, BLOS, T, n, nH2, nHp)
end

function simulation_paths(simu)
    return (
        bx = "$simu/Bx.fits",
        by = "$simu/By.fits",
        bz = "$simu/Bz.fits",
        density = "$simu/density.fits",
        temperature = "$simu/temperature.fits",
        vx = "$simu/Vx.fits",
        vy = "$simu/Vy.fits",
        vz = "$simu/Vz.fits",
        densityH2 = "$simu/densityH2.fits",
        densityHp = "$simu/densityHp.fits",
    )
end

function read_magnetic_components(paths, LOS, conversionB)
    Bx = read_file(paths.bx, conversionB)
    By = read_file(paths.by, conversionB)
    Bz = read_file(paths.bz, conversionB)
    return select_components(LOS, Bx, By, Bz)
end

function read_velocity_components(paths, LOS, conversionV)
    Vx = read_file(paths.vx, conversionV)
    Vy = read_file(paths.vy, conversionV)
    Vz = read_file(paths.vz, conversionV)
    return select_components(LOS, Vx, Vy, Vz)
end

function select_components(LOS, comp_x, comp_y, comp_z)
    if LOS == "z"
        return (comp_x, comp_y, comp_z)
    elseif LOS == "y"
        return (comp_z, comp_x, comp_y)
    elseif LOS == "x"
        return (comp_y, comp_z, comp_x)
    else
        error("LOS must be x/y/z")
    end
end
