include(joinpath(@__DIR__, "lib", "DepolLib.jl"))
using .DepolLib

# Legacy compatibility shim.
const DEFAULT_CONFIG_PATH = DepolLib.DEFAULT_CONFIG_PATH

const load_runtime_config = DepolLib.load_runtime_config
const build_cfg = DepolLib.build_cfg
const cfg_get = DepolLib.cfg_get
const cfg_require = DepolLib.cfg_require

const require_los = DepolLib.require_los
const los_axis = DepolLib.los_axis
const axis_pc = DepolLib.axis_pc
const ticks_pc = DepolLib.ticks_pc

const read_intervals_from_csv = DepolLib.read_intervals_from_csv
const merge_intervals = DepolLib.merge_intervals

const normalize_task_name = DepolLib.normalize_task_name
const standard_output_dir = DepolLib.standard_output_dir
const standard_output_path = DepolLib.standard_output_path

const finite_values = DepolLib.finite_values
const finite_minmax = DepolLib.finite_minmax
const resolve_simulations_root = DepolLib.resolve_simulations_root
const require_existing_files = DepolLib.require_existing_files
const read_fits_f32 = DepolLib.read_fits_f32
const require_ndims = DepolLib.require_ndims
const require_same_size = DepolLib.require_same_size

const read_FITS = DepolLib.read_FITS
const write_FITS = DepolLib.write_FITS
const los_config = DepolLib.los_config
const smooth_moving_average = DepolLib.smooth_moving_average
const sign_eps = DepolLib.sign_eps
const reversal_indices = DepolLib.reversal_indices
