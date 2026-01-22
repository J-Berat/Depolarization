module DepolarizationConstants

using Printf

module Btransition
const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS = "y"
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"

const OUT_PLOTS = "/Users/jb270005/Desktop/Depolarization_canals/sightline_plots"
const SAVE_PLOT = true
const SAVE_EXT = "pdf"

const WRITE_LOG = true
const WRITE_CSV = true

const NPIX = 256
const LBOX_PC = 50.0
const B_SCALE = 1000.0
const RNG_SEED = nothing

const PhiArray = collect(-10.0:0.25:10.0)
const NPHI = length(PhiArray)

const LABSIZE = 25
const TICKSIZE = 25
const ARROW_FONTSIZE = 16

const SMOOTH_B = true
const SMOOTH_WIN = 5

const SIGN_EPS = 0.0
const DERIV_TOL = 0.05

const STOP_AT_NEXT_REVERSAL = true
const MAX_HALF_WIDTH_PIX = 18
const ISOLATE_FALLBACK = true

const B_LINEWIDTH = 3

const BOUNDS_STYLE = :dot
const BOUNDS_WIDTH = 1.2
const BOUNDS_ALPHA = 0.9

const DRAW_LABELS = true
const DRAW_ARROWS = true

const NE_COLOR = (:black, 0.85)
const NE_LW = 2.5
const NE_STYLE = :dashdot

const SCAT_ALPHA = 0.8
const SCAT_MARKER = :x
const SCAT_MSIZE = 7

const TICKALIGN_IN = 1.0
const TICKALIGN_OUT = 0.0

const MAP_MARK_COLOR = (:white, 0.6)
const MAP_MARK_STYLE = :dash
const MAP_MARK_LW = 2.5

const B0_COLOR = :blue
const B0_STYLE = :dash
const B0_LW = 2.0

const FDF_COLOR = :black
end

module ReversalsMap
const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS = "y"
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"

const NPIX = 256
const LBOX_PC = 50.0
const B_SCALE = 1000.0

const SMOOTH_B = true
const SMOOTH_WIN = 5
const SIGN_EPS = 0.0

const TARGET_NREV = 24

const OUT_DIR = "/Users/jb270005/Desktop/Depolarization_canals/sightline_plots"
const SAVE_PLOT = true
const SAVE_EXT = "pdf"

const LABSIZE = 25
const TICKSIZE = 25
end

module PmaxDecilesMeanSpectra
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"
const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS = "y"

const LBOX_PC = 50.0
const SAVE_FIG = true
const OUT_DIR = "/Users/jb270005/Desktop/Depolarization_canals"
const OUT_EXT = "pdf"

const USE_COLRANGE = false
const VMIN = 0.0
const VMAX = 12.0
end

module CutTransition
const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS = "y"
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"

const TRANSITIONS_CSV =
"/Users/jb270005/Desktop/Depolarization_canals/sightline_plots/BLOS_decile1_LOSy_transitions.csv"

const OUT_ROOT =
"/Users/jb270005/Desktop/Depolarization_canals/d1cf05bx10rms18000nograv1024"

const NPIX = 256
const LBOX_PC = 50.0

const MERGE_GAP_PIX = 1

const OVERWRITE = true

const VARS = ["Bx", "By", "Bz", "density", "Temperature", "Vx", "Vy", "Vz"]

const VAR_FILEMAP = Dict(
    "Bx" => "Bx.fits",
    "By" => "By.fits",
    "Bz" => "Bz.fits",
    "density" => "density.fits",
    "Temperature" => "Temperature.fits",
    "Vx" => "Vx.fits",
    "Vy" => "Vy.fits",
    "Vz" => "Vz.fits",
)
end

module MakeChunk
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"
const OUT_ROOT = "/Users/jb270005/Desktop/Depolarization_canals"

const CHUNK_SIZE = 32

const SIMU_LIST = [
    "d1cf00bx10rms18000nograv1024",
    "d1cf02bx10rms09000nograv1024",
    "d1cf02bx10rms18000nograv1024",
    "d1cf02bx10rms36000nograv1024",
    "d1cf02bx10rms72000nograv1024",
    "d1cf05bx10rms09000nograv1024",
    "d1cf05bx10rms18000nograv1024",
    "d1cf05bx10rms36000nograv1024",
    "d1cf05bx10rms72000nograv1024",
    "d1cf08bx10rms09000nograv1024",
    "d1cf08bx10rms18000nograv1024",
    "d1cf08bx10rms36000nograv1024",
    "d1cf08bx10rms72000nograv1024",
    "d1cf10bx10rms18000nograv1024",
]

const FIELDS = [
    "Bx", "By", "Bz",
    "Vx", "Vy", "Vz",
    "density", "temperature",
]
end

module PolarizationDegree
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES"
const SIMU_NAME = "d1cf05bx10rms18000nograv1024"
const LOS = "y"

const NPIX = 256
const LBOX_PC = 50.0

const CMAP = :viridis
const VMIN = 0.0
const VMAX = 1.0

const SYNCHRO_DIR = joinpath(
    SIMU_ROOT, SIMU_NAME, LOS, "Synchrotron", "WithFaraday"
)

const I_FILE = joinpath(SYNCHRO_DIR, "I.fits")
const Q_FILE = joinpath(SYNCHRO_DIR, "Q.fits")
const U_FILE = joinpath(SYNCHRO_DIR, "U.fits")

const OUT_DIR = joinpath(SYNCHRO_DIR, "Plots")
const OUT_PDF = joinpath(OUT_DIR, "p_degree_$(LOS).pdf")
end

module BVersusN
const DATA_DIR = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024"

const FILE_BX = "Bx.fits"
const FILE_BY = "By.fits"
const FILE_BZ = "Bz.fits"
const FILE_N = "density.fits"

const NBINS = 30
const MIN_PER_BIN = 50
const BINNING_IN_LOG = true
const PLOT_LOGLOG = true

const SAVE_FIG = true
const OUT_FIG = "B_vs_n_line.pdf"
end

module GradP
const LOS = "y"

const Q_path = LOS == "x" ?
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/x/Synchrotron/WithFaraday/Qnu.fits" :
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Qnu.fits"

const U_path = LOS == "x" ?
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/x/Synchrotron/WithFaraday/Unu.fits" :
    "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024/y/Synchrotron/WithFaraday/Unu.fits"

const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024"
const DENSITY_PATH = joinpath(SIMU_ROOT, "density.fits")

const Lbox_pc = 50.0
const NPIX = 256
const dx_pc = Lbox_pc / NPIX

const PC_TO_CM = 3.085677581e18
const dz_cm = dx_pc * PC_TO_CM

const P_eps = 1e-12

const ν_min_Hz = 120e6
const ν_max_Hz = 167e6
const dν_Hz = 98e3
const c_ms = 299_792_458.0

const OUTDIR = dirname(Q_path)
const SAVE_CUBES = true
const GRADP_OUT = joinpath(OUTDIR, "gradP_cube.fits")
const GRADRM_OUT = joinpath(OUTDIR, "gradRM_cube.fits")
const NH_OUT = joinpath(OUTDIR, "NH_map.fits")
const GRADNH_OUT = joinpath(OUTDIR, "grad_NH_map.fits")

const DO_PLOT_FOURPANELS = true
const k_plot = 10
const SAVE_PLOT = true
const PLOT_EXT = "pdf"
const PLOTPATH = joinpath(OUTDIR, @sprintf("gradP_gradRM_NH_gradNH_LOS_%s_k%04d.%s", LOS, k_plot, PLOT_EXT))

const GRADP_VMIN = 0.0
const GRADP_VMAX = 40.0
const GRADRM_VMIN = 0.0
const GRADRM_VMAX = 1.0

const NH_VMIN = 1e20
const NH_VMAX = 1e21
const GRADNH_VMIN = nothing
const GRADNH_VMAX = nothing

const LABELSIZE = 28
const TICKLABELSIZE = 22
const COLORBAR_LABELSIZE = 26
const COLORBAR_TICKLABELSIZE = 22
const NU_TITLESIZE = 26

const ticks_pc = (0:10:50, string.(0:10:50))

const DO_PDF_GRADP_NORM = true
const NBINS_PDF = 120
const XRANGE_PDF = (-6.0, 6.0)
const PDF_YLOG = false
const PDF_LABELSIZE = 28
const PDF_TICKLABELSIZE = 22
const PDF_LINEWIDTH = 4
const PDFPATH = joinpath(OUTDIR, @sprintf("PDF_gradPnorm_LOS_%s_k%04d.%s", LOS, k_plot, PLOT_EXT))
end

module HistogramDeciles
const SIMU_ROOT = "/Users/jb270005/Desktop/simu_RAMSES/d1cf05bx10rms18000nograv1024"

const LOS = "y"

const LBOX_PC = 50.0
const NPIX = 256
const PC_PER_PIX = LBOX_PC / NPIX

const VMIN = 0.0
const VMAX = 12.0
const COLOR_RANGE = (VMIN, VMAX)

const LABEL_SIZE = 26
const TICKLABEL_SIZE = 22
const TITLE_SIZE = 26
const LEGEND_SIZE = 22
const COLORBAR_SIZE = 24

const SAVE_PDF = true
const OUT_PDF = joinpath(SIMU_ROOT, "Pmax_deciles_histograms_LOS_$LOS.pdf")
end

module PmaxCut
const USE_FULLCUBE_COLORRANGE = false

const MAP_MARK_COLOR = (:white, 0.6)
const MAP_MARK_STYLE = :dash
const MAP_MARK_LW = 2.5

const i1 = 159
const j1 = 206

const LBOX_PC = 50.0

const TITLE_SIZE = 30
const XLABEL_SIZE = 34
const YLABEL_SIZE = 34
const XTICKLABEL_SIZE = 24
const YTICKLABEL_SIZE = 24
const CB_LABEL_SIZE = 30
const CB_TICKLABEL_SIZE = 24
const CB_WIDTH = 20
end

end
