using PhysicalConstants.CODATA2022: AtomicMassConstant

"""Number of seconds in one year."""
const SECONDS_PER_YEAR = 60 * 60 * 24 * 365

"""Elements in the network"""
const ELEMENTS::Vector{String} = ["H", "D", "He", "C", "O"]

"""Default Cosmic Ray Ionization Rate (s``^{-1}``)"""
const ζ_0 = 1.3e-17

"""Grain radius in m"""
const GRAIN_RADIUS = 1.0e-7

"""Grain density in kg m``^{-3}``"""
const GRAIN_DENSITY = 3.0e3

"""Area of a single grain in m``^2``."""
const GRAIN_AREA = 4 * π * GRAIN_RADIUS^2

"""Volume of a single grain in kg."""
const GRAIN_VOLUME = 4 / 3 * π * GRAIN_RADIUS^3

"""Mass of a single grain in kg."""
const GRAIN_MASS = GRAIN_VOLUME * GRAIN_DENSITY

"""``m_{\\mathrm{gas}}/m_{\\mathrm{dust}}``"""
const GAS_DUST_MASS_DENSITY_RATIO = 100.0

"""``N_{\\mathrm{gas}}/N_{\\mathrm{dust}}``

Note that this is defined as the reciprocal as it is in UCLCHEM."""
const GAS_DUST_NUMBER_DENSITY_RATIO = AtomicMassConstant /
    (GRAIN_MASS * GAS_DUST_MASS_DENSITY_RATIO)

"""Number of sites on grain surface (m``^{-1}``)"""
const GRAIN_SURFACE_SITE_DENSITY = 1.5e19

"""Number of sites on a single grain."""
const NUM_SITES_PER_GRAIN = GRAIN_SURFACE_SITE_DENSITY * GRAIN_AREA

"""Conversion between Habing and Draine"""
const HABING_PER_DRAINE = 1 / 1.7
