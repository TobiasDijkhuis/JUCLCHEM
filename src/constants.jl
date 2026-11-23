using PhysicalConstants.CODATA2022: AtomicMassConstant, k_B, ħ
using Unitful

const AMU_SI = ustrip(u"kg", AtomicMassConstant)
const K_B_SI = ustrip(u"J/K", k_B)
const ħ_SI = ustrip(u"J*s/rad", ħ)

#Grain area per h nuclei, values taken from Cazaux & Tielens 2004 via UCL-PDR to match H2 formation rate (in m^2)
#
# Maybe replace by actual calculation of crossection from grain radius?
const GRAIN_CROSSSECTION_PER_H = 0.5 * (7.908e-26 + 8.473e-26)
const GRAIN_SURFACEAREA_PER_H = 4.0 * GRAIN_CROSSSECTION_PER_H

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
const GAS_DUST_NUMBER_DENSITY_RATIO = AMU_SI / (GRAIN_MASS * GAS_DUST_MASS_DENSITY_RATIO)

"""Number of sites on grain surface (m``^{-1}``)"""
const GRAIN_SURFACE_SITE_DENSITY = 1.5e19

"""Number of sites on a single grain."""
const NUM_SITES_PER_GRAIN = GRAIN_SURFACE_SITE_DENSITY * GRAIN_AREA

"""Conversion between Habing and Draine"""
const HABING_PER_DRAINE = 1 / 1.7

"""Minimum abundance."""
const MIN_ABUNDANCE = 1.0e-30
