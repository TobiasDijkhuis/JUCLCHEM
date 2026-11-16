using PhysicalConstants: AtomicMassConstants, PlankConstant, k_B, ħ

"""
    get_tunneling_probability(E_reac, barrier_width, μ_ij)

``e^{-2a/\\hbar\\sqrt{2\\mu_{ij}E_{\\mathrm{reac}}}}``"""
function get_tunneling_probability(
    E_reac::AbstractFloat,
    barrier_width::AbstractFloat,
    μ_ij::AbstractFloat,
)::AbstractFloat
    return exp(
        -2 * barrier_width * 1.0e-10 / ħ * √(2 * μ_ij * AtomicMassConstants * E_reac * k_B),
    )
end

"""
    get_arrhenius_probability(E_reac, temperature)

``e^{-E/T}``"""
function get_arrhenius_probability(
    E_reac::AbstractFloat,
    temperature::AbstractFloat,
)::AbstractFloat
    return exp(-E_reac / temperature)
end

"""
    get_HH_prefactor(E_bind, mass)

``\\sqrt{2*SURFACE_SITE_DENSITY*E_{\\mathrm{bind}}/(π^2m)}``
"""
function get_HH_prefactor(E_bind::AbstractFloat, mass::AbstractFloat)::AbstractFloat
    return √(2 * SURFACE_SITE_DENSITY * E_bind * k_B / (mass*AtomicMassConstants*π^2))
end

"""
    get_thermal_desorption_rate_constant(E_bind, ν_des, T_dust)

The thermal desorption rate constant,

``\\nu_{\\mathrm{des}}e^{-E_{\\mathrm{bind}}/T_{\\mathrm{dust}}}``
"""
function get_thermal_desorption_rate_constant(
    E_bind::AbstractFloat,
    ν_des::AbstractFloat,
    T_dust::AbstractFloat,
)::AbstractFloat
    return ν_des * get_arrhenius_probability(E_bind, T_dust)
end

"""
    get_thermal_velocity(mass, T)


``v_{\\mathrm{th}}=\\sqrt{\\frac{8k_\\mathrm{B}T}{πm}}``
"""
function get_thermal_velocity(mass::AbstractFloat, T::AbstractFloat)::AbstractFloat
    thermal_vel = √(8*k_B*T/(π*AtomicMassConstants*mass))
    return thermal_vel
end

"""
    get_freezeout_rate_constant(α, β, mass, T_gas, freeze_factor=1)

``\\alpha(1+\\beta*1.671e-4/(r_{\\mathrm{grain}}T_{\\mathrm{gas}})*freeze_factor*v_{\\mathrm{th}}*GRAIN_CROSSSECTION_PER_H``
"""
function get_freezeout_rate_constant(
    α::AbstractFloat,
    β::AbstractFloat,
    mass::AbstractFloat,
    T_gas::AbstractFloat,
    freeze_factor::AbstractFloat = 1.0,
)::AbstractFloat
    if freeze_factor == 0
        return 0
    end

    thermal_vel = get_thermal_velocity(mass, T)
    rate_constant =
        α*(1+β*16.71e-4/(GRAIN_RADIUS*T_gas))*freeze_factor*thermal_vel*GRAIN_CROSSSECTION_PER_H
    return rate_constant
end

"""
    get_sticking_coefficient(S_0, T_critical, T_gas)

Sticking coefficient for freeze out taken from Chaabouni et al. 2012 A&A 538 Equation 1

``S_T = S_0(1+\\beta T_{\\mathrm{gas}}/T_{\\mathrm{crit}})/((1+T_{\\mathrm{gas}}/T_{\\mathrm{crit}})^\\beta)``
"""
function get_sticking_coefficient(
    S_0::AbstractFloat,
    T_critical::AbstractFloat,
    T_gas::AbstractFloat,
)::AbstractFloat
    β=2.5
    temperature_ratio = T_gas/T_critical
    sticking_coefficient = S_0*(1+β*temperature_ratio)/((1+temperature_ratio)^β)
    return sticking_coefficient
end

"""
    get_number_of_monolayers(abundance, num_sites_per_grain, gas_dust_density_ratio)

Calculated by ``n*density/(N_DUST_PARTICLE_DENSITY*NUM_LAYERS_PER_GRAIN)``
"""
function get_number_of_monolayers(
    abundance::AbstractFloat,
    num_sites_per_grain::AbstractFloat,
    gas_dust_density_ratio::AbstractFloat,
)::AbstractFloat
    return abundance*gas_dust_density_ratio / num_sites_per_grain
end

"""
    get_total_elemental_abuncance(abundances, elements_per_species)

"""
function get_total_elemental_abuncance(
    abundances::AbstractVector{AbstractFloat},
    elements_per_species::AbstractMatrix{Int},
)::AbstractArray{AbstractFloat}
    return elements_per_species * abundances
end

include("constants.jl")

"""
    construct_elements_per_species_matrix(constituent_dicts, elements)

"""
function construct_elements_per_species_matrix(
    constituent_dicts::AbstractVector{AbstractDict{AbstractString,Int}},
    elements::AbstractVector{AbstractString},
)::AbstractMatrix{Int}
    num_species = length(constituent_dicts)
    num_elements = length(elements)
    element_matrix = zeros(Int, (num_elements, num_species))
    for (species_idx, species_dict) in enumerate(constituent_dicts)
        for (element_idx, element) in enumerate(elements)
            element_matrix[element_idx, species_idx] = species_dict.get(element, 0)
        end
    end
    return element_matrix
end
