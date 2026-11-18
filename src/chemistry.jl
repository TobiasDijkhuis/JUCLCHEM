"""
    get_freezeout_rate_constant(α, β, mass, T_gas, freeze_factor=1)

``\\alpha(1+\\beta*1.671e-4/(r_{\\mathrm{grain}}T_{\\mathrm{gas}})*freeze_factor*v_{\\mathrm{th}}*GRAIN_CROSSSECTION_PER_H``
"""
function get_freezeout_rate_constant(
    α::Real,
    β::Real,
    mass::Real,
    T_gas::Real,
    freeze_factor::Real = 1.0,
)::Real
    thermal_vel = get_thermal_velocity(mass, T_gas)
    rate_constant = α *
        (1 + β * 16.71e-4 / (GRAIN_RADIUS * T_gas)) *
        freeze_factor *
        thermal_vel *
        GRAIN_CROSSSECTION_PER_H
    return rate_constant
end

"""
    get_sticking_coefficient(S_0, T_critical, T_gas)

Sticking coefficient for freeze out taken from Chaabouni et al. 2012 A&A 538 Equation 1

``S_T = S_0(1+\\beta T_{\\mathrm{gas}}/T_{\\mathrm{crit}})/((1+T_{\\mathrm{gas}}/T_{\\mathrm{crit}})^\\beta)``
"""
function get_sticking_coefficient(S_0::Real, T_critical::Real, T_gas::Real)::Real
    β = 2.5
    temperature_ratio = T_gas / T_critical
    sticking_coefficient = S_0 * (1 + β * temperature_ratio) / ((1 + temperature_ratio)^β)
    return sticking_coefficient
end

"""
    get_number_of_monolayers(abundance)

Calculated by ``N/(N_{\\mathrm{sites-per-grain}}(N_{\\mathrm{dust}}/N_{\\mathrm{gas}})``
"""
function get_number_of_monolayers(abundance::Real)::Real
    return abundance / (NUM_SITES_PER_GRAIN * GAS_DUST_NUMBER_DENSITY_RATIO)
end

"""
    get_total_elemental_abuncance(abundances, elements_per_species)

"""
function get_total_elemental_abuncance(
    abundances::AbstractVector{Real},
    elements_per_species::AbstractMatrix{Int},
)::AbstractArray{Real}
    return elements_per_species * abundances
end

"""
    construct_elements_per_species_matrix(constituent_dicts, elements)

"""
function construct_elements_per_species_matrix(
    constituent_dicts::AbstractVector{AbstractDict{AbstractString, Int}},
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
