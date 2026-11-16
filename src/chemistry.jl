using PhysicalConstants: AtomicMassConstants, PlankConstant, k_B, ħ

"""
    get_tunneling_probability(E_reac, barrier_width, ``\\mu``_ij)

``e^{-2a/ħ√{2μ_{ij}E_{reac}}}``"""
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
