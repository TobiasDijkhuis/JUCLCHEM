using CSV
using Catalyst
using Symbolics

"""
    get_symbolic_rate(reaction_info, reaction_idx, P, masses, species_indices)

Get the symbolic function to calculate the reaction rate of a reaction,
depending on its type.
"""
function get_symbolic_rate(
    reaction_info::ReactionInfo,
    reaction_idx::Int,
    P::NamedTuple,
    masses::Dict{String, Float64},
    species_indices::Dict{String, Int},
)::Union{Symbolics.Num, Nothing}
    reaction_type = reaction_info.type
    i1 = species_indices[reaction_info.r1]
    if !ismissing(reaction_info.r2)
        i2 = species_indices[reaction_info.r2]
    end
    if reaction_type == TWOBODY
        return get_twobody_rate_constant(
            P.α[reaction_idx],
            P.β[reaction_idx],
            P.γ[reaction_idx],
            P.T_gas,
        ) *
            P.n_H
    elseif reaction_type == IONOPOL1
        return get_ionopol1_rate_constant(
            P.α[reaction_idx],
            P.β[reaction_idx],
            P.γ[reaction_idx],
            P.T_gas,
        ) *
            P.n_H
    elseif reaction_type == IONOPOL2
        return get_ionopol2_rate_constant(
            P.α[reaction_idx],
            P.β[reaction_idx],
            P.γ[reaction_idx],
            P.T_gas,
        ) *
            P.n_H
    elseif reaction_type == PHOTON
        return get_photo_rate_constant(
            P.α[reaction_idx],
            P.γ[reaction_idx],
            P.radfield,
            P.A_V,
        )
    elseif reaction_type == CRPHOT
        return get_crp_photo_rate_constant(
            P.α[reaction_idx],
            P.β[reaction_idx],
            P.γ[reaction_idx],
            P.ω_dust,
            P.ζ,
            P.T_gas,
        )
    elseif reaction_type == CRP
        return get_crp_rate_constant(P.α[reaction_idx], P.ζ)
    elseif reaction_type == LH
        return P.α[reaction_idx] *
            get_langmuir_hinshelwood_rate_constant(
                P.γ[reaction_idx],
                P.T_dust,
                P.E_bind[i1],
                P.E_bind[i2],
                P.E_diff[i1],
                P.E_diff[i2],
                masses[reaction_info.r1],
                masses[reaction_info.r2],
            )
    elseif reaction_type == FREEZE
        return get_freezeout_rate_constant(
            P.α[reaction_idx],
            P.β[reaction_idx],
            masses[reaction_info.r1],
            P.T_gas,
            P.freeze_factor,
        ) *
            P.n_H^2
    elseif reaction_type == THERM
        ν_des = get_ν_HH(P.E_bind[i1], masses[reaction_info.r1])
        return get_thermal_desorption_rate_constant(
            P.α[reaction_idx],
            P.E_bind[i1],
            ν_des,
            P.T_dust,
        )
    elseif reaction_type == H2FORM
        return get_CT_H2_form_efficiency(P.T_gas, P.T_dust)
    else
        return nothing
    end
end

function create_catalyst_reaction(
    info,
    reaction_idx,
    sp,
    P,
    species_masses,
    species_indices,
)::Union{Reaction, Nothing}
    rate = get_symbolic_rate(info, reaction_idx, P, species_masses, species_indices)
    if isnothing(rate)
        return nothing
    end
    subs, sst = collapse(get_pure_reactants(info), sp)
    prds, pst = collapse(get_pure_products(info), sp)
    rxn = Reaction(rate, subs, prds, sst, pst)
    return rxn
end

DEFAULT_ω_DUST = 0.5

function create_reaction_system(
    species_list::AbstractVector{SpeciesInfo},
    reactions_list::AbstractVector{ReactionInfo},
)::ReactionSystem
    num_reactions = length(reactions_list)
    num_species = length(species_list)
    @debug "Initialing reaction system with $(num_reactions) reactions and $(num_species) species."

    t = default_t()
    ps = @parameters begin
        α[1:num_reactions] = [r.α for r in reactions_list]
        β[1:num_reactions] = [r.β for r in reactions_list]
        γ[1:num_reactions] = [r.γ for r in reactions_list]
        E_bind[1:num_species] = [s.E_bind for s in species_list]
        E_diff[1:num_species] = [s.E_diff for s in species_list]
        A_V
        T_gas
        T_dust
        n_H
        radfield
        ζ
        ω_dust
        freeze_factor
    end
    P = NamedTuple(Symbolics.getname(p) => p for p in ps)

    function is_positive_ion(name)::Bool
        return endswith(name, "+")
    end

    species_symbols = [only(@species $(Symbol(s.name))(t)) for s in species_list]
    species_indices = Dict(s.name => i for (i, s) in enumerate(species_list))
    masses = Dict(s.name => s.mass for s in species_list)
    masses["E-"] = 1.0 # Use mass 1 for E-, otherwise get mass 0 in freeze-out rate.

    @debug "Initialized $(length(species_symbols)) symbolic species"
    sp = Dict(s.name => x for (s, x) in zip(species_list, species_symbols))

    rxs = filter(
        !isnothing,
        map(
            (r, idx) -> create_catalyst_reaction(r, idx, sp, P, masses, species_indices),
            reactions_list,
            eachindex(reactions_list),
        ),
    )
    @debug "Initialized $(length(rxs)) symbolic reactions."

    @named network = ReactionSystem(rxs, t, species_symbols, collect(values(P)))

    @debug "Initialized reaction system has $(numspecies(network)) species, $(numreactions(network)) reactions and $(numparams(network)) parameters"
    return complete(network)
end
