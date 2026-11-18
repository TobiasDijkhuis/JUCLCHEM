using CSV
using Catalyst
using Symbolics

function collapse(names::AbstractVector{<:AbstractString}, sp::AbstractDict)
    isempty(names) && return nothing, nothing
    u = unique(names)
    return [sp[n] for n in u], [count(==(n), names) for n in u]
end

function get_symbolic_rate(reaction_type::ReactionType, idx::Int, P)::Union{
    Symbolics.Num,
    Nothing,
}
    if reaction_type == TWOBODY
        rate = get_twobody_rate_constant(P.α[idx], P.β[idx], P.γ[idx], P.T_gas)
    elseif reaction_type == IONOPOL1
        rate = get_ionopol1_rate_constant(P.α[idx], P.β[idx], P.γ[idx], P.T_gas)
    elseif reaction_type == IONOPOL2
        rate = get_ionopol2_rate_constant(P.α[idx], P.β[idx], P.γ[idx], P.T_gas)
    elseif reaction_type == PHOTON
        rate = get_photo_rate_constant(P.α[idx], P.β[idx], P.radfield, P.A_V)
    elseif reaction_type == CRPHOT
        rate = get_crp_photo_rate_constant(
            P.α[idx],
            P.β[idx],
            P.γ[idx],
            P.ω_dust,
            P.ζ,
            P.T_gas,
        )
    else
        return nothing
        # throw("Unimplemented reaction type $(info.type)")
    end
    return rate
end

function create_catalyst_reaction(info::ReactionInfo, idx::Int, sp, P)::Union{
    Reaction,
    Nothing,
}
    rate = get_symbolic_rate(info.type, idx, P)
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
    @debug "Initialing reaction system with $(length(reactions_list)) reactions and $(length(species_list)) species."

    t = default_t()
    @parameters begin
        α[1:length(reactions_list)] = [r.α for r in reactions_list]
        β[1:length(reactions_list)] = [r.β for r in reactions_list]
        γ[1:length(reactions_list)] = [r.γ for r in reactions_list]
        A_V = 1.0
        T_gas = 10.0
        T_dust = 10.0
        n_H = 1.0e3
        radfield = 1.0
        ζ = 1.0
        ω_dust = DEFAULT_ω_DUST
    end
    P = (; α, β, γ, A_V, T_gas, T_dust, n_H, radfield, ζ, ω_dust)

    species_symbols = [only(@species $(Symbol(s.name))(t)) for s in species_list]
    @debug "Initialized $(length(species_symbols)) symbolic species"
    sp = Dict(s.name => x for (s, x) in zip(species_list, species_symbols))

    rxs = filter(
        !isnothing,
        map(
            (r, idx) -> create_catalyst_reaction(r, idx, sp, P),
            reactions_list,
            eachindex(reactions_list),
        ),
    )
    @debug "Initialized $(length(rxs)) symbolic reactions."

    @named network = ReactionSystem(
        rxs,
        t,
        species_symbols,
        [A_V, T_gas, T_dust, n_H, radfield, ζ, α, β, γ, ω_dust],
    )

    @debug "Initialized reaction system has $(numspecies(network)) species and $(numreactions(network)) reactions"
    return complete(network)
end
