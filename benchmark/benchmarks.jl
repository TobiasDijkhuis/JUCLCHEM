using BenchmarkTools
using JUCLCHEM: SpeciesInfo, ReactionInfo, REACTION_TYPE_LOOKUP, create_reaction_system

function make_dummy_species(name::String)::SpeciesInfo
    return SpeciesInfo((name = name, mass = 1.0, E_bind = 1.0, E_diff = 0.5, ΔH_F = 1.0))
end
function make_dummy_reaction(
    r1::String,
    r2::Union{String, Missing},
    reaction_type,
    p1,
    p2,
    p3,
)::ReactionInfo
    type = REACTION_TYPE_LOOKUP[reaction_type]
    return ReactionInfo((
        r1 = r1,
        r2 = r2,
        type = type,
        p1 = p1,
        p2 = p2,
        p3 = p3,
        α = 1.0,
        β = 1.0,
        γ = 1.0,
    ))
end
species_infos = [make_dummy_species(name) for name in ["C", "O", "CO"]]
reaction_infos = [
    make_dummy_reaction(r1, r2, type, p1, p2, p3)
    for (r1, r2, type, p1, p2, p3) in [("C", "O", "TWOBODY", "CO", missing, missing)]
]

b = @benchmarkable create_reaction_system(species_infos, reaction_infos)
tune!(b)
run(b)
