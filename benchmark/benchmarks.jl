using BenchmarkTools
using JUCLCHEM:
    SpeciesInfo,
    ReactionInfo,
    REACTION_TYPE_LOOKUP,
    create_reaction_system,
    create_odeproblem,
    simulate,
    create_initial_conditions,
    read_species_file,
    read_reaction_file,
    PhysicsState
using Catalyst
using OrdinaryDiffEq: Rodas5P
using OrdinaryDiffEqBDF: QNDF, FBDF
using Sundials: CVODE_BDF
using ODEInterfaceDiffEq: radau

const SUITE = BenchmarkGroup()

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
    p4,
)::ReactionInfo
    type = REACTION_TYPE_LOOKUP[reaction_type]
    return ReactionInfo((
        r1 = r1,
        r2 = r2,
        type = type,
        p1 = p1,
        p2 = p2,
        p3 = p3,
        p4 = p4,
        α = 1.0,
        β = 1.0,
        γ = 1.0,
        T_min = 0.0,
        T_max = 10000.0,
    ))
end
species_infos = [make_dummy_species(name) for name in ["C", "O", "CO"]]
reaction_infos = [
    make_dummy_reaction(r1, r2, type, p1, p2, p3, p4)
    for (r1, r2, type, p1, p2, p3, p4) in [
        ("C", "O", "TWOBODY", "CO", missing, missing, missing),
    ]
]

SUITE["create_reaction_system"] = @benchmarkable create_reaction_system(
    $species_infos,
    $reaction_infos,
)

let
    network_path = joinpath((@__DIR__), "..", "test", "network")
    reaction_infos = read_reaction_file(joinpath(network_path, "reactions.csv"))
    species_infos = read_species_file(joinpath(network_path, "species.csv"))
    system = create_reaction_system(species_infos, reaction_infos)

    u0 = create_initial_conditions(Catalyst.species(system))
    physics_state = PhysicsState(
        T_gas = 10,
        T_dust = 10,
        A_V = 1,
        n_H = 1.0e4,
        radfield = 1,
        ζ = 1,
    )

    odeproblem = create_odeproblem(
        system,
        species_infos,
        reaction_infos,
        physics_state,
        u0 = u0,
    )
    for solver in [Rodas5P(), QNDF(), FBDF(), CVODE_BDF(), radau()]
        name = string(nameof(typeof(solver)))
        SUITE["model"][name] = @benchmarkable simulate($odeproblem, solver = $solver)
    end
end
