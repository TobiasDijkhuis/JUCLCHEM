using JUCLCHEM:
    read_species_file,
    read_reaction_file,
    create_reaction_system,
    create_odeproblem,
    simulate,
    create_initial_conditions,
    PhysicsState,
    SECONDS_PER_YEAR
using Plots
using Catalyst

test_network_dir = joinpath(@__DIR__(), "..", "test", "network")
species_path = joinpath(test_network_dir, "species.csv")
reactions_path = joinpath(test_network_dir, "reactions.csv")
species_list = read_species_file(species_path)
reactions_list = read_reaction_file(reactions_path)

system = create_reaction_system(species_list, reactions_list)

species_symbols = Catalyst.species(system)
u0 = create_initial_conditions(species_symbols)

ps = PhysicsState(
    T_dust = 10.0,
    T_gas = 10.0,
    A_V = 1.0,
    n_H = 1.0e6,
    radfield = 1.0,
    ζ = 1.0,
)
odeprob = create_odeproblem(system, species_list, reactions_list, ps, u0)

sol = simulate(odeprob)

plot(
    sol,
    xscale = :log10,
    yscale = :log10,
    xlims = (1.0e0 * SECONDS_PER_YEAR, 1.0e6 * SECONDS_PER_YEAR),
    idxs = [:CO, :H2, :H, Symbol("HCO+")],
)
