using DiffEqDevTools
using JUCLCHEM:
    create_reaction_system,
    create_odeproblem,
    create_initial_conditions,
    read_species_file,
    read_reaction_file,
    PhysicsState
using Catalyst
using OrdinaryDiffEq: Rodas5P
using OrdinaryDiffEqBDF: QNDF, FBDF, NordsieckBDF
using OrdinaryDiffEqRosenbrock: Rodas4P
using Sundials: CVODE_BDF
using ODEInterfaceDiffEq: radau, ddebdf
using OrdinaryDiffEqExtrapolation:
    ImplicitEulerExtrapolation,
    ImplicitEulerBarycentricExtrapolation
using Plots

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
    u0,
    jac = true,
    sparse = true,
    remove_conserved = true,
)

# Reference ("true") solution computed with a tight tolerance.
test_sol = TestSolution(solve(
    odeproblem,
    radau();
    abstol = 1.0e-25,
    reltol = 1.0e-14,
    maxiters = 10^8,
))
println("Reference solution finished.")

# Tolerance sweeps (abstol and reltol must have the same length)
abstols = 10.0.^(-20:2:-12)
reltols = 10.0.^(-10:1:-6)

setups = [
    Dict(:alg => FBDF()),
    Dict(:alg => QNDF()),
    Dict(:alg => NordsieckBDF()),
    Dict(:alg => Rodas4P()),
    Dict(:alg => CVODE_BDF()),
    Dict(:alg => ddebdf()),
    # # Dict(:alg => Rodas4()),
    Dict(:alg => Rodas5P()),
    # # Dict(:alg => rodas()),
    Dict(:alg => radau()),
    # Dict(:alg => lsoda()),
    Dict(:alg => ImplicitEulerExtrapolation(
        min_order = 5,
        init_order = 3,
        threading = false,
    )),
    Dict(:alg => ImplicitEulerBarycentricExtrapolation(min_order = 5, threading = false)),
]

wp = WorkPrecisionSet(
    odeproblem,
    abstols,
    reltols,
    setups;
    appxsol = test_sol,
    save_everystep = false,
    maxiters = 10^8,
    numruns = 5,
    error_estimate = :final,
)

plt = plot(wp; title = "JUCLCHEM work-precision", legend = :topright)
savefig(plt, joinpath((@__DIR__), "workprecision.png"))
display(plt)
