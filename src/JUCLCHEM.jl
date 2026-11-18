module JUCLCHEM
using Catalyst
using OrdinaryDiffEq
using SymbolicUtils
using DifferentialEquations

include("constants.jl")
include("network.jl")

function create_initial_conditions(species::AbstractVector)::AbstractDict
    abundances = Dict(spec => MIN_ABUNDANCE for spec in values(species))
    return abundances
end

struct PhysicsState
    T_gas::Real
    T_dust::Real
    A_V::Real
    n_H::Real
    radfield::Real
    ζ::Real
end

function get_pmap(reactions::AbstractVector{ReactionInfo}, ps::PhysicsState)
    return [
        :α => [r.α for r in reactions],
        :β => [r.β for r in reactions],
        :γ => [r.γ for r in reactions],
        :A_V => ps.A_V,
        :T_gas => ps.T_gas,
        :T_dust => ps.T_dust,
        :n_H => ps.n_H,
        :radfield => ps.radfield,
        :ζ => ps.ζ,
        :ω_dust => 0.5,
    ]
end

function create_odeproblem(
    system::ReactionSystem,
    species::AbstractVector{SpeciesInfo},
    reactions::AbstractVector{ReactionInfo},
    physics_state::PhysicsState,
    u0::Dict,
    tspan::Tuple{Real, Real} = (0, 1.0e6 * SECONDS_PER_YEAR),
)::ODEProblem
    pm = get_pmap(reactions, physics_state)
    ode_sys = complete(ode_model(system))

    pmap = Dict(getproperty(system, k) => v for (k, v) in pm)
    return ODEProblem(ode_sys, merge(u0, Dict(pmap)), tspan)
end

const DEFAULT_SAVE_T = collect(logrange(
    1.0e-7 * SECONDS_PER_YEAR,
    1.0e6 * SECONDS_PER_YEAR,
    100,
))

function simulate(
    oprob::ODEProblem,
    save_timepoints::AbstractVector{Float64} = DEFAULT_SAVE_T,
)::ODESolution
    sol = solve(oprob, saveat = save_timepoints)
    return sol
end
end # module JUCLCHEM
