module JUCLCHEM
using Catalyst
using OrdinaryDiffEq
using SymbolicUtils
using DifferentialEquations

include("constants.jl")
include("grain_reactions.jl")
include("gasphase_reactions.jl")
include("reaction.jl")
include("species.jl")
include("network.jl")
include("chemistry.jl")
include("io.jl")

function create_initial_conditions(species::AbstractVector)::Dict{Symbol, Real}
    @debug "Creating initial conditions for $(length(species)) species"

    return Dict(spec => MIN_ABUNDANCE for spec in values(species))
end

struct PhysicsState
    T_gas::Real
    T_dust::Real
    A_V::Real
    n_H::Real
    radfield::Real
    ζ::Real
end

function get_pmap(
    species::AbstractVector{SpeciesInfo},
    reactions::AbstractVector{ReactionInfo},
    ps::PhysicsState,
)
    return [
        :α => [r.α for r in reactions],
        :β => [r.β for r in reactions],
        :γ => [r.γ for r in reactions],
        :E_bind => [s.E_bind for s in species],
        :E_diff => [s.E_diff for s in species],
        :A_V => ps.A_V,
        :T_gas => ps.T_gas,
        :T_dust => ps.T_dust,
        :n_H => ps.n_H,
        :radfield => ps.radfield,
        :ζ => ps.ζ,
        :ω_dust => 0.5,
        :freeze_factor => 1.0,
    ]
end

function create_odeproblem(
    system::ReactionSystem,
    species::AbstractVector{SpeciesInfo},
    reactions::AbstractVector{ReactionInfo},
    physics_state::PhysicsState,
    u0::Dict{Symbol, Real},
    tspan::Tuple{Real, Real} = (0, 1.0e6 * SECONDS_PER_YEAR),
)::ODEProblem
    @debug "Creating ODEFunction from ReactionSystem"
    ode_sys = complete(ode_model(system))

    pm = get_pmap(species, reactions, physics_state)
    pmap = Dict(getproperty(system, k) => v for (k, v) in pm)

    @debug "Creating ODEProblem"
    return ODEProblem(ode_sys, merge(u0, pmap), tspan; sparse = true)
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
