module JUCLCHEM
using Catalyst
using OrdinaryDiffEq
using SymbolicUtils
using DifferentialEquations
using DiffEqCallbacks

include("constants.jl")
include("grain_reactions.jl")
include("gasphase_reactions.jl")
include("reaction.jl")
include("species.jl")
include("network.jl")
include("chemistry.jl")
include("io.jl")

const DEFAULT_ELEMENTAL_ABUNDANCES = Dict{String, Float64}(
    "H" => 1.0,
    "HE" => 0.1,
    "C" => 1.77e-4,
    "N" => 6.18e-5,
    "S" => 3.51e-6,
    "MG" => 2.256e-6,
    "SI" => 1.78e-6,
    "CL" => 3.39e-8,
    "P" => 7.78e-8,
)

function create_initial_conditions(species_symbols)
    @debug "Creating initial conditions for $(length(species_symbols)) species"

    by_name = Dict(String(Symbolics.getname(s)) => s for s in species_symbols)
    initial_conditions = Dict{Any, Float64}(
        spec => MIN_ABUNDANCE
        for spec in species_symbols
    )

    for (element, abundance) in pairs(DEFAULT_ELEMENTAL_ABUNDANCES)
        if !haskey(by_name, element)
            @debug "Skipping setting $element's initial abundance (not in network)."
            continue
        end
        initial_conditions[by_name[element]] = abundance
    end
    return initial_conditions
end

@kwdef struct PhysicsState
    T_gas
    T_dust
    A_V
    n_H
    radfield
    ζ
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
        :ζ => ps.ζ * ζ_0,
        :ω_dust => 0.5,
        :freeze_factor => 0.0,
    ]
end

function create_odeproblem(
    system::ReactionSystem,
    species::AbstractVector{SpeciesInfo},
    reactions::AbstractVector{ReactionInfo},
    physics_state::PhysicsState,
    u0,
    tspan = (0, 1.0e6 * SECONDS_PER_YEAR),
)::ODEProblem
    @debug "Creating ODEFunction from ReactionSystem"
    ode_sys = complete(ode_model(system))

    pm = get_pmap(species, reactions, physics_state)
    pmap = Dict(getproperty(system, k) => v for (k, v) in pm)

    @debug "Creating ODEProblem"
    return ODEProblem(
        ode_sys,
        merge(u0, pmap),
        tspan;
        sparse = true,
        checks = true,
        remove_conserved = true,
        jac = true,
    )
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
    sol = solve(
        oprob,
        Rodas5P();
        saveat = save_timepoints,
        abstol = 1.0e-25,
        reltol = 1.0e-7,
        callback = PositiveDomain(),
    )
    return sol
end
end # module JUCLCHEM
