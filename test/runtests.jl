using Test: @test, @testset
using JUCLCHEM: get_pure_reactants, ReactionInfo, TWOBODY, PHOTON, CRP, _get_reaction_type
using Catalyst

@testset "get_pure_reactants" begin
    function make_reaction_info(; r1::String, r2::Union{String, Missing})::ReactionInfo
        return ReactionInfo((
            r1 = r1,
            r2 = r2,
            type = TWOBODY,
            p1 = missing,
            p2 = missing,
            p3 = missing,
            α = 0,
            β = 0,
            γ = 0,
        ))
    end
    @test get_pure_reactants(make_reaction_info(r1 = "C", r2 = "O")) == ["C", "O"]
    @test get_pure_reactants(make_reaction_info(r1 = "CO", r2 = "PHOTON")) ==
        ["CO"] # keyword dropped
    @test get_pure_reactants(make_reaction_info(r1 = "H2", r2 = "CRP")) == ["H2"]
    @test get_pure_reactants(make_reaction_info(r1 = "C", r2 = missing)) ==
        ["C"] # missing dropped
end

@testset "_get_reaction_type" begin
    @test _get_reaction_type(["CO", "PHOTON", "NAN"]) == PHOTON
    @test _get_reaction_type(["H2", "CRP"]) == CRP
    @test _get_reaction_type(["C", "O"]) == TWOBODY
    @test _get_reaction_type(["C", "O", "NAN"]) == TWOBODY
    @test _get_reaction_type(["C", "unknown_kw"]) == TWOBODY
end

@testset "create_catalyst_reaction" begin
    rx = create_catalyst_reaction(info, 1, sp, P)
    @test rx isa Reaction
    @test isequal(rx.substrates, [C, O])
    @test rx.substoich == [1, 1]
    @test rx.prodstoich == [1]
end

let
    idx = 1
    @parameters begin
        α[1:2]
        β[1:2]
        γ[1:2]
        T_gas
        ω_dust
        ζ
    end
    P = (; α, β, γ, T_gas, ω_dust, ζ)
    EXPECTED_FUNCTIONS = Dict("TWOBODY" => get_twobody_rate_constant(
        α[idx],
        β[idx],
        γ[idx],
        T_gas,
    ))

    @testset "reaction to symbolic rate $reaction_type" for reaction_type in [
            "TWOBODY",
            "CRPHOT",
        ]
        type = JUCLCHEM.REACTION_TYPE_LOOKUP[reaction_type]
        rate = JUCLCHEM.get_symbolic_rate(type, idx, P)
        @test !isnothing(rate)
        @test isequal(rate, EXPECTED_FUNCTIONS[reaction_type])
    end
end
