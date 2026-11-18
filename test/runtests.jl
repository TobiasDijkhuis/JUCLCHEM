using Test: @test, @testset, @test_throws
using JUCLCHEM:
    get_pure_reactants,
    ReactionInfo,
    TWOBODY,
    PHOTON,
    CRP,
    _get_reaction_type,
    create_catalyst_reaction,
    get_twobody_rate_constant,
    get_crp_photo_rate_constant,
    REACTION_TYPE_LOOKUP,
    get_symbolic_rate
using Catalyst
using Symbolics

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
    @test _get_reaction_type(["H2", "CRP", "NAN"]) == CRP
    @test _get_reaction_type(["C", "O", "NAN"]) == TWOBODY
    @test_throws ArgumentError _get_reaction_type(["C", "O"])
    @test_throws KeyError _get_reaction_type(["C", "O2", "unknown_kw"])
end

let
    function make_reaction_info(; r1, r2, p1, p2)::ReactionInfo
        return ReactionInfo((
            r1 = r1,
            r2 = r2,
            type = TWOBODY,
            p1 = p1,
            p2 = p2,
            p3 = missing,
            α = 0,
            β = 0,
            γ = 0,
        ))
    end
    ps = @parameters begin
        α[1:2]
        β[1:2]
        γ[1:2]
        T_gas
        ω_dust
        ζ
    end
    P = NamedTuple(Symbolics.getname(p) => p for p in ps)

    sp = @species C O CO C2
    species_dict = Dict(String(Symbolics.getname(spec)) => spec for spec in sp)

    @testset "create_catalyst_reaction $r1 + $r2 -> $p1 + $p2" for (
            r1,
            r2,
            p1,
            p2,
            expected_reactants,
            expected_reac_stoich,
            expected_products,
            expected_prod_stoich,
        ) in zip(
            ["C", "C", "C2"],
            ["O", "C", missing],
            ["CO", "C2", "C"],
            [missing, missing, "C"],
            [[C, O], [C], [C2]],
            [[1, 1], [2], [1]],
            [[CO], [C2], [C]],
            [[1], [1], [2]],
        )
        info = make_reaction_info(r1 = r1, r2 = r2, p1 = p1, p2 = p2)

        rx = create_catalyst_reaction(info, 1, species_dict, P)
        @test rx isa Reaction
        @test isequal(rx.substrates, expected_reactants)
        @test isequal(rx.products, expected_products)
        @test rx.substoich == expected_reac_stoich
        @test rx.prodstoich == expected_prod_stoich
    end
end

let
    idx = 1
    ps = @parameters begin
        α[1:2]
        β[1:2]
        γ[1:2]
        T_gas
        ω_dust
        ζ
    end
    P = NamedTuple(Symbolics.getname(p) => p for p in ps)
    EXPECTED_FUNCTIONS = Dict(
        "TWOBODY" => get_twobody_rate_constant(α[idx], β[idx], γ[idx], T_gas),
        "CRPHOT" => get_crp_photo_rate_constant(α[idx], β[idx], γ[idx], ω_dust, ζ, T_gas),
    )

    @testset "reaction to symbolic rate $reaction_type" for reaction_type in [
            "TWOBODY",
            "CRPHOT",
        ]
        type = REACTION_TYPE_LOOKUP[reaction_type]
        rate = get_symbolic_rate(type, idx, P)
        @test !isnothing(rate)
        @test isequal(rate, EXPECTED_FUNCTIONS[reaction_type])
    end
end

# @testset "reaction system" begin
#     rs = create_reaction_system(species_infos, reactions)
#     @test rs isa ReactionSystem
#     @test ModelingToolkitBase.iscomplete(rs)
#     @test length(Catalyst.species(rs)) == length(species_infos)    # no dropped species
#     @test Set(Symbol.(ModelingToolkitBase.parameters(rs))) ⊇ Set([:T_gas, :A_V, ...])
# end
