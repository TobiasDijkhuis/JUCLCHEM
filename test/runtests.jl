using Test: @test, @testset, @test_throws
using JUCLCHEM:
    get_pure_reactants,
    ReactionInfo,
    PhysicsState,
    SpeciesInfo,
    TWOBODY,
    PHOTON,
    BULKSWAP,
    LH,
    CRP,
    _get_reaction_type,
    create_catalyst_reaction,
    get_twobody_rate_constant,
    get_crp_photo_rate_constant,
    REACTION_TYPE_LOOKUP,
    get_symbolic_rate,
    create_reaction_system,
    read_species_file,
    read_reaction_file,
    create_odeproblem,
    create_initial_conditions,
    simulate,
    initialize_reaction
using Catalyst
using Symbolics
using SciMLBase
using SymbolicIndexingInterface

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

let
    function finish_row(row)
        return Dict(
            "REACTANT 1" => row[1],
            "REACTANT 2" => row[2],
            "REACTANT 3" => row[3],
            "PRODUCT 1" => row[4],
            "PRODUCT 2" => row[5],
            "PRODUCT 3" => "NAN",
            "PRODUCT 4" => "NAN",
            "ALPHA" => 0,
            "BETA" => 0,
            "GAMMA" => 0,
        )
    end
    function make_dummy_info(r1, r2, type, p1, p2, p3)
        return ReactionInfo((
            r1 = r1,
            r2 = r2,
            type = type,
            p1 = p1,
            p2 = p2,
            p3 = p3,
            α = 0,
            β = 0,
            γ = 0,
        ))
    end
    @testset "initialize_reaction $row" for (row, expected_info) in [
            (
                ["C", "O", "", "CO", "NAN"],
                make_dummy_info("C", "O", TWOBODY, "CO", missing, missing),
            ),
            (
                ["#C", "BULKSWAP", "", "@C", "NAN"],
                make_dummy_info("#C", missing, BULKSWAP, "@C", missing, missing),
            ),
            (
                ["#C", "#O", "LH", "#CO", "NAN"],
                make_dummy_info("#C", "#O", LH, "#CO", missing, missing),
            ),
            (
                ["#CH3OH", "#OH", "LH", "#CH3O", "#H2O"],
                make_dummy_info("#CH3OH", "#OH", LH, "#CH3O", "#H2O", missing),
            ),
        ]
        @test isequal(initialize_reaction(finish_row(row)), expected_info)
    end
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
        n_H
    end
    P = NamedTuple(Symbolics.getname(p) => p for p in ps)

    sp = @species C O CO C2
    species_dict = Dict(String(Symbolics.getname(spec)) => spec for spec in sp)
    masses_dict = Dict(String(Symbolics.getname(spec)) => 1.0 for spec in sp) # dummy mass
    species_indices = Dict(
        String(Symbolics.getname(spec)) => idx
        for (idx, spec) in enumerate(sp)
    )

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

        rx = create_catalyst_reaction(
            info,
            1,
            species_dict,
            P,
            masses_dict,
            species_indices,
        )
        @test rx isa Reaction
        @test isequal(rx.substrates, expected_reactants)
        @test isequal(rx.products, expected_products)
        @test rx.substoich == expected_reac_stoich
        @test rx.prodstoich == expected_prod_stoich
    end
end

let
    function make_reaction_info(reaction_type::String)::ReactionInfo
        type = REACTION_TYPE_LOOKUP[reaction_type]
        return ReactionInfo((
            r1 = "a",
            r2 = missing,
            type = type,
            p1 = missing,
            p2 = missing,
            p3 = missing,
            α = 0,
            β = 0,
            γ = 0,
        ))
    end
    idx = 1
    ps = @parameters begin
        α[1:2]
        β[1:2]
        γ[1:2]
        T_gas
        ω_dust
        ζ
        n_H
    end
    P = NamedTuple(Symbolics.getname(p) => p for p in ps)
    masses = Dict{String, Float64}()
    species_indices = Dict{String, Int64}("a" => 1)
    EXPECTED_FUNCTIONS = Dict(
        "TWOBODY" => get_twobody_rate_constant(α[idx], β[idx], γ[idx], T_gas) * n_H,
        "CRPHOT" => get_crp_photo_rate_constant(α[idx], β[idx], γ[idx], ω_dust, ζ, T_gas),
    )

    @testset "reaction to symbolic rate $reaction_type" for reaction_type in [
            "TWOBODY",
            "CRPHOT",
        ]
        info = make_reaction_info(reaction_type)
        rate = get_symbolic_rate(info, idx, P, masses, species_indices)
        @test !isnothing(rate)
        @test isequal(rate, EXPECTED_FUNCTIONS[reaction_type])
    end
end

@testset "reaction system" begin
    function make_dummy_species(name::String)::SpeciesInfo
        return SpeciesInfo((
            name = name,
            mass = 1.0,
            E_bind = 1.0,
            E_diff = 0.5,
            ΔH_F = 1.0,
        ))
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

    rs = create_reaction_system(species_infos, reaction_infos)
    @test rs isa ReactionSystem
    @test ModelingToolkitBase.iscomplete(rs)
    @test length(Catalyst.species(rs)) == length(species_infos) # no dropped species
    @test Set(Symbol.(ModelingToolkitBase.parameters(rs))) ⊇
        Set([:T_gas, :A_V, :T_dust, :n_H])
end

@testset "UCLCHEM network" begin
    species_path = joinpath(@__DIR__(), "network", "species.csv")
    reactions_path = joinpath(@__DIR__(), "network", "reactions.csv")
    @test isfile(species_path)
    @test isfile(reactions_path)
    species_list = read_species_file(species_path)
    reactions_list = read_reaction_file(reactions_path)
    @test !isempty(species_list)
    @test !isempty(reactions_list)

    system = create_reaction_system(species_list, reactions_list)
    @test ModelingToolkitBase.iscomplete(system)

    species_symbols = Catalyst.species(system)
    u0 = create_initial_conditions(species_symbols)

    ps = PhysicsState(
        T_gas = 10.0,
        T_dust = 10.0,
        ζ = 1.0,
        radfield = 1.0,
        A_V = 1.0,
        n_H = 1.0e4,
    )
    odeprob = create_odeproblem(system, species_list, reactions_list, ps, u0)

    @test isequal(Set(keys(u0)), Set(species_symbols))

    du = similar(odeprob.u0)
    odeprob.f(du, odeprob.u0, odeprob.p, 0.0)

    print(observables(system))

    @test all(isfinite, du)

    # any NaN/Inf parameter?
    for p in ModelingToolkitBase.parameters(system)
        v = odeprob.ps[p]
        @test all(isfinite, v)
    end

    rxs = Catalyst.reactions(system)
    rates = SymbolicIndexingInterface.getu(odeprob, [r.rate for r in rxs])(odeprob)
    bad = findall(k -> !isfinite(k) || k < 0, rates)
    for i in bad
        println(rxs[i], "   rate = ", rates[i])
    end

    sol = simulate(odeprob)
    @test SciMLBase.successful_retcode(sol)
    @test minimum(sol.u[end]) ≥ 0
end
