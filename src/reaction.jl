"""Reaction type"""
@enum ReactionType begin
    TWOBODY
    IONOPOL1
    IONOPOL2
    FREEZE
    BULKSWAP
    THERM
    LH
    LHDES
    ER
    ERDES
    H2FORM
    CRP
    CRPHOT
    DESCR
    DESOH2
    DEUVCR
    ED
    GAR
    PHOTON
    SURFSWAP
end

const ReactionInfo = @NamedTuple begin
    r1::String
    r2::Union{String, Missing}
    type::ReactionType
    p1::Union{String, Missing}
    p2::Union{String, Missing}
    p3::Union{String, Missing}
    α
    β
    γ
    T_min::Union{Missing, Any}
    T_max::Union{Missing, Any}
end

const REACTION_TYPE_LOOKUP = Dict(string(x) => x for x in instances(ReactionType))

"""
    _get_reaction_type(reactants)

Get the reaction type from a list of reactants.
"""
function _get_reaction_type(reactants::AbstractVector{String})::ReactionType
    if length(reactants) ≠ 3
        throw(ArgumentError(
            "reactants should have length 3, but had length $(length(reactants))",
        ))
    end
    if !isempty(reactants[end]) && reactants[end] ≠ "NAN"
        rt = get(REACTION_TYPE_LOOKUP, reactants[end], nothing)
        if isnothing(rt)
            throw(KeyError("Unknown reaction type $(reactants[end])"))
        end
        return rt
    end
    for name in (reactants[end], reactants[end - 1])
        rt = get(REACTION_TYPE_LOOKUP, name, nothing)
        isnothing(rt) || return rt
    end
    return TWOBODY
end

function collapse(names::AbstractVector{String}, sp::Dict)
    isempty(names) && return nothing, nothing
    u = unique(names)
    return [sp[n] for n in u], [count(==(n), names) for n in u]
end

function _is_real_species(field)::Bool
    return !ismissing(field) &&
        !isempty(field) &&
        field ≠ "NAN" &&
        !haskey(REACTION_TYPE_LOOKUP, field)
end

"""
    get_pure_reactants(r)


Get the pure reactants (without 'NAN' or reaction types)"""
function get_pure_reactants(r::ReactionInfo)::AbstractVector{String}
    return [x for x in (r.r1, r.r2) if _is_real_species(x)]
end

"""
    get_pure_products(r)


Get the pure products (without 'NAN')"""
function get_pure_products(r::ReactionInfo)::AbstractVector{String}
    return [x for x in (r.p1, r.p2, r.p3) if _is_real_species(x)]
end
