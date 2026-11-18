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
    r1::AbstractString
    r2::Union{AbstractString, Missing}
    type::ReactionType
    p1::Union{AbstractString, Missing}
    p2::Union{AbstractString, Missing}
    p3::Union{AbstractString, Missing}
    α::Real
    β::Real
    γ::Real
end

const REACTION_TYPE_LOOKUP = Dict(string(x) => x for x in instances(ReactionType))

function _get_reaction_type(reactants::Vector{String})::ReactionType
    if !isempty(reactants[end]) && reactants[end] ≠ "NAN"
        rt = get(REACTION_TYPE_LOOKUP, reactants[end], nothing)
        if isnothing(rt)
            error("Unknown reaction type $(reactants[end])")
        end
    end
    for name in (reactants[end], reactants[end - 1])
        rt = get(REACTION_TYPE_LOOKUP, name, nothing)
        isnothing(rt) || return rt
    end
    return TWOBODY
end

function get_pure_reactants(r::ReactionInfo)::AbstractVector{AbstractString}
    return [
        x
        for x in (r.r1, r.r2)
        if !ismissing(x) && x ≠ "NAN" && !haskey(REACTION_TYPE_LOOKUP, x)
    ]
end

function get_pure_products(r::ReactionInfo)::AbstractVector{AbstractString}
    return [
        x
        for x in (r.p1, r.p2, r.p3)
        if !ismissing(x) && x ≠ "NAN" && !haskey(REACTION_TYPE_LOOKUP, x)
    ]
end
