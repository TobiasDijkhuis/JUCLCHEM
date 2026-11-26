"""
    initialize_species(row)

Initialize a species from a row of a CSV file."""
function initialize_species(row::CSV.Row)::SpeciesInfo
    return SpeciesInfo((
        name = row.NAME,
        mass = row.MASS,
        E_bind = row.BINDING_ENERGY,
        E_diff = 0.5 * row.BINDING_ENERGY,
        ΔH_F = row.ENTHALPY,
    ))
end

"""Mapping from the column name to the column data type."""
const SPECIES_COLUMN_TYPES = Dict(
    :NAME => String,
    :MASS => Float64,
    :ENTHALPY => Union{Missing, Float64},
    :BINDING_ENERGY => Union{Missing, Float64},
    :DIFFUSION_BARRIER => Union{Missing, Float64},
)

"""
    read_species_file(file)

Read a list of UCLCHEM species."""
function read_species_file(file::String)::Vector{SpeciesInfo}
    rows = CSV.Rows(file, types = SPECIES_COLUMN_TYPES)
    species = map(initialize_species, rows)
    return species
end

"""Mapping from the column name to the column data-type."""
const REACTION_COLUMN_TYPES = Dict(
    "REACTANT 1" => String,
    "REACTANT 2" => String,
    "REACTANT 3" => String,
    :ALPHA => Float64,
    :BETA => Float64,
    :GAMMA => Float64,
    :"PRODUCT 1" => Union{Missing, String},
    :"PRODUCT 2" => Union{Missing, String},
    :"PRODUCT 3" => Union{Missing, String},
    :"PRODUCT 4" => Union{Missing, String},
    :T_MIN => Union{Missing, Float64},
    :T_MAX => Union{Missing, Float64},
)

"""
    initialize_reaction(row)

Initialize a reaction from a row of a CSV file."""
function initialize_reaction(row)::ReactionInfo
    reaction_type = _get_reaction_type([
        row["REACTANT 1"],
        row["REACTANT 2"],
        row["REACTANT 3"],
    ])
    return ReactionInfo((
        r1 = row["REACTANT 1"],
        r2 = _is_real_species(row["REACTANT 2"]) ? row["REACTANT 2"] : missing,
        type = reaction_type,
        p1 = _is_real_species(row["PRODUCT 1"]) ? row["PRODUCT 1"] : missing,
        p2 = _is_real_species(row["PRODUCT 2"]) ? row["PRODUCT 2"] : missing,
        p3 = _is_real_species(row["PRODUCT 3"]) ? row["PRODUCT 3"] : missing,
        p4 = _is_real_species(row["PRODUCT 4"]) ? row["PRODUCT 4"] : missing,
        α = row["ALPHA"],
        β = row["BETA"],
        γ = row["GAMMA"],
        T_min = row["T_MIN"],
        T_max = row["T_MAX"],
    ))
end

"""
    read_reaction_file(file)

Read a list of UCLCHEM reactions."""
function read_reaction_file(file::String)::Vector{ReactionInfo}
    rows = CSV.Rows(file, types = REACTION_COLUMN_TYPES)
    reactions = map(initialize_reaction, rows)
    return reactions
end
