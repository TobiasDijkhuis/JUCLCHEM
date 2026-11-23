const SpeciesInfo = @NamedTuple begin
    name::String
    mass
    E_bind
    ΔH_F
    E_diff
end

"""
    is_positive_ion(name)

Whether the species is a positive ion, determined by whether its name ends with "+"

"""
function is_positive_ion(name)
    return endswith(name, "+")
end
