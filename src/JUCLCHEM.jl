module JUCLCHEM
using Catalyst: species, Reaction
using CSV: File
include("species.jl")
using .SpeciesModule: Species

"""Read the species from a file `file`."""
function read_species_list(file::AbstractString)::Vector{species}
    species_vector = Vector{Species}(nothing)
    for row in CSV.File(file)
        species = Species(row)
        push!(species_vector, species)
    end
    return
end

end # module JUCLCHEM
