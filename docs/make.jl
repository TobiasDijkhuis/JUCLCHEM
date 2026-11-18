using Documenter
using JUCLCHEM

pages = Any["Home" => "index.md", "API" => "api.md"]

makedocs(
    sitename = "JUCLCHEM",
    authors = "Tobias Dijkhuis",
    pages = pages,
    clean = true,
    pagesonly = true,
)
deploydocs(repo = "github.com/TobiasDijkhuis/JUCLCHEM.git", push_preview = true)
