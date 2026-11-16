using Documenter
using JUCLCHEM


makedocs(sitename = "JUCLCHEM.jl", authors = "Tobias Dijkhuis")
deploy_docs(repo = "github.com/TobiasDijkhuis/JUCLCHEM.git", push_preview = true)
