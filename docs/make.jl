# Builds the documentation with Documenter.jl; CI also deploys it to GitHub Pages.
using Documenter
using CompReS

makedocs(;
    sitename = "CompReS",
    modules = [CompReS],
    repo = Remotes.GitHub("GregMRR", "CompReS"),
    pages = ["Home" => "index.md", "API reference" => "api.md"],
)

deploydocs(; repo = "github.com/GregMRR/CompReS.git", devbranch = "main")
