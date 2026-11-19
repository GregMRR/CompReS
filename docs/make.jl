# Builds the documentation with Documenter.jl; CI also deploys it to GitHub Pages.
using Documenter
using NS2D

makedocs(;
    sitename = "NS2D",
    modules = [NS2D],
    remotes = nothing,
    format = Documenter.HTML(; edit_link = nothing, repolink = nothing),
    pages = ["Home" => "index.md"],
)

# GitHub Actions sets GITHUB_REPOSITORY; local builds skip the deployment.
if haskey(ENV, "GITHUB_REPOSITORY")
    deploydocs(; repo = "github.com/$(ENV["GITHUB_REPOSITORY"]).git", devbranch = "main")
end
