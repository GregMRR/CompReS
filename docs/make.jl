# Builds the documentation with Documenter.jl; CI also deploys it to GitHub Pages.
using Documenter
using CompReS

# Plotting functions form their own API section, wherever the source file defines them.
const PLOT_NAMES = (:plot_field, :plot_spectrum, :plot_global, :animate_snapshots)
isplot(x) = x isa Union{Function,Type} && nameof(x) in PLOT_NAMES
notplot(x) = !isplot(x)

const DOCS_DIR = @__DIR__
const ADR_URL = "https://github.com/GregMRR/CompReS/blob/main/docs/adr"

# Copy the specification next to the pages, so the site links to it.
for f in ("ns2d_equations.pdf",)
    cp(joinpath(DOCS_DIR, f), joinpath(DOCS_DIR, "src", f); force = true)
end

# Markdown list of every ADR file with its title; new files appear without editing the page.
function adr_index()
    io = IOBuffer()
    for file in sort(readdir(joinpath(DOCS_DIR, "adr")))
        m = match(r"^(\d{4})-.*\.md$", file)
        m === nothing && continue
        lines = readlines(joinpath(DOCS_DIR, "adr", file))
        i = findfirst(startswith("# "), lines)
        title = i === nothing ? file : strip(lines[i][3:end])
        println(io, "- ADR-", m[1], ": [", title, "](", ADR_URL, "/", file, ")")
    end
    return String(take!(io))
end
const ADR_INDEX = adr_index()

makedocs(;
    sitename = "CompReS",
    modules = [CompReS],
    repo = Remotes.GitHub("GregMRR", "CompReS"),
    format = Documenter.HTML(; edit_link = "main"),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Equations and conventions" => "equations.md",
        "Output format" => "output.md",
        "Tests and GPU" => "tests_and_gpu.md",
        "API reference" => [
            "Overview" => "api.md",
            "Grid and transforms" => "api/grid_transforms.md",
            "Models and equations of state" => "api/models.md",
            "Initial conditions" => "api/initial_conditions.md",
            "Forcing" => "api/forcing.md",
            "Time stepping" => "api/timestepping.md",
            "Diagnostics" => "api/diagnostics.md",
            "Output and restart" => "api/output.md",
            "Simulation" => "api/simulation.md",
            "Plotting functions" => "api/plotting.md",
        ],
        "ADR index" => "adr.md",
    ],
)

deploydocs(; repo = "github.com/GregMRR/CompReS.git", devbranch = "main")
