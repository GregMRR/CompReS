# Visualisation in package extensions; deferred scope

Date: 2026-11-19. Status: accepted.

All plotting lives in package extensions, so a simulation runs without loading or precompiling a plotting library. `NS2DMakieExt` triggers on GLMakie (interactive) or CairoMakie (headless cluster nodes) and uses only the Makie API. `NS2DPlotsExt` triggers on Plots, the light option. The core package declares the plotting function names; the extensions add the methods.

Deferred to a late stage: intermittency analysis (velocity increments, structure functions, increment PDFs). It runs offline on FIELDS snapshots in a separate analysis module, not inside the solver loop.

Tooling: docs use Documenter.jl and deploy to GitHub Pages; GitHub Actions runs the CPU tests and builds the docs. `NS2D.jl` is a placeholder package name.
