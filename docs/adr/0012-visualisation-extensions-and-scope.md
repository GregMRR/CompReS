# Visualisation in package extensions; deferred scope

Date: 2026-11-19. Status: accepted.

All plotting lives in package extensions, so a simulation runs without loading or precompiling a plotting library. `NS2DMakieExt` triggers on Makie, which GLMakie (interactive) and CairoMakie (headless cluster nodes) both load, and uses only the Makie API. An extension loads only when all its triggers are loaded, so it cannot trigger on "GLMakie or CairoMakie" (amended by ADR-0014). `NS2DPlotsExt` triggers on Plots, the light option. The core package declares the plotting function names; the extensions add the methods.

Deferred to a late stage: intermittency analysis (velocity increments, structure functions, increment PDFs). It runs offline on FIELDS snapshots in a separate analysis module, not inside the solver loop.

Tooling: docs use Documenter.jl and deploy to GitHub Pages; GitHub Actions runs the CPU tests and builds the docs. `NS2D.jl` is a placeholder package name.
