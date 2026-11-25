# Documentation site: page structure, API split and generated ADR index

Date: 2026-11-20. Status: accepted.

WP-17 fixes the structure of the Documenter site. Changes:

- **Pages.** Home, Getting started, Equations and conventions, Output format, Tests and GPU, API reference (one overview page and one page per section), ADR index. The pages cite the specification, the ADRs and the docstrings; they do not restate derivations.
- **API split.** One `@autodocs` block per page, selected by `Pages` (source files). One page of all docstrings is 169 KiB of HTML and trips the 100 KiB warning of Documenter. The sections follow the pipeline: grid and transforms, models and equations of state, initial conditions, forcing, time stepping, diagnostics, output and restart, simulation, plotting functions. A new source file needs a line in `docs/src/api/*.md`; the default `checkdocs = :all` fails the build if a docstring is on no page.
- **Plotting functions by name.** The plotting section selects `plot_field`, `plot_spectrum`, `plot_global` and `animate_snapshots` by name with a `Filter`, and every other section excludes these names. The selection does not depend on the source file that declares them (ADR-0012).
- **ADR index generated.** `docs/make.jl` lists `design/adr/NNNN-*.md`, reads the first `# ` line of each file as its title and links to the file on GitHub. A new ADR appears without an edit of the page.
- **Specification PDF.** `docs/make.jl` copies `docs/ns2d_equations.pdf` to `docs/src/` so the site serves it. The copy is in `.gitignore`.
