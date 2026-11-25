<!-- verified: 013551e 2026-11-20 -->
# CompReS
Julia pseudospectral solver for forced-dissipated 2D turbulence, compressible (low Mach) and incompressible, to study intermittency against Mach number. One source runs on CPU threads, NVIDIA and AMD GPUs.

## Status
Operative. The solver is complete; open work is in `PLAN.md` (Next, Ideas). The intermittency analysis waits for a specification by the user.

## Work plan
`PLAN.md` is the single list of work: **Next** (specified items), **Ideas** (unspecified) and **Done** (numbered, in order of completion). `design/` holds the design record; nothing else in the repo does. Follow these rules on every feature request:
1. Find the item in `PLAN.md`. If there is none, add it under Next (goal, ADRs to read, acceptance checks with their tier) and delete any Ideas entry it replaces. Ask the user for missing decisions before you implement.
2. If the specification exceeds about 10 lines, write it to `design/plans/<slug>.md` and link it from the item.
3. Read the ADRs the item lists. Write tests first, on small hand-checkable examples (single Fourier modes, `N = 8` or `16`). Make them pass on CPU and CUDA.
4. In the commit that completes the item, move it to Done as `N. Title — YYYY-MM-DD — ADRs — plan file`, with `N` one above the last Done number.
5. Record each new design decision as the next `design/adr/NNNN-slug.md` (header: `# Title`, then `Date: YYYY-MM-DD. Status: accepted.`). An ADR that changes an earlier one names it. Record each new term in `GLOSSARY.md`.
6. Add each idea that comes up and stays out of scope under Ideas (cost or defect, change, source).
Bug fixes and refactors change `PLAN.md` only when they change a design decision; that decision gets an ADR.

Tiers of an acceptance check: **blocking** passes in `Pkg.test` on CPU and CUDA (the default); **non-blocking** runs in `test/nonblocking/` and reports values; **user** needs a `pass` from the user in `test/acceptance/RESULTS.md`; **benchmark** meets its bound; **docs** builds with `docs/make.jl`.

Where design files go:
- `design/adr/` - design decisions; the docs site lists them (`docs/make.jl`).
- `design/plans/<slug>.md` - long specifications of items.
- `design/reviews/<target>/` - adversarial-review reports (`rN.md`) and their evidence scripts. This is the review location for `the-enemy`: for example, a review of `src/io.jl` goes to `design/reviews/src-io/`.
- `design/origin/` - the original prompt and the build roadmap (items 0 to 17 of Done). Read only; do not edit.

## Where to read first
- `PLAN.md`: open work and the record of implemented items.
- `GLOSSARY.md`: the project vocabulary. Use its terms in code, docstrings, commits and docs; its _Avoid_ lists name the words to replace.
- `docs/ns2d_equations.tex`: the binding equations, invariants, spectra and fluxes (ADR-0002). Cite equation names from it in docstrings.
- `design/adr/`: design decisions in chronological order. Read the ADRs an item lists before editing.

## Layout
Data flow of a run: `Simulation` (`src/simulation.jl`) builds a `Problem` (grid, model, forcing, linear rates; `src/nonlinear_incompressible.jl`), a `State` and a `Workspace` (`src/state.jl`); `run!` calls `step!` (`src/timestepping.jl`), which calls `update_forcing!`, `record_injection!` and, per stage, `rhs!` (`src/nonlinear_incompressible.jl`, `src/nonlinear_compressible.jl`) and `propagate!` or `linear_tendency!` (`src/linear.jl`); the writers of `src/io.jl` evaluate the diagnostics of `src/diagnostics/`.
- `src/diagnostics/` - `registry.jl` maps an output class and a name to a generic function `f(state, problem, workspace)` (ADR-0022); `incompressible.jl` and `compressible.jl` add the model methods and the default lists.
- `src/plotting.jl` - stubs and file readers of the plotting functions; `ext/` holds the Makie and Plots drawing code (ADR-0029).
- `benchmark/` - its own environment; `results/<machine>.toml` per machine (ADR-0030).
- `test/acceptance/` - user-inspected runs and `RESULTS.md`; `test/nonblocking/` - reported, not asserted.
- `design/` - design record: ADRs, plans, reviews, origin (see "Work plan").
- `runs/` - simulation output, ignored by git.

## Where to change X
- Add an equation of state: a struct in `src/eos.jl` with methods `pressure`, `enthalpy`, `internal_energy`, `sound_speed²`, `density_from_enthalpy`, `convert_float`.
- Add a forcing: a type in `src/forcing.jl` with `update_forcing!`, `forcing_state`, `set_forcing_state!`, `convert_float`; `ModeForcing` subtypes reuse `add_forcing!`; `check_forcing` rejects a model it cannot drive.
- Add a diagnostic: a method per model in `src/diagnostics/incompressible.jl` or `compressible.jl`, then `register_diagnostic!`; add the name to `INCOMPRESSIBLE_DEFAULTS` or `COMPRESSIBLE_DEFAULTS` to write it by default.
- Add a work array: one more field of `Workspace`, built in its constructor (`src/state.jl`).
- Add a Runge-Kutta scheme: a subtype of `RungeKuttaScheme` and its `butcher_table` in `src/timestepping.jl`.

## Commands
- All blocking tests: `julia --project -e 'using Pkg; Pkg.test()'` (CPU always, CUDA when `CUDA.functional()`).
- One test file: `julia --project=test test/runtests.jl <name>` (`<name>` is the file name without `.jl`; `test/Project.toml` loads CompReS from `..`).
- Acceptance animation (about 2 min on CUDA): `julia --project=test test/acceptance/crossed_kolmogorov.jl`.
- Non-blocking tier: `julia --project=test test/nonblocking/runtests.jl`.
- AMD tier (needs AMD hardware): `julia --project=test test/amdgpu/runtests.jl`.
- Docs: `julia --project=docs docs/make.jl`.
- Benchmarks (full sweep about 4 min): `julia --project=benchmark -t auto benchmark/run.jl [--key=list ...]` writes `benchmark/results/<machine>.toml`; `benchmark/predict.jl --N=... --model=... --device=... --T=... --t_end=... --dt=...` predicts a run (see `docs/src/tests_and_gpu.md`).
- A new file in `src/` needs an `@autodocs` line in one of `docs/src/api/*.md`; otherwise the docs check `checkdocs = :all` fails.

## Conventions
- Kernels use KernelAbstractions.jl only; no `CUDA.@cuda` or AMD-specific code outside `ext/`. Reason: one source for all devices (ADR-0001).
- Every function takes the float type from its arguments; literals are `T(…)`. Reason: Float32 runs must not promote (ADR-0006).
- Allocate all work arrays once, in a workspace; the time loop allocates nothing. Reason: speed on large grids.
- Every test is `run_tests(backend)` so CPU, CUDA and AMD share it (ADR-0011).
- Output variable names follow the `_inc` / `_comp` convention of ADR-0009.
- Plotting code lives only in `ext/` (ADR-0012).
- Between two equally fast designs, take the simpler one.
- Docstrings on every public symbol, in Documenter.jl format, with the equation from the specification in LaTeX.
- Write prose, comments and commit messages in plain technical English: one idea per sentence, active voice, one term per concept.

## Gotchas
- The real FFT stores half of the `kx` axis: Parseval sums double every mode except `kx = 0` and the Nyquist column.
- The `kx = 0` column stores both `(0, ky)` and `(0, −ky)`: any spectral field generated directly must pass through `enforce_hermitian!`, or the inverse FFT silently drops part of it.
- Fluxes vanish at `k_max` only for the incompressible model; the compressible ones carry an aliasing floor (specification, §4).
- In the compressible model the invariant is `Omega = ½⟨ω²/ρ⟩`, not `Z`.
- The RTX 5070 on the development machine runs Float64 at about 1/64 of its Float32 rate; GPU timings there do not predict cluster timings.
