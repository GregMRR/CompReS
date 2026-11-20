<!-- verified: 013551e 2026-11-20 -->
# CompReS
Julia pseudospectral solver for forced-dissipated 2D turbulence, compressible (low Mach) and incompressible, to study intermittency against Mach number. One source runs on CPU threads, NVIDIA and AMD GPUs.

## Status
The roadmap is implemented: WP-0 to WP-17 are done, except WP-15, which waits for the user's pass of the acceptance animation (`test/acceptance/RESULTS.md`). WP-18 (intermittency analysis) waits for a specification by the user; it is the first entry of `docs/BACKLOG.md`.

## Start new work
1. Take the task from the user, or an entry of `docs/BACKLOG.md` that the user names.
2. Add it to `docs/ROADMAP.md` as the next WP (dependencies, ADRs to read, acceptance tests) and delete the backlog entry.
3. Follow "Workflow for a work package" below.

## Where to read first
- `GLOSSARY.md`: the project vocabulary. Use its terms in code, docstrings, commits and docs; its _Avoid_ lists name the words to replace.
- `docs/ns2d_equations.tex`: the binding equations, invariants, spectra and fluxes (ADR-0002). Cite equation names from it in docstrings.
- `docs/adr/`: design decisions in chronological order. Read the ADRs a work package lists before editing.
- `docs/ROADMAP.md`: work packages, dependencies, interfaces, acceptance tests, target user API and package layout.
- `docs/BACKLOG.md`: possible future changes (performance, robustness) that no WP covers yet. Add every new idea there; move an entry into a WP when it becomes work.

## Layout
Data flow of a run: `Simulation` (`src/simulation.jl`) builds a `Problem` (grid, model, forcing, linear rates; `src/nonlinear_incompressible.jl`), a `State` and a `Workspace` (`src/state.jl`); `run!` calls `step!` (`src/timestepping.jl`), which calls `update_forcing!`, `record_injection!` and, per stage, `rhs!` (`src/nonlinear_incompressible.jl`, `src/nonlinear_compressible.jl`) and `propagate!` or `linear_tendency!` (`src/linear.jl`); the writers of `src/io.jl` evaluate the diagnostics of `src/diagnostics/`.
- `src/diagnostics/` - `registry.jl` maps an output class and a name to a generic function `f(state, problem, workspace)` (ADR-0022); `incompressible.jl` and `compressible.jl` add the model methods and the default lists.
- `src/plotting.jl` - stubs and file readers of the plotting functions; `ext/` holds the Makie and Plots drawing code (ADR-0029).
- `benchmark/` - its own environment; `results/<machine>.toml` per machine (ADR-0030).
- `test/acceptance/` - user-inspected runs and `RESULTS.md`; `test/nonblocking/` - reported, not asserted.

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

## Workflow for a work package
1. Read the WP in `docs/ROADMAP.md` and its listed ADRs.
2. Write tests first, on small hand-checkable examples (single Fourier modes, `N = 8` or `16`).
3. Make the tests pass on CPU and CUDA.
4. Set the WP status line to `done <sha>`.
5. Record any new design decision as the next `docs/adr/NNNN-slug.md`; record any new term in `GLOSSARY.md`.

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
