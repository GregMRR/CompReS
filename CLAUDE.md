<!-- verified: initial 2026-11-19 -->
# CompReS
Julia pseudospectral solver for forced-dissipated 2D turbulence, compressible (low Mach) and incompressible, to study intermittency against Mach number. One source runs on CPU threads, NVIDIA and AMD GPUs.

## Status
Implementation in progress by work packages in `docs/ROADMAP.md`; each WP status line there says what is done.

## Where to read first
- `GLOSSARY.md`: the project vocabulary. Use its terms in code, docstrings, commits and docs; its _Avoid_ lists name the words to replace.
- `docs/ns2d_equations.tex`: the binding equations, invariants, spectra and fluxes (ADR-0002). Cite equation names from it in docstrings.
- `docs/adr/`: design decisions in chronological order. Read the ADRs a work package lists before editing.
- `docs/ROADMAP.md`: work packages, dependencies, interfaces, acceptance tests, target user API and package layout.

## Commands
- All blocking tests: `julia --project -e 'using Pkg; Pkg.test()'` (CPU always, CUDA when `CUDA.functional()`).
- One test file: `julia --project=test test/runtests.jl <name>` (WP-0 implements the name filter; `test/Project.toml` loads CompReS from `..`).
- Non-blocking tier: `julia --project=test test/nonblocking/runtests.jl`.
- AMD tier (needs AMD hardware): `julia --project=test test/amdgpu/runtests.jl`.
- Docs: `julia --project=docs docs/make.jl`.

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
