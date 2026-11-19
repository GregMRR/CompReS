# Roadmap: work packages

Each work package (WP) is sized for one subagent. Read `CLAUDE.md`, `GLOSSARY.md` and the ADRs listed in the WP before starting. A WP is done when its acceptance tests pass on CPU and CUDA (`Pkg.test`), its public symbols have docstrings, and its status line below says `done` with the commit sha.

Notation: `WP-n ← WP-a, WP-b` means WP-n needs WP-a and WP-b done first. WPs with no mutual dependency can run in parallel.

## Target user API

A run script must read like this. Every WP keeps its interface compatible with it.

```julia
using NS2D

grid    = Grid(; Lx = 2π, Ly = 2π, Nx = 512, Ny = 512, T = Float64, backend = CPU())
model   = Compressible(; eos = Isothermal(cs = 10.0), μ = 1e-4, ζ = 0.0, κ = 0.0,
                         large_scale = LargeScaleDissipation(α = 0.1, m = 0))
forcing = OUForcing(band(4, 6); ε = 1.0, τ_f = 0.05, seed = 1)

sim = Simulation(grid, model; forcing, scheme = RK4(), integrating_factor = true, Δt = 1e-4,
                 outputs = (GLOBAL   = Output(period = 0.01),
                            SPECTRAL = Output(period = 0.1),
                            FIELDS   = Output(period = 1.0)),
                 outdir = "runs/example", max_walltime = 23.5 * 3600)

set_initial_condition!(sim, random_velocity(grid; spectrum = k -> k^4 * exp(-k^2 / 8),
                                            E_kin = 0.5, seed = 2))
run!(sim; t_end = 10.0)
# Restart: sim = Simulation(...same...); restart!(sim, "runs/example/fields/snap_000010.nc")
```

`Incompressible(; ν, large_scale)` replaces `Compressible(...)` for the incompressible model.

## Package layout (target)

```
Project.toml                 weak deps: CUDA, AMDGPU, GLMakie, CairoMakie, Plots
src/NS2D.jl                  module, exports, includes in dependency order
src/grid.jl                  WP-1
src/transforms.jl            WP-2
src/eos.jl                   WP-3
src/state.jl                 WP-4
src/initial_conditions.jl    WP-4
src/linear.jl                WP-5
src/nonlinear_incompressible.jl  WP-6
src/nonlinear_compressible.jl    WP-10
src/timestepping.jl          WP-7
src/diagnostics/registry.jl  WP-8
src/diagnostics/incompressible.jl  WP-8
src/diagnostics/compressible.jl    WP-11
src/forcing.jl               WP-9
src/io.jl                    WP-12
src/simulation.jl            WP-13
src/balanced_density.jl      WP-14
ext/NS2DCUDAExt.jl, ext/NS2DAMDGPUExt.jl     WP-0
ext/NS2DMakieExt.jl, ext/NS2DPlotsExt.jl     WP-15
test/runtests.jl, test/<module>.jl           blocking tier
test/nonblocking/runtests.jl                 non-blocking tier
test/acceptance/                             user-inspected animation
test/amdgpu/runtests.jl                      stored, not run
benchmark/                                   WP-16
docs/make.jl, docs/src/                      WP-0, WP-17
```

## Dependency graph

```
WP-0 → WP-1 → WP-2
WP-0 → WP-3
WP-2 → WP-4, WP-5, WP-9
WP-4 + WP-5 → WP-6
WP-5 + WP-6 → WP-7
WP-6 → WP-8
WP-3 + WP-6 → WP-10
WP-8 + WP-10 → WP-11
WP-10 → WP-14
WP-7 + WP-8 + WP-9 → WP-12
WP-7 + WP-12 → WP-13
WP-13 → WP-15, WP-16, WP-17, WP-18 (late)
```

## Work packages

### WP-0 Package skeleton and test harness
Status: todo. ADRs: 0001, 0011, 0012.
- `Project.toml` with deps (KernelAbstractions, AbstractFFTs, FFTW, NCDatasets, TimerOutputs, Random, TOML) and weak deps; empty extensions that load.
- `test/runtests.jl`: a `backends()` function that returns `[CPU()]` plus `CUDABackend()` when `CUDA.functional()`; a `@testset` loop that calls each test file's `run_tests(backend)`.
- `test/amdgpu/runtests.jl` calling the same functions with `ROCBackend()`.
- `test/nonblocking/runtests.jl` that reports results and always exits 0.
- `docs/make.jl` with Documenter; `.github/workflows/CI.yml` (CPU tests, docs build and deploy).
- Accept: `Pkg.test()` passes with one trivial test on CPU and CUDA; `julia --project=docs docs/make.jl` builds.

### WP-1 Grid
Status: todo. ← WP-0. ADRs: 0003, 0006.
- `Grid(; Lx, Ly, Nx, Ny, T, backend)`. Rejects `Lx/Nx ≠ Ly/Ny` and odd `N`.
- Device arrays for the real-FFT layout `(Nx÷2+1, Ny)`: `kx`, `ky`, `k2`, `truncation` (Bool, `|k| ≤ k_max`), `shell_index`; host vectors `shell_k`, `shell_count`.
- Accept: for `Lx = 2π, Nx = 8`, hand-listed `kx`, `ky`, `k_max`, retained modes and shell counts; rectangle `Lx = 4π, Ly = 2π` example; rejection of unequal spacing.

### WP-2 Transforms and spectral operators
Status: todo. ← WP-1. ADRs: 0001.
- Plans created once per grid (`plan_rfft`, `plan_irfft`) for CPU (FFTW threads), CUDA, AMD.
- `to_spectral!(f̂, f, grid)`, `to_physical!(f, f̂, grid)`, `truncate!(f̂, grid)`, `∂x!`, `∂y!`, `helmholtz_split!(ûs, ûc, û, grid)`, `project_solenoidal!`.
- Normalisation matches the specification: `f̂(k) = ⟨f e^{−ik·x}⟩`.
- Accept: round trip to round-off; derivatives of `sin(3x)cos(2y)` exact; projection idempotent and `k·û^s = 0`; Parseval `⟨fg⟩ = Σ f̂* ĝ` with the half-plane weights of the real FFT.

### WP-3 Equation of state
Status: todo. ← WP-0. ADRs: 0008.
- `Isothermal(cs)`, `Polytropic(K, γ)`; methods `pressure`, `enthalpy`, `sound_speed²`, `density_from_enthalpy`, `internal_energy`.
- Accept: analytic values at chosen `ρ`; `h'(ρ) = p'(ρ)/ρ` and `e'(ρ) = p/ρ²` by finite difference; `density_from_enthalpy ∘ enthalpy = id`; a kernel calling each method runs on the device.

### WP-4 State and initial conditions
Status: todo. ← WP-2. ADRs: 0002, 0008.
- `State` holds spectral fields in a `NamedTuple` (`(:ux, :uy)` or `(:ux, :uy, :rho)`), time `t`, step `n`.
- `random_velocity(grid; spectrum, E_kin, seed)` (solenoidal, truncated), `kolmogorov(grid; A, n)`, `crossed_kolmogorov(grid; A1, k1, A2, k2)`, `uniform_density`.
- Accept: `E_kin` hits its target; the field is solenoidal and inside the truncation; same seed gives the same field on CPU and CUDA.

### WP-5 Linear operators and integrating factor
Status: todo. ← WP-2. ADRs: 0004, 0005.
- Rate arrays for small-scale dissipation, large-scale dissipation (`k = 0` excluded) and density diffusivity.
- `linear_tendency!` (explicit use) and `propagator(Δt)` (IF use); compressible propagator diagonal in the Helmholtz basis.
- Interface leaves room for a later per-mode 2×2 acoustic propagator (ADR-0005, open question).
- Accept: single modes decay at `νk²`, `μk²`, `(μ+ζ)k²`, `α k^{−2m}`; `k = 0` untouched.

### WP-6 Incompressible nonlinear term
Status: todo. ← WP-4, WP-5. ADRs: 0002, 0003.
- Rotational form `ω ẑ × u` with projection; 5 FFTs per evaluation; result truncated.
- `rhs!(dû, state, model, forcing, workspace)` for the incompressible model.
- Accept: shear flow `sin(ky) x̂` gives zero to round-off; cellular flow `sin x sin y` gives zero; `Σ Re[û*·N̂] = 0` and `Σ Re[ω̂*·(u·∇ω)^] = 0` to round-off on a random field.

### WP-7 Time stepping
Status: todo. ← WP-5, WP-6. ADRs: 0005.
- Butcher tables `Euler`, `RK2`, `SSPRK3`, `RK4`; one generic `step!` with and without IF; preallocated stage buffers.
- Accept on scalar ODEs: convergence order 1, 2, 3, 4. **Blocking**: conservation test and Kolmogorov decay test (ADR-0011).

### WP-8 Diagnostics registry and incompressible diagnostics
Status: todo. ← WP-6. ADRs: 0009.
- `register_diagnostic!(class, name, f)`; per-class selection by name.
- Shell spectra and cumulative fluxes on device, reduced per shell; GLOBAL and SPECTRAL defaults of ADR-0009 for the incompressible model.
- Accept: spectrum of a single mode lands in the right shell with the right value; `Σ shells = E_kin`. **Blocking**: `Pi_kin(k_max) = Pi_Z(k_max) = 0` on a random projected field.

### WP-9 Forcing
Status: todo. ← WP-2. ADRs: 0007.
- `band(k_start, k_end)`, `modes`, `support`; `KolmogorovForcing`, `OUForcing`, `RandomPhaseForcing`; compressive fraction; `amp(k)`; per-step seeding by `(seed, step)`; `forcing_state` and `set_forcing_state!` for restart.
- Accept: hand-counted support on `N = 16`; empty or out-of-truncation support rejected; Kolmogorov field matches the formula; OU mean injection within 5% of `ε` over a long sample with white-noise velocity; phases change only at multiples of `τ_r`; solenoidal output when the compressive fraction is 0.

### WP-10 Compressible right-hand side
Status: todo. ← WP-3, WP-6. ADRs: 0002, 0004.
- Continuity, velocity equation with `∇h(ρ)`, explicit viscous remainder, `κ` term; all products in physical space, tendencies truncated.
- Accept: uniform state at rest stays at rest; a small acoustic wave has frequency `c_s k`; mass `⟨ρ⟩` conserved to round-off; with `ρ = 1` and solenoidal `u`, the solenoidal part of the tendency equals the incompressible one.

### WP-11 Compressible diagnostics
Status: todo. ← WP-8, WP-10. ADRs: 0009, 0011.
- `w = √ρ u`, `s = ω/√ρ`; spectra and fluxes with the `_inc` / `_comp` split; `E_int`, `Omega`, `p_theta`, `Ma`; `κ` terms in the budgets.
- Accept: with `ρ = 1`, compressible diagnostics equal the incompressible ones. **Non-blocking**: `Pi_kin` over all grid modes = 0 (rotational form); `Pi_Omega(k_max)` floor drops with `N`.

### WP-12 NetCDF I/O and restart
Status: todo. ← WP-7, WP-8, WP-9. ADRs: 0009.
- Writers for `global.nc`, `spectral.nc`, `fields/snap_NNNNNN.nc`; parameters as global attributes; `params.toml`; `read_snapshot`.
- Accept: written values read back exactly; a run of `2n` steps equals `n` steps + restart + `n` steps bit for bit, with OU forcing.

### WP-13 Simulation driver
Status: todo. ← WP-7, WP-12. ADRs: 0009, 0010.
- `Simulation`, `set_initial_condition!`, `run!`, `restart!`; per-class sampling periods; `run.log` (TimerOutputs, `profile` flag); `max_walltime`; NaN and CFL checks every `n_check` steps.
- Accept: a tiny forced run produces all files; a walltime of 0 writes one snapshot and stops; an injected NaN writes a snapshot and raises an error.

### WP-14 Balanced-density routine
Status: todo. ← WP-10. ADRs: 0013.
- `balanced_density(u, eos, grid)`.
- Accept: `θ` and `∂_tθ` vanish at `t = 0` to the expected order; at `Ma = 0.1` the emitted `E_kin_comp` is at least 10 times lower than with uniform `ρ`.

### WP-15 Visualisation extensions and acceptance animation
Status: todo. ← WP-13. ADRs: 0011, 0012.
- `plot_field`, `plot_spectrum`, `plot_global`, `animate_snapshots` in both extensions.
- Choose `ν` once for the crossed Kolmogorov run (ADR-0011) and fix it in `test/acceptance/crossed_kolmogorov.jl`. Produce the animation (vorticity plus spectrum tail) and add an entry `pending` to `test/acceptance/RESULTS.md`.
- Accept: the user records pass in `RESULTS.md`.

### WP-16 Benchmark suite
Status: todo. ← WP-13. ADRs: 0010.
- `benchmark/run.jl`: time per step and memory over `N` and `(model, scheme, device, T)`; fit; `benchmark/predict.jl` for a target run; TOML output.
- Accept: fit residual under 10% on this machine for `N = 256 … 4096`.

### WP-17 Documentation site
Status: todo. ← WP-13. ADRs: 0012.
- Documenter pages: getting started (the target script), equations (link to the specification), API from docstrings, output format, tests and GPU instructions, ADR index.
- Accept: docs build in CI with no missing docstrings.

### WP-18 Intermittency analysis (late)
Status: todo. ← WP-13. ADRs: 0012.
- Offline module on FIELDS snapshots: increments, structure functions, increment PDFs. Specify in detail before starting.
