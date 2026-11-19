# Roadmap: work packages

Each work package (WP) is sized for one subagent. Read `CLAUDE.md`, `GLOSSARY.md` and the ADRs listed in the WP before starting. A WP is done when all three conditions hold:
1. Every acceptance item is met in its tier: **blocking** items pass in `Pkg.test` on CPU and CUDA; **non-blocking** items run in `test/nonblocking/` and report their values; **user** items have a `pass` entry from the user in `test/acceptance/RESULTS.md`; **benchmark** items meet their stated bound; **docs** items build with `julia --project=docs docs/make.jl` without error. An item with no tier label is blocking.
2. Its public symbols have docstrings.
3. Its status line below says `done` with the commit sha.

Notation: `WP-n ← WP-a, WP-b` means WP-n needs WP-a and WP-b done first. WPs with no mutual dependency can run in parallel.

## Target user API

A run script must read like this. Every WP keeps its interface compatible with it.

```julia
using CompReS

grid    = Grid(; Lx = 2π, Ly = 2π, Nx = 512, Ny = 512, T = Float64, backend = CPU())
model   = Compressible(; eos = Isothermal(cs = 10.0), μ = 1e-4, ζ = 0.0, κ = 0.0,
                         large_scale = LargeScaleDissipation(α = 0.1, m = 0))
forcing = OUForcing(grid, band(4, 6); ε = 1.0, τ_f = 0.05, seed = 1)   # band in physical |k|

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

Parameters are written as plain Julia numbers. `Simulation` converts every parameter struct (model, EOS, forcing) to the float type `T` of the grid, so a Float32 run never promotes to Float64 (ADR-0006). The initial density is uniform, `ρ = 1`, unless the initial condition sets `rho`.

## Package layout (target)

```
Project.toml                 weak deps: CUDA, AMDGPU, Makie, Plots
test/Project.toml            test environment: CompReS via [sources] path "..", Test, CUDA, AMDGPU
docs/Project.toml            docs environment: CompReS via [sources] path "..", Documenter
src/CompReS.jl                  module, exports, includes in dependency order
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
ext/CompReSCUDAExt.jl, ext/CompReSAMDGPUExt.jl     WP-0
ext/CompReSMakieExt.jl (trigger: Makie), ext/CompReSPlotsExt.jl (trigger: Plots)   WP-15
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
WP-2 → WP-4
WP-2 + WP-3 → WP-5
WP-4 + WP-5 → WP-6
WP-6 → WP-9
WP-6 + WP-9 → WP-8
WP-5 + WP-6 + WP-8 → WP-7
WP-3 + WP-6 → WP-10
WP-8 + WP-10 → WP-11
WP-7 + WP-11 → WP-14
WP-7 + WP-8 + WP-9 → WP-12
WP-7 + WP-12 → WP-13
WP-13 → WP-15, WP-16, WP-17, WP-18 (late)
```

## Work packages

### WP-0 Package skeleton and test harness
Status: done a7daa08. ADRs: 0001, 0011, 0012.
- `Project.toml` with deps (KernelAbstractions, AbstractFFTs, FFTW, NCDatasets, TimerOutputs, Random, TOML) and weak deps; empty extensions that load. Tests use their own environment `test/Project.toml` (CompReS through `[sources]` with `path = ".."`, plus Test, CUDA, AMDGPU), so `julia --project=test test/<file>.jl` works outside `Pkg.test`.
- `test/runtests.jl`: a `backends()` function that returns `[CPU()]` plus `CUDABackend()` when `CUDA.functional()`; a `@testset` loop that calls each test file's `run_tests(backend)`. A helper `is_cpu(backend)` lets a test function run a CPU-only item (ADR-0011 exception) and its device counterpart on the other backends.
- `test/amdgpu/runtests.jl` calling the same functions with `ROCBackend()`.
- `test/nonblocking/runtests.jl` that reports results and always exits 0.
- `docs/make.jl` and `docs/Project.toml` (CompReS through `[sources]`, Documenter); `.github/workflows/CI.yml` (CPU tests, docs build and deploy).
- Accept: `Pkg.test()` and `julia --project=test test/runtests.jl` pass with one trivial test on CPU and CUDA. Accept (docs): the docs build.

### WP-1 Grid
Status: done 2a68a48. ← WP-0. ADRs: 0003, 0006.
- `Grid(; Lx, Ly, Nx, Ny, T, backend)`. Rejects odd `N` and unequal spacing: `|Lx/Nx − Ly/Ny| > 1e-12·Lx/Nx`, evaluated in Float64 for every `T`.
- Galerkin truncation: keep `|k| < k_c` (strict), `k_c = (Nx/3)(2π/Lx)`. The strict bound keeps the disc alias-free also when 3 divides `N` (the self-triad at `|k| = N/3` would fold back inside). Evaluate it on the integer mode indices `(i, j)`: keep the mode iff `9(i²Ny² + j²Nx²) < Nx²Ny²` (Int64), never on float wavenumbers. `k_max` is the largest retained `|k|`.
- Device arrays for the real-FFT layout `(Nx÷2+1, Ny)`: `kx`, `ky`, `k2`, `truncation` (Bool), `shell_index`; host vectors `shell_k` (centres `nΔk`), `shell_edge` (upper edges `(n+½)Δk`), `shell_count` (retained modes) and `shell_count_all` (all grid modes, for the compressible spectra of ADR-0020); `shell_index` covers every grid mode, and the shells extend to the grid corner. Assign shells on integer indices too: mode `(i, j)` is in the smallest shell `n` with `4·Nm²(i²Ny² + j²Nx²) < (2n+1)²Nx²Ny²`, `Nm = max(Nx, Ny)` (Int128, because the products exceed Int64 for `N ≳ 1300`). A mode on an edge then goes to the upper shell, as the specification's half-open shells require.
- Accept: for `Lx = 2π, Nx = 8`, hand-listed `kx`, `ky`, `k_max`, retained modes and shell counts; `N = 12` has no retained mode at `|k| = 4`; on `Lx = 5·2π, Ly = 4·2π, Nx = 60, Ny = 48`, mode `(6, −14)` lies on an edge and goes to shell 19; rectangle `Lx = 4π, Ly = 2π` example; rejection of unequal spacing; acceptance of `Lx = 0.3, Nx = 30, Ly = 0.1, Ny = 10` with `T = Float64` and `Float32`; `Lx = 0.3, Nx = 30, Ly = 0.9, Ny = 90` does not keep mode `(0, 30)`.

### WP-2 Transforms and spectral operators
Status: done 7b469ea. ← WP-1. ADRs: 0001.
- Plans created once per grid (`plan_rfft`, `plan_irfft`) for CPU (FFTW threads), CUDA, AMD.
- `Transforms(grid)` holds the plans and the derivative wavenumbers. `to_spectral!(f̂, f, tr)`, `to_physical!(f, f̂, tr)` (never changes `f̂`), `∂x!(df̂, f̂, tr)`, `∂y!`, `helmholtz_split!(ûs, ûc, û, tr)`, `project_solenoidal!(û, tr)` take the `Transforms`; `truncate!(f̂, grid)` and `enforce_hermitian!(f̂, grid)` take the grid. Vector fields are pairs of arrays. The forward normalisation is a separate pass; later WPs may fuse it into `truncate!` where a truncation follows.
- Normalisation matches the specification: `f̂(k) = ⟨f e^{−ik·x}⟩`.
- Derivative multipliers `ikx`, `iky` are zero on the Nyquist column and row, so the discrete derivative is skew-adjoint and keeps Hermitian symmetry (needed by the grid identity of WP-11).
- `enforce_hermitian!(f̂, grid)`: on the `kx = 0` column copy from the source half `ky > 0` to `f̂(0, −ky) = conj(f̂(0, ky))` (one thread per source mode, so no race on GPU), make `f̂(0, 0)` real, and set the Nyquist column and row to zero. Every routine that generates spectral values directly (initial conditions, forcing) calls it, because the real-FFT layout stores both `(0, ky)` and `(0, −ky)`.
- Accept: round trip to round-off; derivatives of `sin(3x)cos(2y)` exact; the multiplier arrays are zero on the Nyquist column and row (checked on the arrays, and on the spectral output of `∂x!` for `cos(Nx x/2)`, not after an inverse FFT, which hides the defect); `enforce_hermitian!` output, from random untruncated input, survives an `irfft`/`rfft` round trip unchanged to round-off; `Σ f·∂x g = −Σ g·∂x f` on random grid fields (not truncated) to round-off; projection idempotent, `k·û^s = 0`, `k = 0` mode finite; Parseval `⟨fg⟩ = Σ f̂* ĝ` with the half-plane weights of the real FFT.

### WP-3 Equation of state
Status: done c38fdc4. ← WP-0. ADRs: 0008.
- `Isothermal(; cs)`, `Polytropic(; K, γ)`; the five methods of ADR-0008: `pressure`, `enthalpy`, `internal_energy`, `sound_speed²`, `density_from_enthalpy`; `convert_float(T, eos)`.
- Accept: analytic values at chosen `ρ`; `h'(ρ) = p'(ρ)/ρ` and `e'(ρ) = p/ρ²` by finite difference; `density_from_enthalpy ∘ enthalpy = id`; a kernel calling each method runs on the device.

### WP-4 State and initial conditions
Status: done 179a778. ← WP-2. ADRs: 0002, 0008.
- `State` holds spectral fields in a `NamedTuple` (`(:ux, :uy)` or `(:ux, :uy, :rho)`), time `t`, step `n`.
- `Workspace(grid, state)`: the `Transforms` and all preallocated work arrays. Later WPs add their buffers to it as one typed field each, built in the constructor; no other module allocates in the time loop.
- `random_velocity(grid; spectrum, E_kin, seed)` (solenoidal, truncated; `spectrum(k)` is the target shell spectrum `E_kin_inc(k)` up to a factor that `E_kin` fixes; random phases), `kolmogorov(grid; A, n)`, `crossed_kolmogorov(grid; A1, k1, A2, k2)`, `uniform_density`.
- Accept: `E_kin` hits its target, computed both spectrally and as the grid mean of the physical field; the field is solenoidal, inside the truncation and Hermitian on the `kx = 0` column (`û(0, −ky) == conj(û(0, ky))`); same seed gives the same field on CPU and CUDA.

### WP-5 Linear operators and integrating factor
Status: done d96cb38. ← WP-2, WP-3. ADRs: 0004, 0005.
- Model types `Incompressible(; ν, large_scale, nonlinear = true)`, `Compressible(; eos, μ, ζ, κ, large_scale, nonlinear = true)`; `nonlinear = false` is a test switch that drops every term except the linear ones and the forcing (incompressible: the advection; compressible: advection, pressure, viscous remainder and the continuity flux, so `ρ` stays fixed) (EOS field is parametric; the tests use the WP-3 EOS types), `LargeScaleDissipation(; α, m)`; `convert_float(T, x)` for each parameter struct.
- Rate arrays for small-scale dissipation, large-scale dissipation (`k = 0` excluded) and density diffusivity.
- `linear_rates(grid, model)`; `linear_tendency!(dû, û, rates)` (explicit use) and `propagate!(out, û, rates, τ)` (IF use, computes `exp(−στ)` per mode for any `τ = cΔt`); compressible propagator diagonal in the Helmholtz basis.
- Interface leaves room for a later per-mode 2×2 acoustic propagator (ADR-0005, open question).
- Accept: single modes decay at `νk²`, `μk²`, `(μ+ζ)k²`, `α k^{−2m}`, `κk²` (density); `k = 0` untouched; `convert_float(Float32, model)` holds only Float32.

### WP-6 Incompressible nonlinear term
Status: done 923fc45. ← WP-4, WP-5. ADRs: 0002, 0003.
- Lamb form `(u·∇)u = ∇(½|u|²) + ω ẑ × u`: compute `ω ẑ × u` and project (the projection removes the gradient); 5 FFTs per evaluation; result truncated.
- `Problem(grid, model, forcing)`: bundles grid, model, forcing and the linear rates of WP-5; WP-7 advances it, and WP-8 and WP-11 add `record_injection!` methods for it.
- `rhs!(dû, state, problem, workspace)` for the incompressible model; it honours `nonlinear = false` (WP-5) by skipping the advection. `rhs!` returns the nonlinear term plus the forcing only; `step!` (WP-7) adds the linear terms with `linear_tendency!` or applies them with `propagate!`. It calls the forcing hook `add_forcing!(dû, forcing, state, workspace)`; WP-6 defines `NoForcing` and its no-op method, WP-9 adds the others.
- Accept: shear flow `sin(ky) x̂` gives zero to round-off; cellular flow with streamfunction `ψ = sin x sin y` (`u = (sin x cos y, −cos x sin y)`) gives zero; `Σ Re[û*·N̂] = 0` and `Σ Re[ω̂*·(u·∇ω)^] = 0` to round-off on a random field.

### WP-7 Time stepping
Status: done 16a063f. ← WP-5, WP-6, WP-8. ADRs: 0005, 0011, 0023.
- Butcher tables `Euler`, `RK2`, `SSPRK3`, `RK4`; one generic `step!` with and without IF; preallocated stage buffers. `step!(state, problem, workspace, Δt)` is generic: `problem` bundles everything the right-hand side needs (model, forcing, linear rates, grid) and provides `rhs!(dû, state, problem, workspace)`, so the same code advances a full grid and a one-mode test problem.
- `step!` calls `update_forcing!` (WP-9) once, before the first stage; the forcing stays constant inside the stages (ADR-0005). `update_forcing!` writes the new `f̂` into the workspace. Right after it, `step!` calls `record_injection!(problem, state, workspace, Δt)`, which reads `f̂` from the workspace and stores the injection rates of WP-8 there. Its default method does nothing, so the one-mode test problem needs no change. The stored rates start as NaN. `run!` writes a sample of every output class at the start time `t₀` (start or restart), so that sample holds NaN for `inj_*`, written as a plain value (NCDatasets reads it back as NaN).
- Accept: through the production `step!`, on a one-mode problem `y' = −2y + y²` (`L = −2`, nonlinear part `y²`), convergence orders 1, 2, 3, 4 **both without and with the IF**. Conservation test and Kolmogorov decay test (ADR-0011), with `E_kin` and `Z` from WP-8. Injection test: an incompressible run with forcing only (nonlinear term, `ν` and `α` off, OU forcing) gives `|inj_E_kin − (E_kin(t_{n+1}) − E_kin(t_n))/Δt| ≤ 10·eps(T)·E_kin(t_{n+1})/Δt`, and the same for `Z`, for `T = Float64` and `Float32`, over 1000 steps, for every scheme (the scale is `E_kin/Δt` because the energy difference cancels) (with forcing only, every RK scheme gives `u_{n+1} = u_n + Δt f_n`).

### WP-8 Diagnostics registry and incompressible diagnostics
Status: done 90b6d1b. ← WP-6, WP-9. ADRs: 0009, 0022.
- `register_diagnostic!(class, name, f)`; per-class selection by name.
- Shell spectra and cumulative fluxes on device, reduced per shell; GLOBAL and SPECTRAL defaults of ADR-0009 for the incompressible model, including the injection rates (from the forcing of WP-9).
- Flux coordinate: `Pi(K_n)` sums all modes in shells `0…n`, that is `|q| < K_n` with `K_n = (n+½)Δk` (`shell_edge`). It equals the specification's `Π(K)` for every `K` between the largest `|q|` in shells `0…n` and `K_n`, excluding `K_n` itself (on a rectangle a mode can lie exactly on `K_n`; it belongs to shell `n+1`). It is reported at `K_n`; spectra are reported at the shell centres.
- Injection rate: for forcing held constant over a step, `inj_E_kin = Σ w_k (Re[û*·f̂] + ½Δt|f̂|²)` per unit time, with `û` and `f̂` both taken at the start of the step, and `w_k` the half-plane weight (1 on the `kx = 0` column and the Nyquist column, 2 elsewhere); `inj_Z` uses weights `w_k|k|²`. The `½Δt|f̂|²` term removes the `O(Δt/τ_f)` bias of the plain product. A diagnostic sees the state after the step, so WP-8 provides the `record_injection!` method for the incompressible model, which `step!` (WP-7) calls at the start of each step; the `inj_*` diagnostics read the stored values.
- Accept: spectrum of a single mode lands in the right shell with the right value; `Σ shells = E_kin`. **Blocking**: `Pi_kin(k_max) = Pi_Z(k_max) = 0` on a random projected field.

### WP-9 Forcing
Status: done 32d0299. ← WP-6. ADRs: 0007.
- `band(k_start, k_end)` in physical wavenumber units (`|k|`, not mode index), `modes` (integer mode indices), `support(pred)`; constructors take the grid first: `KolmogorovForcing(grid; F, n)`, `OUForcing(grid, support; ε, τ_f, seed)`, `RandomPhaseForcing(grid, support; f_rms, τ_r, seed)`; keyword `compressive_fraction` (default 0); `amp(k)`; per-step seeding by `(seed, step)`; `forcing_state` and `set_forcing_state!` for restart.
- `update_forcing!(forcing, state, workspace, Δt)`: advance the forcing state once per step (OU update, phase redraw) and write the new `f̂` into the workspace; called by `step!` (WP-7).
- `add_forcing!` methods for each type: add `f̂` to the tendency. The incompressible model always projects `f̂`; a nonzero compressive fraction is allowed only with the compressible model, and the constructor or `Simulation` rejects it otherwise.
- Generated forcing values pass through `enforce_hermitian!` (WP-2). `convert_float(T, forcing)`.
- Accept: hand-counted support on `N = 16`; `band` on `Lx = 4π, Ly = 2π` selects by physical `|k|`; empty or out-of-truncation support rejected at construction; Kolmogorov field matches the formula.
- Accept (OU normalisation, CPU only): `N = 32`, `band(4, 6)`, `Δt = 0.01τ_f`, seed fixed, `M = 4000` realisations. Drive the forced modes with `u' = f` from `u = 0` (no other term). At the start of a step (`u_n` with `f_n`, after `update_forcing!`), measure the injection as the grid mean `⟨u·f⟩` of the physical fields plus `½Δt⟨|f|²⟩`. Sample it once per `τ_f` over `t ∈ [20τ_f, 40τ_f]` (samples closer than `τ_f` add no information). Average the samples of each realisation first; the `M` realisation means are independent. Pass when `|mean − ε| ≤ 3.3·std/√M` (99.9% interval; the fixed seed makes the test deterministic). A factor-2 variance error and a non-Hermitian generator both lie outside this interval for every seed. Budget: under 60 s on CPU. On CUDA, the test checks instead that 100 steps with the same seed give the CPU forcing sequence to round-off.
- Accept: every generated `f̂` is Hermitian on the `kx = 0` column.
- Accept: phases change only at multiples of `τ_r`; solenoidal output when the compressive fraction is 0.

### WP-10 Compressible right-hand side
Status: done 32d0299. ← WP-3, WP-6. ADRs: 0002, 0004.
- Continuity, velocity equation with `∇h(ρ)`, explicit viscous remainder, `κ` term; all products in physical space, tendencies truncated. `rhs!` calls `add_forcing!` (same hook as WP-6), so the compressive forcing fraction applies here. With `nonlinear = false` it keeps only the linear terms and the forcing; with `κ = 0` the density then stays fixed.
- Accept: uniform state at rest stays at rest; a small acoustic wave has frequency `c_s k`; mass `⟨ρ⟩` conserved to round-off; with `ρ = 1` and solenoidal `u`, the solenoidal part of the tendency equals the incompressible one.

### WP-11 Compressible diagnostics
Status: done 865c126. ← WP-8, WP-10. ADRs: 0009, 0011, 0025.
- `w = √ρ u`, `s = ω/√ρ`. `√ρ` and `1/ρ` are not band-limited, so `ŵ` and `ŝ` are not truncated: spectra and fluxes sum over all grid modes (shells up to the grid corner), and GLOBAL `E_kin`, `Omega`, `E_int` are grid means of the physical fields. By Parseval the two agree. Spectra and fluxes carry the `_inc` / `_comp` split; `E_int`, `Omega`, `p_theta`, `Ma`; `κ` terms in the budgets.
- `record_injection!` for the compressible model: rates from the global budgets of the specification, `⟨ρu·f⟩ + ½Δt⟨ρ|f|²⟩` for `E_kin` and `⟨q curl f⟩ + ½Δt⟨|curl f|²/ρ⟩` for `Omega`, at the start of the step.
- Accept: with `ρ = 1`, compressible diagnostics, including the injection rates, equal the incompressible ones.
- Accept: forcing-only compressible run (`nonlinear = false`, `μ = ζ = κ = α = 0`, `OUForcing(grid, band(4, 6); ε = 1, τ_f = 0.05, seed = 1, compressive_fraction = 0.5)`, `Δt = 1e-3`), `N = 32`, initial `u` from `random_velocity` (not at rest), fixed nonuniform `ρ = 1 + 0.1 cos(x + y)`, 100 steps, `T = Float64` and `Float32`: `|inj_E_kin − ΔE_kin/Δt| ≤ 10·eps(T)·E_kin(t_{n+1})/Δt`, and the same for `Omega`. A multi-mode band is required: `ρ` couples only modes that differ by `±(1, 1)`, so a single forced mode hides a formula that drops `ρ`. On this set-up such formulas fail by several percent (the exact value depends on the seed), far above the tolerance. **Non-blocking**: `Pi_kin` over all grid modes = 0 (Lamb form, specification §4 note 2); `Pi_Omega(k_max)` floor drops with `N`.

### WP-12 NetCDF I/O and restart
Status: done 865c126. ← WP-7, WP-8, WP-9. ADRs: 0009, 0024.
- `Output(; period, diagnostics)` per output class; writers for `global.nc`, `spectral.nc`, `fields/snap_NNNNNN.nc`; parameters as global attributes; `params.toml`; `read_snapshot!`.
- Accept: written values read back exactly (compare with `isequal`, because `inj_*` at `t₀` is NaN); the state rebuilt from a snapshot satisfies `‖û_restart − û‖₂/‖û‖₂ ≤ 10·eps(T)` for `T = Float64` and `Float32`; a run of `2n` steps and a run of `n` steps + restart + `n` steps, with OU forcing, use the identical forcing sequence and agree to `‖·‖₂`-relative `100·eps(T)` after `n = 10` steps (ADR-0009).

### WP-13 Simulation driver
Status: done ef3417b. ← WP-7, WP-12. ADRs: 0009, 0010, 0027.
- `Simulation`, `set_initial_condition!`, `run!`, `restart!`; parameter conversion to `T`; per-class sampling periods; `run.log` (TimerOutputs, `profile` flag); `max_walltime`; NaN and CFL checks every `n_check` steps; the NaN check inspects the state fields only, not diagnostics or stored rates.
- CFL number `C = Δt·max_x[(|u_x| + c)/Δx + (|u_y| + c)/Δy]`, with `c = 0` (incompressible) or `c = max_x c_s(ρ)` (compressible). The run stops when `C > C_stop`; default `C_stop = 1`, user-settable. Euler and RK2 are unstable on undamped oscillating modes at any `C`, so they are for tests, not production; the driver warns when a forced run uses them.
- Accept: a tiny forced run produces all files; a walltime of 0 writes one snapshot and stops; an injected NaN writes a snapshot and raises an error.

### WP-14 Balanced-density routine
Status: done ef3417b. ← WP-7, WP-11. ADRs: 0013, 0026.
- `balanced_density(u, eos, grid)`: solve `∇²h = −∂_i∂_j(u_iu_j)` for `h` up to a constant; find the constant `h₀` by a scalar root find so that `⟨density_from_enthalpy(h + h₀)⟩ = 1`; return `ρ = density_from_enthalpy(h + h₀)` truncated to the disc. (Shifting `ρ` after the inversion changes `∇h` for a nonlinear `h(ρ)`.)
- Accept, `μ = ζ = κ = 0`, cellular flow, `Ma = 0.1`, `Isothermal` and `Polytropic(γ = 1.4)`: `θ(0) = 0` to round-off; `∂_tθ(0)` relative to `max|k_ik_j(u_iu_j)^|` (the value for uniform `ρ`) is below `1e-9` at `N = 64`. Truncating `ρ̂` leaves a residual that falls spectrally with `N` (about `1e-8` at `N = 16`); below that sits a round-off floor that grows with `N` (about `2e-12` at `N = 64`, `1.5e-11` at `N = 256`). The wrong procedure (shifting `ρ`) gives about `1e-5`. With viscosity, `∂_tθ(0) = O(μ Ma²)`, because `∂_i(1/ρ)∇²u_i ≠ 0`.
- Accept: at `Ma = 0.1`, the maximum of `E_kin_comp` over `t ∈ [0, 10 max(Lx, Ly)/c_s]` is at least 10 times lower than with uniform `ρ`.

### WP-15 Visualisation extensions and acceptance animation
Status: todo. ← WP-13. ADRs: 0011, 0012.
- `plot_field`, `plot_spectrum`, `plot_global`, `animate_snapshots` in both extensions.
- Choose `ν` once for the crossed Kolmogorov run (ADR-0011) and fix it in `test/acceptance/crossed_kolmogorov.jl`. Produce the animation (vorticity plus spectrum tail) and add an entry `pending` to `test/acceptance/RESULTS.md`.
- Accept (user): the user records pass in `RESULTS.md`.

### WP-16 Benchmark suite
Status: todo. ← WP-13. ADRs: 0010.
- `benchmark/run.jl`: time per step and memory over `N` and `(model, scheme, device, T)`; fit `t_step = a·N² log N + b·N² + c` (the constant `c` absorbs kernel-launch latency on GPU); `benchmark/predict.jl` for a target run; TOML output.
- Accept (benchmark): fit residual under 10% on this machine for `N = 256 … 4096`, CPU and CUDA.

### WP-17 Documentation site
Status: todo. ← WP-13. ADRs: 0012.
- Documenter pages: getting started (the target script), equations (link to the specification), API from docstrings, output format, tests and GPU instructions, ADR index.
- Accept (docs): docs build in CI with no missing docstrings.

### WP-18 Intermittency analysis (late)
Status: todo. ← WP-13. ADRs: 0012.
- Offline module on FIELDS snapshots: increments, structure functions, increment PDFs. Specify in detail before starting.
