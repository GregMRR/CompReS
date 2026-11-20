# Test policy

Date: 2026-11-19. Status: accepted.

Every test is a function of the backend. `test/runtests.jl` runs each test on the CPU and, when `CUDA.functional()`, on CUDA. `test/amdgpu/runtests.jl` runs the same functions on `ROCBackend()`; it is stored and not run until AMD hardware is available.

Five tiers (a work package names the tier of each acceptance item; amended by ADR-0015):

1. **Blocking** (`Pkg.test`, must pass on CPU and CUDA). Exception: a statistical test that is too slow on GPU runs on CPU only, and its CUDA counterpart checks that the device reproduces the CPU sequence for the same seed (amended by ADR-0016, ADR-0017):
   - Integrating factor: every scheme reaches its order on `y' = −2y + y²` with the IF on, through the production `step!`.
   - Incompressible fluxes: `Pi_kin(k_max) = Pi_Z(k_max) = 0` to round-off on a random field projected to `θ = 0`.
   - Conservation: unforced, `ν = 0`, incompressible. Measure `|E_kin(T) − E_kin(0)|/E_kin(0)` and the same for `Z` at fixed `T` for four values of `Δt`. Fit the slope of the log–log plot. Pass when the slope is at least `r − 0.2` for Euler (`r = 1`), RK2 (`r = 2`) and RK4 (`r = 4`). A lower bound accepts superconvergence (RK schemes can converge one order faster on quadratic invariants). Choose `Δt` so that every error stays above `1e-12`.
   - Kolmogorov decay: a single-mode Kolmogorov initial condition decays as `exp(−νk²t)`; with the IF it matches to round-off.
   - All module unit tests.
2. **Non-blocking** (`test/nonblocking/runtests.jl`, reports only): compressible `Pi_kin` summed over all grid modes is 0 to round-off (Lamb form, specification §4 note 2); compressible `Pi_Omega` at `k_max` drops as `N` grows (aliasing floor).
3. **User acceptance** (`test/acceptance/`): an animation of the crossed Kolmogorov initial condition `A₁ sin(k₁y) x̂ + A₂ sin(k₂x) ŷ`, `k₁ ≠ k₂`, unforced, small `ν`. The user inspects it and records pass or reject in `test/acceptance/RESULTS.md`. The value of `ν` is chosen once, offline, as the smallest value with `k_max/k_ν ≥ 1.5` at the palinstrophy peak (`k_ν = (2νP/ν³)^{1/6}`) and `Re = U/(k₁ν) ≥ 100`, and is then fixed in the script (amended by ADR-0029: the resolution test is the grid-scale vorticity, below 1% of `max|ω|`).

4. **Benchmark** (`benchmark/`): a stated bound, for example the fit residual of WP-16.
5. **Docs**: `julia --project=docs docs/make.jl` builds without error.

Each module has unit tests on small, hand-checkable examples (single Fourier modes, `N = 8` or `16` grids).
