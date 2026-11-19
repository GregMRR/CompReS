# CompReS

Pseudospectral simulation of forced-dissipated 2D turbulence, compressible (low Mach) and incompressible, in a periodic rectangle. The reference equations are in `docs/ns2d_equations.tex`.

## Physical model

**Model**:
The set of evolved equations: either *incompressible* (velocity only, divergence-free) or *compressible* (velocity and density, barotropic).
_Avoid_: mode, case, regime

**Equation of state (EOS)**:
The user-chosen barotropic law `p(ρ)`, for example isothermal or polytropic.
_Avoid_: pressure law, closure

**Mach number**:
`Ma = u_rms / c_s(ρ₀)`, with `c_s² = p'(ρ)`.

**Small-scale dissipation**:
The Navier–Stokes viscous term: `ν∇²u` (incompressible) or `(1/ρ)∂_jτ_ij` with shear viscosity `μ` and bulk viscosity `ζ` (compressible).
_Avoid_: viscous term, viscosity (alone), hyperviscosity

**Dissipation rate**:
The rate at which a dissipation term removes a quadratic quantity, for example `2νZ` for `E_kin` in the incompressible model, or `⟨D⟩` with `D = τ_ij∂_ju_i` in the compressible model. Each rate is labelled small-scale or large-scale after the term that produces it.
_Avoid_: dissipation (alone, when the rate is meant)

**Injection rate**:
The rate at which the forcing adds a quadratic quantity, for example `inj_E_kin`, measured at the start of a step with the forcing of that step.
_Avoid_: input, production

**Large-scale dissipation**:
The user-set damping `−α(−∇²)^{−m}u` that removes energy at large scales; `m = 0` is linear friction, `m ≥ 1` is hypoviscosity.
_Avoid_: drag, friction (alone)

**Density diffusivity**:
The coefficient `κ` of the regularising term `κ∇²ρ` in the continuity equation.

## Forcing

**Forcing**:
An acceleration `f` (force per unit mass) added to the velocity equation.
_Avoid_: source, stirring

**Kolmogorov forcing**:
`f = F sin(k_f y) x̂`, a single-mode shear forcing.

**OU forcing**:
Forcing in a wavenumber band where each mode follows an independent Ornstein–Uhlenbeck process with correlation time `τ_f`.

**Random-phase forcing**:
Forcing in a wavenumber band with fixed amplitudes and phases redrawn every `τ_r`, constant in between.

**Forcing support**:
The set of wavevectors on which a forcing acts. Kolmogorov forcing has support `(0, ±k_f)`.
_Avoid_: forcing range, forced shell

**Forcing band**:
A forcing support of the form `k_start ≤ |k| ≤ k_end`, in physical wavenumbers.

## Spectral quantities

**Galerkin truncation**:
The set of retained Fourier modes, the disc `|k| < k_c` with `k_c = (Nx/3)(2π/Lx)`; `k_max` is the largest retained `|k|`.
_Avoid_: dealiasing mask, cutoff

**Shell**:
The set of modes `k` with `n − ½ ≤ |k|/Δk < n + ½`, labelled by `k_n = nΔk`.
_Avoid_: bin, ring

**Spectrum**:
The shell sum `½ Σ_{q ∈ k_n} |â(q)|²` of a field `a`.

**Flux**:
`Π_a(K) = Σ_{|q| ≤ K} Re[â*(q) N̂(q)]`, the nonlinear transfer of the quadratic quantity of `a` from `|q| ≤ K` to `|q| > K`. The kinetic-energy flux splits exactly as `Π_kin = Π_kin^inc + Π_kin^comp`; only the sum vanishes at the Galerkin truncation.
_Avoid_: transfer (for the cumulative quantity)

**Solenoidal part / compressive part**:
The Helmholtz components of a vector field: `P_ij(k)â_j` (divergence-free) and `k_ik_jâ_j/k²` (curl-free). A quantity built on the solenoidal part has the superscript `inc`; one built on the compressive part has the superscript `comp`.
_Avoid_: rotational/irrotational, vortical/acoustic, longitudinal/transverse

**Kinetic energy**:
`E_kin = ½⟨ρ|u|²⟩ = ½⟨|w|²⟩` with `w = √ρ u`, split as `E_kin = E_kin^inc + E_kin^comp`. In the incompressible model `E_kin = E_kin^inc`.
_Avoid_: E_K, E (alone)

**Internal energy**:
`E_int = ⟨ρ e(ρ)⟩`, compressible model only.

**Total energy**:
`E_tot = E_kin + E_int`.

**Enstrophy**:
`Z = ½⟨ω²⟩`. Conserved only in the incompressible model.

**Palinstrophy**:
`P = ½⟨|∇ω|²⟩`. In the incompressible model the small-scale dissipation removes `Z` at the rate `2νP`.

**Potential enstrophy**:
`Ω = ½⟨ω²/ρ⟩`, the quadratic invariant of the compressible model.
_Avoid_: enstrophy (in the compressible model)

## Time stepping

**Scheme**:
An explicit Runge–Kutta method given by its Butcher table: `Euler`, `RK2`, `SSPRK3` or `RK4`.
_Avoid_: integrator, method (alone)

**Stage**:
One evaluation of the right-hand side inside a step, at the stage value `û^(i)` and node `c_i` of the Butcher table.

**Integrating factor (IF)**:
The exact integration of the linear terms by the propagator `e^{Lτ}` inside each step (Lawson form); the scheme then integrates only the nonlinear term and the forcing.
_Avoid_: exponential integrator, ETD

**Time stepper**:
The `TimeStepper` in the workspace: the scheme, the IF switch and the stage buffers.

## Runs and output

**Run**:
One simulation from an initial condition or restart to its end time, with one set of parameters.

**Output class**:
A named group of diagnostics sampled with one common period: GLOBAL, SPECTRAL or FIELDS.

**Diagnostic**:
A named function `(state, problem, workspace) → value` registered under the output class GLOBAL (one number) or SPECTRAL (one value per shell).

**Snapshot**:
One FIELDS sample: all independent physical-space fields at one time, plus the forcing state. A snapshot is also a valid restart input.
_Avoid_: dump, checkpoint

**Runtime log**:
The per-run record of wall-clock timings of the solver tasks.
