# Equations and conventions

The binding equations are in the specification:
[`ns2d_equations.pdf`](ns2d_equations.pdf) (LaTeX source: [`ns2d_equations.tex`](https://github.com/GregMRR/CompReS/blob/main/docs/ns2d_equations.tex)).
This page summarises the models and the conventions.
It does not repeat the derivations; the docstrings cite the equation names of the specification.
ADR-0002 records the choice of the specification as the binding source.

## Models

Both models live on a periodic box with space average ``\langle f\rangle = L^{-2}\int f\,d^2x``.
All global quantities are per unit area.
The forcing ``\mathbf f`` is an acceleration.

**Incompressible** (specification, Sec. 2).
The density is the constant ``\rho_0``, and ``\nu = \mu/\rho_0``.
The velocity is solenoidal.
The invariants of the inviscid, unforced, truncated system are the energy ``E`` and the enstrophy ``Z = \tfrac12\langle\omega^2\rangle``.
Parameters: `Incompressible(; ν, large_scale)`.

**Compressible, barotropic** (specification, Sec. 3).
The density ``\rho`` obeys a continuity equation with the density diffusivity ``\kappa``.
The pressure is a function of the density, ``p = p(\rho)``, given by the equation of state.
The viscous stress uses the shear viscosity ``\mu`` and the bulk viscosity ``\zeta``.
`Isothermal(cs)` and `Polytropic(…)` are the built-in equations of state; the enthalpy ``h``, the internal energy ``e`` and the sound speed ``c_s`` follow from ``p(\rho)``.
The vorticity-like invariant is ``\Omega = \tfrac12\langle\omega^2/\rho\rangle``, not ``Z``.
The Mach number is ``\mathrm{Ma} = u_\mathrm{rms}/c_s(\rho_0)``.
Parameters: `Compressible(; eos, μ, ζ, κ, large_scale)`.

**Large-scale dissipation.** `LargeScaleDissipation(α, m)` adds the rate ``\alpha k^{-2m}`` on every mode except ``\mathbf k = 0``, in both models (ADR-0004).

## Fourier convention

- The Fourier series is ``f(\mathbf x) = \sum_{\mathbf k}\hat f(\mathbf k)e^{i\mathbf k\cdot\mathbf x}`` with ``\hat f(\mathbf k) = \langle f e^{-i\mathbf k\cdot\mathbf x}\rangle``.
  `to_spectral!` returns these coefficients, so the ``\mathbf k = 0`` mode is the space average.
- Parseval reads ``\langle fg\rangle = \sum_{\mathbf k}\hat f^*(\mathbf k)\hat g(\mathbf k)`` for real ``f`` and ``g``.
- The code uses the real FFT. It stores half of the ``k_x`` axis, so an array has size `(Nx÷2+1, Ny)`.
  Parseval sums double every mode except `kx = 0` and the Nyquist column.
- The column `kx = 0` stores both ``(0, k_y)`` and ``(0, -k_y)``.
  A spectral field generated directly must pass through `enforce_hermitian!`.

## Truncation

The Galerkin truncation keeps the modes with ``|\mathbf k| < k_c``, ``k_c = (N_x/3)(2\pi/L_x)``, with a strict inequality (ADR-0003).
`k_max` is the largest retained ``|\mathbf k|``.
The test for a mode uses integer mode indices, never float wavenumbers.
The grid must have equal spacing in `x` and `y`, so `Lx/Nx = Ly/Ny`; a rectangular box is allowed.

## Shells

The shell width is ``\Delta k = 2\pi/L``.
A mode ``\mathbf q`` is in shell ``k_n = n\Delta k`` if ``n - \tfrac12 \le |\mathbf q|/\Delta k < n + \tfrac12``.
A mode on an edge goes to the upper shell.
The shell spectrum of a field ``a`` is ``A(k_n) = \tfrac12\sum_{\mathbf q\in k_n}|\hat a(\mathbf q)|^2``.
The shells extend to the grid corner.
The grid holds the shell centres `shell_k`, the upper edges `shell_edge` and the mode counts `shell_count` (retained modes) and `shell_count_all` (all modes).

## Fluxes

If ``\partial_t\hat a = -\hat N + \dots``, the flux through ``K`` is ``\Pi_a(K) = \sum_{|\mathbf q|\le K}\mathrm{Re}[\hat a^*(\mathbf q)\hat N(\mathbf q)]``.
A positive ``\Pi_a`` means transfer to ``|\mathbf q| > K``.

- The incompressible fluxes vanish at `k_max`.
- The compressible fluxes carry an aliasing floor at `k_max` (specification, Sec. 4).
- The compressible spectra sum over all grid modes, because ``\sqrt\rho\,\mathbf u`` and ``\omega/\sqrt\rho`` are not band-limited (ADR-0020).

## Time stepping

`step!` advances the state with an explicit Runge-Kutta scheme.
With `integrating_factor = true`, the propagator ``e^{L\tau}`` integrates the linear terms exactly and the scheme integrates the nonlinear term and the forcing (ADR-0005, ADR-0023).
