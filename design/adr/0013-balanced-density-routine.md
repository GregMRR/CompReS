# Balanced-density routine

Date: 2026-11-19. Status: accepted.

A uniform `ρ` with a solenoidal `u` emits acoustic waves of amplitude `O(Ma²)`. The routine `balanced_density(u, eos, grid)` returns the `ρ` that makes `θ = 0` and `∂_tθ = 0` at `t = 0` for `μ = ζ = κ = 0`, up to a truncation residual that falls spectrally with `N` (with viscosity, `∂_tθ(0) = O(μ Ma²)`): it solves `∇²h = −∂_i∂_j(u_iu_j)` spectrally, adds the constant `h₀` that gives `⟨density_from_enthalpy(h + h₀)⟩ = 1` (scalar root find), and returns `ρ = density_from_enthalpy(h + h₀)` truncated to the disc. The constant goes on `h`, not on `ρ`: shifting `ρ` after the inversion changes `∇h` for a nonlinear `h(ρ)` (amended by ADR-0015). Higher-order balancing (a compressive initial `u^c`) is out of scope.

The routine is standalone. No run and no initial condition calls it by default; the user calls it for chosen initial conditions.
