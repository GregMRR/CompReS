# Balanced-density routine

Date: 2026-11-19. Status: accepted.

A uniform `ρ` with a solenoidal `u` emits acoustic waves of amplitude `O(Ma²)`. The routine `balanced_density(u, eos, grid)` returns the `ρ` that makes `θ = 0` and `∂_tθ = 0` at `t = 0`: it solves `∇²h(ρ) = −∂_i∂_j(u_iu_j)` spectrally, applies `density_from_enthalpy`, and shifts the constant so that `⟨ρ⟩ = 1`. Higher-order balancing (a compressive initial `u^c`) is out of scope.

The routine is standalone. No run and no initial condition calls it by default; the user calls it for chosen initial conditions.
