# Dissipation terms

Date: 2026-11-19. Status: accepted.

**Small-scale dissipation** is the Navier–Stokes term only (no hyperviscosity option). Incompressible: `ν∇²u`. Compressible: `(1/ρ)(μ∇²u + ζ∇θ)` with dynamic shear and bulk viscosities, because only this form gives a sign-definite dissipation rate `D ≥ 0` (specification, note 5).

The compressible term splits into a linear part with `ρ₀ = 1`, `μ∇²u + ζ∇θ`, and an explicit remainder `(1/ρ − 1)(μ∇²u + ζ∇θ)` of order `Ma²`. The integrating factor (ADR-0005) treats the linear part exactly. Per mode, the linear part is diagonal in the Helmholtz basis: rate `μk²` on the solenoidal part and `(μ + ζ)k²` on the compressive part.

**Large-scale dissipation** is `−α(−∇²)^{−m}u` on all modes except `k = 0`: `m = 0` is linear friction, `m ≥ 1` is hypoviscosity. It acts on `u` in both models.

**Density diffusivity** `κ∇²ρ` exists with default `κ = 0`. Its terms `R^w`, `R^s` enter the budget diagnostics.

All linear rates are one precomputed spectral multiplier array per field, so a new linear term costs one line.
