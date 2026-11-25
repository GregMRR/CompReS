# Balanced density: interface, truncation of the enthalpy and root find

Date: 2026-11-19. Status: accepted.

WP-14 fixes the parts of the balanced-density routine that ADR-0013 leaves open. Changes:

- **Interface.** `balanced_density(u, eos, grid)` takes the `NamedTuple` `(; ux, uy)` of the initial-condition functions and returns `(; rho)`, as `uniform_density` does. The user merges the two: `merge(u, balanced_density(u, eos, grid))` gives the fields of a compressible `State` and the input of `set_initial_condition!`. The function converts `eos` to the float type of the grid (ADR-0006).
- **The enthalpy `ĥ` is truncated to the disc before the inversion.** The grid products `u_iu_j` of a velocity inside the truncation have aliases only outside the disc (the 2/3 rule), so `ĥ` is exact on the disc and the modes outside hold aliases. `ρ = ρ(h + h₀)` is then computed from the band-limited `h`, and `ρ̂` is truncated again. The measured `∂_tθ(0)` residual falls spectrally with `N` and agrees with the values of ROADMAP WP-14.
- **Root find.** Newton steps on `f(h₀) = ⟨ρ(h + h₀)⟩ − 1`, with `f' = ⟨ρ/c_s²(ρ)⟩`, start at `h₀ = h(1)`. A step outside the bracket `[h(1) − max h, h(1) − min h]`, where `f` changes sign, becomes a bisection step. The iteration stops when the step is below `eps(T)·c_s²(1)`, which changes `⟨ρ⟩` by less than `eps(T)`. The grid means are device reductions; the root find costs a few reductions once per run.
- **No new workspace field.** The routine runs once, before the time loop, so it allocates its own `Transforms` and work arrays instead of taking a `Workspace`.
