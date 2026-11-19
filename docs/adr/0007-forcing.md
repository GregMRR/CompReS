# Forcing

Date: 2026-11-19. Status: accepted.

A forcing is an acceleration `f` that acts on an explicit **forcing support**: a set of wavevectors built by `band(k_start, k_end)` (`k_start ≤ |k| ≤ k_end`, physical wavenumbers), `modes([(kx, ky), …])` or `support(pred)`. Every forcing constructor takes the grid as its first argument (amended by ADR-0014). The constructor reports the mode count and rejects an empty support or one outside the Galerkin truncation. The forcing stores an index list, so its cost scales with the number of forced modes, not with `N²`.

Types:

- `KolmogorovForcing`: `F sin(k_f y) x̂` with `k_f = n·2π/Ly`; support `(0, ±k_f)`.
- `OUForcing`: one independent complex Ornstein–Uhlenbeck process per mode, parameters `(ε, τ_f)`. The variance makes the mean injection rate equal `ε` for short `τ_f`. Exact OU update once per step.
- `RandomPhaseForcing`: fixed amplitudes with rms `f_rms`; all phases are redrawn together every `τ_r` and are constant in between.

The amplitude is uniform over the support, with an optional profile `amp(k)`. In the compressible model the forcing is solenoidal by default; a parameter sets the compressive fraction.

A new forcing type is one struct plus its methods (same pattern as the EOS, ADR-0008).

Random numbers come from a generator seeded by `(seed, step)` at each step. A restart then needs only the seed, the step and the forcing state (OU values, current phases and the time of the last reset), all of which the snapshot stores (ADR-0009).
