# Amendments from review round 2

Date: 2026-11-19. Status: accepted.

Review round 2 of `design/origin/ROADMAP-build.md` (report: `design/reviews/roadmap/r2.md`) found defects that remained after ADR-0014. This ADR records the resulting design changes; the details live in the roadmap work packages.

- **Hermitian `kx = 0` column (ADR-0007).** The real-FFT layout stores `(0, ky)` and `(0, −ky)`. A generator that fills them independently makes spectral diagnostics disagree with the physical field (4.4% injection error in the reviewer's probe). Every generated spectral field passes through `enforce_hermitian!`, and the tests measure injection on the physical grid.
- **Balanced density (ADR-0013).** The normalising constant goes on `h`, found by a scalar root find; `ρ` is truncated to the disc, so the `∂_tθ(0)` residual falls spectrally with `N` and the round-off test runs at `N = 64`.
- **Integer truncation test (ADR-0003).** The strict disc is evaluated on integer mode indices; a float comparison can keep the self-triad mode on rectangles.
- **Injection rate with the `½Δt|f̂|²` term.** For forcing held constant over a step, the plain product `Re[û*·f̂]` is biased by `O(Δt/τ_f)`.
- **Restart tolerance in a norm (ADR-0009).** `‖û_restart − û‖₂/‖û‖₂ ≤ 10·eps(T)`, valid for Float32 too.
- **Test environment.** `test/Project.toml` holds the test dependencies, so the per-file commands of `CLAUDE.md` work outside `Pkg.test`.
- **Five test tiers (ADR-0011).** Blocking, non-blocking, user, benchmark, docs.
- **Compressive forcing only in the compressible model (ADR-0007).**
