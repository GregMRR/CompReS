# Amendments from review round 3

Date: 2026-11-19. Status: accepted.

Review round 3 of `design/origin/ROADMAP-build.md` (report: `design/reviews/roadmap/r3.md`) found two acceptance tests that a correct implementation fails, and smaller gaps. Changes:

- **Injection rates are computed in `step!` (ADR-0009).** A diagnostic sees `u_{n+1}`, but the injection rate needs `u_n` with `f_n`; sampling after the step biases it by `+Δt/τ_f`. `step!` stores the rates in the workspace.
- **Balanced-density bound.** The round-off floor of `∂_tθ(0)` grows with `N` (`h(1)` dominates `h`), so the WP-14 bound is `1e-9` at `N = 64`, well above the floor and well below the wrong-procedure value.
- **Shell assignment on integer indices (ADR-0003).** Float rounding misassigns modes that lie exactly on a shell edge.
- **`enforce_hermitian!`** reads the `ky > 0` half and also zeroes the Nyquist column and row.
- **OU normalisation test runs on CPU only (ADR-0011)**, sampled once per `τ_f`; CUDA checks the same forcing sequence.
- **`docs/Project.toml`** loads CompReS through `[sources]`. WP-5 depends on WP-3.
