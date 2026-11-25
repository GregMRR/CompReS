# Amendments from review round 6

Date: 2026-11-19. Status: accepted.

Review round 6 of `design/origin/ROADMAP-build.md` (report: `design/reviews/roadmap/r6.md`) found that the forcing-only compressible test of ADR-0018 can miss a `ρ`-free formula on a single forced mode. Changes:

- **The test set-up is fixed:** a multi-mode band, a nonzero initial velocity, `N = 32`, 100 steps, both float types. The density `1 + 0.1 cos(x + y)` couples only modes that differ by `±(1, 1)`.
- **Tolerance scale `E_kin(t_{n+1})/Δt`.** At the first step from rest `E_kin(t_n) = 0`.
- **`Problem` is owned by WP-6.** `update_forcing!` takes the workspace.
- **Samples at `t₀`.** `run!` writes every output class at the start time; `inj_*` there is NaN as a plain value, not a NetCDF fill value. The NaN safeguard inspects only the state fields (ADR-0010).
