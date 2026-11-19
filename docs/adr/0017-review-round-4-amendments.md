# Amendments from review round 4

Date: 2026-11-19. Status: accepted.

Review round 4 of `docs/ROADMAP.md` (report: `reviews/docs/ROADMAP/r4.md`) found that no work package owned the start-of-step injection rates of ADR-0016. Changes:

- **`record_injection!` hook in `step!` (WP-7).** `step!` calls it after `update_forcing!`, before the first stage. Its default does nothing, so `step!` stays generic. WP-8 adds the incompressible method, WP-11 the compressible one.
- **Injection test.** With forcing only, every RK scheme gives `u_{n+1} = u_n + Δt f_n`, so `inj_E_kin` must equal `ΔE_kin/Δt` to round-off.
- **OU test interval.** Realisation means are independent, time samples are not; the test uses a 99.9% interval over realisation means with a fixed seed.
