# Amendments from review round 5

Date: 2026-11-19. Status: accepted.

Review round 5 of `docs/ROADMAP.md` (report: `reviews/docs/ROADMAP/r5.md`) found that the compressible injection rates were tested only at `ρ = 1`. Changes:

- **`nonlinear = false` test switch on the model types.** It keeps only the linear terms and the forcing; in the compressible model `ρ` stays fixed. Forcing-only tests then check `inj = ΔE/Δt` exactly, also with nonuniform `ρ`.
- **Injection tolerance scaled by `E_kin/Δt`.** The energy difference cancels, so a tolerance relative to `ΔE/Δt` fails correct code.
- **`problem` and the stored `f̂`.** `step!(state, problem, workspace, Δt)`; `update_forcing!` stores `f̂` in the workspace for `record_injection!`. Stored rates start as NaN, which GLOBAL writes as the fill value before the first step of a run segment.
- **Compressible `rhs!` calls `add_forcing!`.**
