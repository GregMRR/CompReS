# Amendments from review round 7 (confirmation round)

Date: 2026-11-19. Status: accepted.

The confirmation round on `design/origin/ROADMAP-build.md` (report: `design/reviews/roadmap/r7.md`) found one more defect. Changes:

- **Compressible global quantities are grid means (ADR-0009).** `w = √ρ u` and `s = ω/√ρ` are not band-limited. Truncating `ŵ` or `ŝ` drops weight (34 eps for `E_kin`, 1300 eps for `Omega` at `N = 32`), so the injection identity fails. Compressible spectra and fluxes therefore sum over all grid modes, and GLOBAL values are grid means; by Parseval they agree.
- **Forcing keyword `compressive_fraction`;** the forcing-only test states all its parameters.
- **The NetCDF read-back test uses `isequal`,** because `inj_*` at `t₀` is NaN.
