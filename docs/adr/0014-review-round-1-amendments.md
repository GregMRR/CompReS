# Amendments from review round 1

Date: 2026-11-19. Status: accepted.

An adversarial review of `docs/ROADMAP.md` (report: `reviews/docs/ROADMAP/r1.md`) found defects in the design. This ADR records the resulting changes; the amended ADRs point here.

- **Restart tolerance (ADR-0009).** A snapshot stores physical-space fields, and the physical-to-spectral round trip changes the last bit of most modes. A restart therefore reproduces the saved state to relative `1e-14`, not bit for bit. The forcing sequence stays exact. The user chose this over storing the spectral state as well, which would double the snapshot size.
- **Strict truncation (ADR-0003).** The disc `|k| ≤ N/3` aliases the self-triad at `|k| = N/3` when 3 divides `N`. The truncation is now `|k| < k_c`, `k_c = N/3` (units `2π/L`); `k_max` is the largest retained `|k|`.
- **Forcing constructors take the grid (ADR-0007).** The constructor must count modes and check them against the truncation, so it needs the grid. `band` uses physical wavenumbers.
- **Five EOS methods (ADR-0008).** `internal_energy` is needed for `E_int`.
- **Makie extension trigger (ADR-0012).** An extension loads only when all its triggers are loaded; `NS2DMakieExt` triggers on Makie.
- **IF order test (ADR-0005, ADR-0011).** The decay and conservation tests leave the integrating-factor stage algebra untested; every scheme now gets an order test with the IF on.
- **CFL definition (ADR-0005, ADR-0010).** Defined in ROADMAP WP-13, including the sound speed; Euler and RK2 serve tests only.
- **Balanced density (ADR-0013).** `∂_tθ(0) = 0` holds for the inviscid case; viscosity leaves `O(μ Ma²)`.
