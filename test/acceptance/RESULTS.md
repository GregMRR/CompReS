# Acceptance results

Records the user-tier acceptance items. The user writes `pass` in the Result column after inspecting the item.

| Work package | Item | Result | Date |
|--------------|------|--------|------|
| WP-15 | Crossed Kolmogorov animation, `test/acceptance/crossed_kolmogorov.mp4` (vorticity and `E_kin_inc(k)`, 121 frames, t = 0 to 30) | pass | 2026-11-25 |

### WP-15: crossed Kolmogorov animation

- Run: incompressible, unforced, `u = sin(2y) x̂ + sin(3x) ŷ` (`A₁ = A₂ = 1`, `k₁ = 2`, `k₂ = 3`), `Lx = Ly = 2π`, `N = 512` (`k_max = 170`), Float64, RK4 with the integrating factor, `Δt = 10⁻³`, `t_end = 30`, `ν = 4·10⁻⁴`.
- Measured: palinstrophy peak `P = 316` at `t = 3.2`; `k_max/k_ν = 4.3`; `Re = A₁/(k₁ν) = 1250`; `E_kin_inc(k = 160)` at most `4·10⁻¹¹` of the spectrum peak; `E_kin` falls from 0.500 to 0.446 and `Z` from 3.25 to 1.92 by `t = 30`. ADR-0029 gives the choice of `ν`.
- Look for: the shear layers roll up into vortices; the vorticity stays free of grid-scale stripes; the spectrum tail falls steeply towards `k_max` (dashed line) at all times.
- Reproduce: `julia --project=test test/acceptance/crossed_kolmogorov.jl` (about 2 min on the RTX 5070 with CUDA). The output files go to `runs/acceptance_crossed_kolmogorov` (490 MB, ignored by git).
