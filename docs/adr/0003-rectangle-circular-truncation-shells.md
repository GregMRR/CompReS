# Rectangular box, circular Galerkin truncation, shell width

Date: 2026-11-19. Status: accepted.

The domain is a periodic rectangle `Lx × Ly` on an `Nx × Ny` grid with equal spacing, `Lx/Nx = Ly/Ny`. The grid constructor rejects other inputs.

The Galerkin truncation is a circular disc `|k| ≤ k_max`, with `k_max = (Nx/3)(2π/Lx)`. The disc removes all quadratic aliasing, it is isotropic, and the energy and enstrophy fluxes vanish exactly at `k_max`. A square 2/3 rule keeps about 21% more modes, but its fluxes vanish only at `√2·k_max`, and it is anisotropic.

Shells have width `Δk = 2π/max(Lx, Ly)` and use the rounding rule of the specification. On a rectangle, some low shells hold modes of one axis only, so SPECTRAL output records the mode count per shell.

## Considered options

- Square 2/3 rule: rejected; see above.
- Unequal grid spacing: rejected; the disc would waste the modes outside it along the finer axis.
