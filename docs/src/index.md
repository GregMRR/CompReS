# CompReS

CompReS is a Julia pseudospectral solver for forced-dissipated 2D turbulence.
It covers the incompressible and the compressible (low Mach) model.
One source runs on CPU threads, NVIDIA GPUs and AMD GPUs.

- The binding equations are in `docs/ns2d_equations.tex` (specification).
- The design decisions are in `docs/adr/` (architecture decision records).
- The work packages are in `docs/ROADMAP.md`.

```@docs
CompReS
```
