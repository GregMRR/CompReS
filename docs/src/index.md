# CompReS

CompReS is a Julia pseudospectral solver for forced-dissipated 2D turbulence.
It covers the incompressible and the compressible (low Mach) model.
Its purpose is to study intermittency as a function of the Mach number.
One source runs on CPU threads, NVIDIA GPUs and AMD GPUs.

A run is a short Julia script against the package API.
There is no input-file format.

## Pages

- [Getting started](getting_started.md): one complete run script, with each argument explained, and the restart line.
- [Equations and conventions](equations.md): the two models, the FFT normalisation, the truncation and the shells.
- [Output format](output.md): the files of a run and how to read them from Julia and Python.
- [Tests and GPU](tests_and_gpu.md): test tiers, commands, and how to run on NVIDIA and AMD GPUs.
- [API reference](api.md): every documented symbol of the package.
- [Architecture decision records](adr.md): the design decisions, in chronological order.

## Sources

- The binding equations are in the specification: [`ns2d_equations.pdf`](ns2d_equations.pdf) (source: [`ns2d_equations.tex`](https://github.com/GregMRR/CompReS/blob/main/docs/ns2d_equations.tex)).
- The project vocabulary is in [`GLOSSARY.md`](https://github.com/GregMRR/CompReS/blob/main/GLOSSARY.md).
- The work plan (open items, ideas and implemented items) is in [`PLAN.md`](https://github.com/GregMRR/CompReS/blob/main/PLAN.md).

```@docs
CompReS
```
