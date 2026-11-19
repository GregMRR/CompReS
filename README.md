# CompReS

CompReS is a Julia pseudospectral solver for forced-dissipated 2D turbulence, compressible (low Mach) and incompressible.
It studies intermittency as a function of the Mach number.
One source runs on CPU threads, NVIDIA GPUs and AMD GPUs.

Documentation: <https://gregmrr.github.io/CompReS/dev/> (getting started, equations, output format, API, design decisions).

## Quick start

```julia
using CompReS
grid  = Grid(; Lx = 2π, Ly = 2π, Nx = 64, Ny = 64, T = Float64, backend = CPU())
model = Incompressible(; ν = 1e-3, large_scale = LargeScaleDissipation(α = 0.1, m = 0))
sim   = Simulation(grid, model; Δt = 1e-3, outputs = (GLOBAL = Output(period = 0.1),), outdir = "runs/quick")
set_initial_condition!(sim, random_velocity(grid; spectrum = k -> k^4 * exp(-k^2 / 8), E_kin = 0.5, seed = 2))
run!(sim; t_end = 1.0)
```

Run the tests with `julia --project -e 'using Pkg; Pkg.test()'`.
Build the docs with `julia --project=docs docs/make.jl`.
