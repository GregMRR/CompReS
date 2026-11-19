# Getting started

## Install

CompReS is not in the General registry.
Clone the repository, then start Julia in its directory:

```
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Use `julia -t auto` to run the CPU backend on all threads.

## A compressible run

This script runs the compressible model with Ornstein-Uhlenbeck forcing.
The values `Nx = Ny = 512` and `t_end = 10.0` give a production run; for a first test use `Nx = Ny = 32` and `t_end = 0.2`.

```julia
using CompReS

grid    = Grid(; Lx = 2π, Ly = 2π, Nx = 512, Ny = 512, T = Float64, backend = CPU())
model   = Compressible(; eos = Isothermal(cs = 10.0), μ = 1e-4, ζ = 0.0, κ = 0.0,
                         large_scale = LargeScaleDissipation(α = 0.1, m = 0))
forcing = OUForcing(grid, band(4, 6); ε = 1.0, τ_f = 0.05, seed = 1)

sim = Simulation(grid, model; forcing, scheme = RK4(), integrating_factor = true, Δt = 1e-4,
                 outputs = (GLOBAL   = Output(period = 0.01),
                            SPECTRAL = Output(period = 0.1),
                            FIELDS   = Output(period = 1.0)),
                 outdir = "runs/example", max_walltime = 23.5 * 3600)

set_initial_condition!(sim, random_velocity(grid; spectrum = k -> k^4 * exp(-k^2 / 8),
                                            E_kin = 0.5, seed = 2))
run!(sim; t_end = 10.0)
```

To restart, build the same `Simulation` again and read a snapshot:

```julia
sim = Simulation(grid, model; forcing, scheme = RK4(), integrating_factor = true, Δt = 1e-4,
                 outputs = (GLOBAL   = Output(period = 0.01),
                            SPECTRAL = Output(period = 0.1),
                            FIELDS   = Output(period = 1.0)),
                 outdir = "runs/example", max_walltime = 23.5 * 3600)
restart!(sim, "runs/example/fields/snap_000010.nc")
run!(sim; t_end = 20.0)
```

The script was run at `N = 32` with `t_end = 0.2`, followed by a restart from `snap_000000.nc` and `run!(sim; t_end = 0.3)`.
It needs no change from the target script.

## The arguments

### Grid

`Grid(; Lx, Ly, Nx, Ny, T, backend)` fixes the periodic box, the number of grid points and the float type.

- `Lx`, `Ly`: box lengths. The grid spacings `Lx/Nx` and `Ly/Ny` must be equal.
- `Nx`, `Ny`: grid points per direction. Both must be even.
- `T`: float type of all fields. A `Float32` run never promotes to `Float64` (ADR-0006).
- `backend`: `CPU()` here. For GPUs see [Tests and GPU](tests_and_gpu.md).

The Galerkin truncation keeps the modes with `|k| < Nx/3 · 2π/Lx`.

### Model

`Compressible(; eos, μ, ζ, κ, large_scale)` is the compressible model; `Incompressible(; ν, large_scale)` replaces it for the incompressible model.

- `eos`: the equation of state, `Isothermal(cs = …)` or `Polytropic(; K, γ)`.
- `μ`, `ζ`: shear and bulk viscosity.
- `κ`: density diffusivity.
- `large_scale`: `LargeScaleDissipation(α, m)` adds the damping rate `α k^(−2m)` to every mode except `k = 0`. `m = 0` is linear friction; `m ≥ 1` is hypoviscosity.

Write the parameters as plain Julia numbers.
`Simulation` converts them to the float type `T` of the grid.

### Forcing

`OUForcing(grid, band(4, 6); ε, τ_f, seed)` is an Ornstein-Uhlenbeck forcing.

- `band(4, 6)`: the forced modes, in physical wavenumber `|k|`.
- `ε`: the mean injection rate of energy.
- `τ_f`: the correlation time of the forcing.
- `seed`: the seed of the forcing noise. A restart reproduces the forcing sequence.

Other forcings are `KolmogorovForcing` and `RandomPhaseForcing`; `NoForcing()` is the default.

### Simulation

- `scheme`: the Runge-Kutta scheme, `Euler()`, `RK2()`, `SSPRK3()` or `RK4()`.
- `integrating_factor`: when `true`, the step integrates the linear terms exactly (Lawson form).
- `Δt`: the time step. It has no default.
- `outputs`: one `Output(period = …)` per output class. A class that is absent writes no file. See [Output format](output.md).
- `outdir`: the run directory. It has no default. A new run throws an `ArgumentError` if `outdir` already holds output files; pass `overwrite = true` to delete them.
- `max_walltime`: wall-clock limit in seconds. At the limit, `run!` writes a snapshot and returns.

Further keywords are `n_check`, `C_stop` and `profile` (see [`Simulation`](@ref)).
Every `n_check = 10` steps the driver checks that the state is finite, that the density is positive and that the CFL number is at most `C_stop = 1`.

### Initial condition and run

- `random_velocity(grid; spectrum, E_kin, seed)`: a random solenoidal velocity field with the given spectrum and kinetic energy. Other initial conditions are `kolmogorov`, `crossed_kolmogorov` and `uniform_density`. The initial density is uniform, `ρ = 1`, unless the initial condition sets `rho`.
- `run!(sim; t_end)` integrates to the step nearest to `t_end`.
- `restart!(sim, path)` loads the state, the time and the forcing state from a snapshot. The next `run!` continues the run in the same directory.

The run directory holds the files described in [Output format](output.md), and the runtime log `run.log`.
