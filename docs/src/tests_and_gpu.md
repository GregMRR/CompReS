# Tests and GPU

## Test tiers

The test policy is in ADR-0011.
Every test is a function `run_tests(backend)`, so the CPU, CUDA and AMD runs share it.

1. **Blocking.** Must pass on the CPU and, when `CUDA.functional()`, on CUDA.
   Contents: unit tests per module on small hand-checkable examples (`N = 8` or `16`), and the integration tests of ADR-0011.
2. **Non-blocking.** Reports results and always exits 0.
   Contents: the compressible fluxes (sum over all grid modes, aliasing floor at `k_max`).
3. **User acceptance** (`test/acceptance/`). An animation that the user inspects; the result goes to `test/acceptance/RESULTS.md`.
4. **AMD.** The blocking tests on `ROCBackend()`. Stored; it needs AMD hardware.
5. **Docs.** `julia --project=docs docs/make.jl` builds without error.

## Commands

| Tier | Command |
|---|---|
| All blocking tests | `julia --project -e 'using Pkg; Pkg.test()'` |
| One test file | `julia --project=test test/runtests.jl <name>` |
| Non-blocking | `julia --project=test test/nonblocking/runtests.jl` |
| AMD (needs AMD hardware) | `julia --project=test test/amdgpu/runtests.jl` |
| Docs | `julia --project=docs docs/make.jl` |

`<name>` is the file name of `test/<name>.jl` without the extension, for example `grid` or `io`.

## Running on a GPU

One source runs on all devices (ADR-0001).
The package has no GPU dependency.
Loading the GPU package activates the package extension of CompReS for that device.
Pass the backend to `Grid`; every array of the run then lives on the device.

NVIDIA:

```julia
using CUDA, CompReS
grid = Grid(; Lx = 2π, Ly = 2π, Nx = 2048, Ny = 2048, T = Float32, backend = CUDABackend())
```

AMD:

```julia
using AMDGPU, CompReS
grid = Grid(; Lx = 2π, Ly = 2π, Nx = 2048, Ny = 2048, T = Float32, backend = ROCBackend())
```

The rest of the script does not change.
The output files are the same on every backend.
The AMD test run is not part of CI and has not been run on hardware.

## Float64 on consumer GPUs

The development machine has an RTX 5070.
It runs Float64 at about 1/64 of its Float32 rate.
Timings from this card do not predict timings on a cluster GPU with full-rate Float64.
Use `T = Float32` for a quick GPU check and `T = Float64` for production runs on a data-centre GPU or on the CPU.

## Benchmarks

The benchmark suite times `step!` and measures the memory of a run (ADR-0010, ADR-0030).
It has its own environment, `benchmark/`.

```
julia --project=benchmark -t auto benchmark/run.jl
julia --project=benchmark -t auto benchmark/run.jl --targets=cuda:Float32 --models=compressible --N=1024,2048
```

The keys of `run.jl` are `targets` (`device:Type` pairs, `cpu` or `cuda`), `models`, `schemes`, `ifs`, `N` and `out`.
The default set is both models with RK4 and the integrating factor on `cpu:Float64`, `cuda:Float64` and `cuda:Float32`, for `N = 256 … 4096`.
The full sweep takes about 4 minutes on the development machine.
Run nothing else on the machine during a sweep.

The run fits `t_step = a·N² log N + b·N² + c` to each configuration and writes the points and the coefficients to `benchmark/results/<machine>.toml`. The machine label is `--machine=`, else the environment variable `COMPRES_MACHINE`, else the hostname; `predict.jl` takes the same option. The development machine uses `COMPRES_MACHINE=DesktopNvdia`.
The residual of a fit is the maximum over the points of `|fit − t| / t`.
The result file of the development machine is in the repository.
A fit does not transfer to another machine: run `run.jl` on the cluster node before you extrapolate.

`predict.jl` reads a result file and prints the time per step, the walltime and the memory of a target run:

```
julia --project=benchmark benchmark/predict.jl --N=2048 --model=compressible --device=cuda --T=Float32 \
    --t_end=100 --dt=1e-4
```

The prediction covers `step!` only; forcing, diagnostics and output add to it.
`benchmark/cpu_kernel.jl` compares one kernel under the KernelAbstractions `CPU()` backend with a `Threads.@threads` loop.
