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
