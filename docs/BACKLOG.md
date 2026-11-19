# Backlog: possible future implementations

This file lists changes that no work package in `docs/ROADMAP.md` covers yet. Each entry names the cost or defect, the proposed change, and where it came from. When an entry becomes work, move it into a WP and delete it here.

## Performance

### Share the advection between `Pi_kin` and `Pi_Z`
- Cost: each flux diagnostic calls `advection!` again (5 FFTs). Writing both fluxes doubles the cost. At N = 2048, Float64: `Pi_kin` takes 51 ms on CPU and 7.6 ms on CUDA.
- Change: compute the advection once per SPECTRAL sample and let both fluxes read it, for example through a per-sample cache in the workspace.
- Source: WP-8 report, `src/diagnostics/incompressible.jl` (`cumulative_flux`).

### Skip the viscous FFTs in an inviscid compressible run
- Cost: the compressible `rhs!` always computes the viscous remainder. With `μ = ζ = 0`, it wastes 2 FFTs per call.
- Change: skip the viscous remainder when `μ = ζ = 0`.
- Source: WP-10 review, `src/nonlinear_compressible.jl`.

### Fuse the forward-FFT normalisation into the truncation
- Cost: the normalisation and the truncation are two passes over the same spectral array.
- Change: merge them into one kernel where a truncation follows a forward FFT.
- Source: WP-6 and WP-10 reviews, `src/transforms.jl`.

### Measure the CPU kernels under the KernelAbstractions POCL backend
- Cost: unknown. KernelAbstractions 0.10 runs `CPU()` kernels through `POCLBackend`.
- Change: WP-16 measures the kernel time against plain threaded loops. If POCL is slow, add a CPU path.
- Source: WP-1 to WP-10 test logs.

## Robustness

### Guard against a non-positive density
- Defect: `ρ ≤ 0` raises a `DomainError` from `log` on CPU and gives NaN on GPU.
- Change: check `minimum(ρ) > 0` together with the NaN check of WP-13, and stop with a clear message.
- Source: WP-10 review.

### Seed the OU forcing through the public random API
- Defect: `src/forcing.jl` calls the internal 4-word `Xoshiro` constructor. A Julia release can change it without notice.
- Change: derive the per-step generator from `(seed, step)` through public API only. Keep the forcing sequence reproducible across a restart.
- Source: WP-9 review.
