# Backlog: possible future implementations

This file lists changes that no work package in `docs/ROADMAP.md` covers yet. Each entry names the cost or defect, the proposed change, and where it came from. When an entry becomes work, move it into a WP and delete it here.

## Performance

### Fuse the stage combinations of `step!`
- Cost: one RK4 step with IF at N = 2048, Float64, takes 260 ms on CPU and 33 ms on CUDA; the four `rhs!` calls take 176 ms and 27 ms of it. The rest comes from 13 separate add and propagate passes, the two host-synchronising reductions of `record_injection!` and the host-to-device copy of the OU noise.
- Change: fuse each add-and-propagate sequence into one kernel that reads the linear rates; generate the OU noise on the device.
- Source: WP-7 report, `src/timestepping.jl` (`apply_ops!`).

### Reduce the stage storage of `step!`
- Cost: the stepper holds one tendency set per stage plus one stage set: 336 MB for RK4 at N = 2048, Float64, incompressible; about 5.4 GB at N = 8192, close to the 12 GB of the development GPU.
- Change: accumulate the final combination during the stages (low-storage form), so finished tendencies are freed for reuse.
- Source: WP-7 report.

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

### Check the SSPRK3 integrating-factor growth factor
- Defect: with the IF, SSPRK3 has `c_3 < c_2`, so its third stage multiplies by `exp(σΔt/2)`. It overflows when `σΔt/2 > log(floatmax(T))` (88 for Float32).
- Change: `Simulation` (WP-13) checks `max σ·Δt` once at construction and rejects the combination.
- Source: WP-7 report.

### Test `step!` on the nonlinear compressible model
- Defect: `step!` is tested on the compressible model only with `nonlinear = false`.
- Change: add a short nonlinear compressible run with a budget check to WP-11 or WP-13.
- Source: WP-7 report.

### Guard against a non-positive density
- Defect: `ρ ≤ 0` raises a `DomainError` from `log` on CPU and gives NaN on GPU.
- Change: check `minimum(ρ) > 0` together with the NaN check of WP-13, and stop with a clear message.
- Source: WP-10 review.

### Seed the OU forcing through the public random API
- Defect: `src/forcing.jl` calls the internal 4-word `Xoshiro` constructor. A Julia release can change it without notice.
- Change: derive the per-step generator from `(seed, step)` through public API only. Keep the forcing sequence reproducible across a restart.
- Source: WP-9 review.
