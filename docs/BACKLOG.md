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

### Share the transfer terms between the compressible fluxes and spectra
- Cost: `Pi_kin`, `Pi_kin_inc` and `Pi_kin_comp` each rebuild `ŵ` and `N̂^w` (14 FFTs), and `E_kin_inc_k`, `E_kin_comp_k`, `E_kin_inc`, `E_kin_comp` each rebuild `ŵ`. At N = 2048, Float64, the compressible SPECTRAL defaults take 741 ms on CPU and 102 ms on CUDA, about 4/5 of it in the three `Pi_kin` variants; the GLOBAL defaults take 491 ms and 69 ms.
- Change: compute `ŵ`, `N̂^w` and the physical `u`, `ρ`, `ω` once per sample, for example through a per-sample cache in the workspace, and let every diagnostic of that sample read them.
- Source: WP-11 report, `src/diagnostics/compressible.jl`.

### Cheaper compressible injection rates
- Cost: the compressible `record_injection!` takes 7 inverse FFTs per step: 66 ms on CPU and 9.6 ms on CUDA at N = 2048, Float64, about 10% of one RK4 step with IF (623 ms and 82 ms).
- Change: reuse the physical `u`, `ρ`, `ω` of the first stage of `rhs!`, which `step!` evaluates on the same state right after, or evaluate `⟨ρu·f⟩` on the forcing support in spectral space.
- Source: WP-11 report.

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

### Guard against a non-positive density
- Defect: `ρ ≤ 0` raises a `DomainError` from `log` on CPU and gives NaN on GPU.
- Change: check `minimum(ρ) > 0` together with the NaN check of WP-13, and stop with a clear message.
- Source: WP-10 review.

### Seed the OU forcing through the public random API
- Defect: `src/forcing.jl` calls the internal 4-word `Xoshiro` constructor. A Julia release can change it without notice.
- Change: derive the per-step generator from `(seed, step)` through public API only. Keep the forcing sequence reproducible across a restart.
- Source: WP-9 review.

## Output

### Write snapshots atomically
- Defect: `write_snapshot` writes straight to `snap_NNNNNN.nc`. A kill during the write leaves a partial file that a restart can pick up.
- Change: write to a temporary name in the same directory and rename after `close`.
- Source: WP-12 report, `src/io.jl`.

### Restart in another float type
- Defect: `read_snapshot!` expects a snapshot of the float type of the state; it does not check the attribute `float_type`.
- Change: compare `float_type` with the grid and either reject or convert with a warning.
- Source: WP-12 report.
