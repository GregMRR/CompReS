# Backlog: possible future implementations

This file lists changes that no work package in `docs/ROADMAP.md` covers yet. Each entry names the cost or defect, the proposed change, and where it came from. When an entry becomes work, move it into a WP and delete it here.

## Analysis

### Intermittency analysis (WP-18)
- Defect: the solver has no analysis of the FIELDS snapshots, which is the purpose of the project.
- Change: specify WP-18 in detail with the user (increments, structure functions, increment PDFs, their dependence on the Mach number), then implement it as an offline module.
- Source: the user postponed WP-18 on 2026-11-20; `docs/ROADMAP.md`, WP-18.

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

### Fit the time per step piecewise in `N`
- Cost: the fit `a·N² log N + b·N² + c` of `benchmark/run.jl` misses the cache and thread-overhead transitions. On the development machine the maximum relative error is 15.2% (CPU, compressible) and 15.5% (CUDA Float32, incompressible), against the 10% of WP-16. A minimax fit gives 11.1% and 11.7%, so the form itself is too stiff. The time per mode is not monotonic: 65 ns at N = 256 (thread overhead), 21 ns at 512, 78 ns at 4096 on the CPU; on CUDA Float32 it doubles from N = 2048 to 4096, where one real array (64 MiB) probably exceeds the L2 cache.
- Change: fit one regime per cache size, or interpolate the measured points in `log N` and keep the three-term form for extrapolation beyond the largest `N`.
- Source: WP-16 benchmark, `benchmark/run.jl`, `benchmark/results/`.

### Time the FFTs in the profile mode
- Cost: the profile mode of `run!` times the tasks of `step!` (`rhs!`, stage combination, forcing, injection), not the FFTs and the pointwise kernels inside `rhs!`. ADR-0010 lists FFT timings.
- Change: pass the profile timer to `rhs!` and time each transform and kernel with `@profiled`.
- Source: WP-13 report, `src/timestepping.jl`.

## Robustness

### Seed the OU forcing through the public random API
- Defect: `src/forcing.jl` calls the internal 4-word `Xoshiro` constructor. A Julia release can change it without notice.
- Change: derive the per-step generator from `(seed, step)` through public API only. Keep the forcing sequence reproducible across a restart.
- Source: WP-9 review.

### Reject a polytropic balanced density with a negative enthalpy
- Defect: for `Polytropic`, `balanced_density` needs `h + h₀ > 0` at every grid point. At a high Mach number the root find evaluates a negative enthalpy: `density_from_enthalpy` raises a `DomainError` on CPU and gives NaN on GPU.
- Change: check `minimum(h) + h₀ > 0` (and the lower bracket end) and throw an `ArgumentError` that names the Mach number limit.
- Source: WP-14 report.

## Output

### Write snapshots atomically
- Defect: `write_snapshot` writes straight to `snap_NNNNNN.nc`. A kill during the write leaves a partial file that a restart can pick up.
- Change: write to a temporary name in the same directory and rename after `close`.
- Source: WP-12 report, `src/io.jl`.

### Restart in another float type
- Defect: `read_snapshot!` expects a snapshot of the float type of the state; it does not check the attribute `float_type`.
- Change: compare `float_type` with the grid and either reject or convert with a warning.
- Source: WP-12 report.

### Keep the old records of a restart without loading them
- Cost: a continued run reads the kept records of `spectral.nc` into memory and writes them into a new file. At N = 8192 with 10⁴ samples and 8 diagnostics, this is about 3.7 GB.
- Change: copy the records in chunks of time, or truncate the unlimited dimension in place.
- Source: WP-13 report, `src/simulation.jl` (`kept_records`).

## Plotting

### Plot the fluxes
- Defect: `plot_spectrum` rejects the SPECTRAL variables with the coordinate `shell_edge` (`Pi_kin`, `Pi_Z`, `Pi_Omega`), because they take both signs and do not fit log–log axes.
- Change: add `plot_flux(path; names, time)` with a logarithmic `k` axis and a linear value axis, and an optional flux panel in `animate_snapshots`.
- Source: WP-15 report, ADR-0029.

### Animate large runs faster
- Cost: `animate_snapshots` reads every snapshot twice (once for the shared colour range, once for the frame), and `field = :omega` builds FFT plans and computes the vorticity on the CPU per snapshot. At N = 512 this is under 1 s per snapshot; at N = 8192 the reads and FFTs dominate.
- Change: store the field extrema in the snapshot attributes, or build the `Transforms` once per call.
- Source: WP-15 report, `src/plotting.jl`.
