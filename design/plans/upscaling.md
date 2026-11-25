# Upscaling restart

Item of `PLAN.md`. Decisions: ADR-0031. Specified with the user on 2026-11-25.

## Interface

`restart!(sim, path; upscale = true)` forwards to `read_snapshot!(state, problem, workspace, path; upscale = true)`. Without `upscale`, both functions behave as before: a snapshot of another grid size, forcing type or seed throws an `ArgumentError`.

## Fields

1. The snapshot and the grid of `sim` must have the same `Lx` and `Ly`. Otherwise the function throws an `ArgumentError`.
2. The function reads each physical field of the snapshot on the host, converts it to the float type `T` of `sim`, and applies a real FFT with the normalisation of `to_spectral!`.
3. It copies each mode `(i, j)` (`i ≥ 0`, signed `j`) that the snapshot grid stores and the new grid stores to the same wavevector on the new grid. Every other mode of the new grid is zero (Fourier padding). The Nyquist column and row of either grid are left out.
4. It copies the spectral fields to the device, applies `truncate!` (a smaller grid removes the modes outside its disc) and `enforce_hermitian!`.

The source fields lie inside the source disc, so a refinement keeps every mode, the energy and the mean density exactly. Upscaling covers both directions; a coarsening truncates.

## Time, step and forcing

- `t`, `n` and the snapshot numbering (`snapshot_index + 1`) continue from the source.
- The forcing state is restored when the forcing type is the same and the forcing-state length matches. The order of the stored forced modes does not depend on the grid size, so the state maps mode by mode.
- A different seed gives a warning; the function restores the forcing state, and the run continues with the new seed.
- A different forcing type, or a forcing state of another length, gives a warning; the forcing keeps the initial state of its constructor.
- Limit: a support of the same type and the same mode count but other wavevectors is not detected, because the snapshot stores no wavevectors.

## Output directory and provenance

- `restart!(...; upscale = true)` throws an `ArgumentError` when `sim.outdir` holds a NetCDF file (`global.nc`, `spectral.nc`, a snapshot) of another grid size. The SPECTRAL records of two grids have different shell counts.
- The run table (`run_parameters`) gets `upscaled_from` (the snapshot path), `upscaled_from_Nx` and `upscaled_from_Ny`. They appear in `params.toml`, in the attributes of every file (`run_upscaled_from`, ...) and in the header of `run.log`.
- A plain `restart!` from a snapshot that has these attributes keeps them. `set_initial_condition!` clears them.
- The new run's parameters (model, forcing, grid) are those of the new `Simulation`.

## Acceptance checks (blocking, CPU and CUDA)

1. Single mode, refinement: an `N = 8` snapshot of one Fourier mode, read on `N = 16`. The mode equals the source to `10·eps(T)`; every other mode is zero.
2. Round trip: `16 → 32 → snapshot → 16` gives the source spectral fields to `10·eps(T)`.
3. Coarsening `32 → 16`: the modes inside the `N = 16` disc are unchanged; a mode outside it is zero.
4. Invariants: refinement keeps `E_kin` (incompressible and compressible) and the mean of `rho`.
5. Guards: without `upscale`, a grid mismatch throws. With `upscale`, a different box throws; an `outdir` with output of another grid throws.
6. Forcing: a different seed warns and restores the forcing state; a different type warns and keeps the initial state.
7. Provenance and continuation: `run_upscaled_from` and the source grid in GLOBAL, SPECTRAL and `params.toml`; `t`, `n` and the next snapshot index continue.
8. Run: `run!` after an upscaling restart gives finite fields.
