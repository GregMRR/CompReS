# Upscaling restart

Date: 2026-11-25. Status: accepted.

A run can continue from a FIELDS snapshot on another grid size of the same box: `restart!(sim, path; upscale = true)`, which calls `read_snapshot!(...; upscale = true)`. Plan: `design/plans/upscaling.md`. This ADR amends ADR-0007 (exact continuation of the forcing sequence) and ADR-0024 (restart rules) for this case only; a plain restart keeps their rules.

- **Keyword, not a conversion tool.** The upscaling happens in memory during the restart. No converted snapshot file is written, so the run directories hold no extra files. Without the keyword, a snapshot of another grid size still throws.
- **Fourier padding.** The function applies a host real FFT to each physical field of the snapshot, with the normalisation of `to_spectral!`. It copies each mode that both grids store to the same wavevector, sets the other modes to zero, and applies `truncate!` and `enforce_hermitian!`. The Nyquist column and row of either grid are left out. A refinement keeps every mode, so it keeps the energy and the mean density exactly. A coarsening truncates to the new disc. The box (`Lx`, `Ly`) must be the same. The float type may change.
- **Host only.** The transform of the snapshot runs once per restart, on the host. At `N = 8192`, three Float64 fields take about 1.6 GB.
- **Forcing.** The order of the stored forced modes does not depend on the grid size, so the forcing state maps mode by mode. Another seed gives a warning, and the state is restored. Another forcing type or another number of forced modes gives a warning, and the forcing keeps its initial state. The forcing sequence of an upscaled run therefore does not continue the source sequence exactly. A support of the same type and mode count but with other wavevectors is not detected (Ideas entry in `PLAN.md`).
- **Continuation.** The time `t`, the step `n` and the snapshot numbering continue from the source.
- **Parameters and provenance.** The new `Simulation` sets every parameter; the output files record them. The run table gets `upscaled_from` (absolute snapshot path), `upscaled_from_Nx` and `upscaled_from_Ny`, so every file has the attributes `run_upscaled_from`, ... and `params.toml` has them under `run`. A plain restart from a snapshot of the upscaled run keeps them; `set_initial_condition!` clears them.
- **Output directory.** An upscaling restart throws an `ArgumentError` when `outdir` holds a NetCDF file of another grid size, because the kept SPECTRAL records would have another shell count.
