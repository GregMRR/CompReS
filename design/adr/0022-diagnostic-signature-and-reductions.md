# Diagnostic signature, registry keys and shell reductions

Date: 2026-11-19. Status: accepted.

WP-8 fixes the parts of the diagnostics registry that ADR-0009 leaves open. Changes:

- **A diagnostic takes the problem (amends ADR-0009).** The signature is `f(state, problem, workspace)`, not `f(state, workspace)`. The dissipation rates need the model parameters (`ν`, `α`, `m`), and the fluxes need the nonlinear term of the model. The `problem` holds both; the argument order follows `rhs!`.
- **One name per quantity, one generic function per name.** The registry maps `(class, name)` to a generic function, for example `:E_kin => kinetic_energy`. Each model adds its method, so WP-11 adds compressible methods for the shared names and registers only its new names. `default_diagnostics(class, model)` gives the ADR-0009 defaults per model.
- **Classes are the symbols `:GLOBAL` and `:SPECTRAL`**, the keys of `outputs` in the target user API. FIELDS holds snapshots, not diagnostics.
- **A SPECTRAL diagnostic names its coordinate.** `:shell_k` (shell centres) for spectra, `:shell_edge` (upper shell edges `K_n`) for fluxes. The NetCDF writer (WP-12) reads it from the registry entry.
- **SPECTRAL values cover every shell of the grid** (`length(grid.shell_k)`, up to the grid corner), for both models. The compressible model needs these shells (ADR-0020); one length lets one file layout serve both models. In the incompressible model the shells outside the Galerkin truncation hold zero spectra and a constant flux.
- **Shell reduction without atomics.** The workspace holds the linear indices of all modes sorted by shell. One kernel thread per shell sums its modes in this fixed order. The result does not depend on thread scheduling, and the CPU and GPU sums agree. The cumulative flux sums the shell vector on the host after the copy; it has one entry per shell.
- **Stored injection rates are a `Dict{Symbol,T}` keyed by diagnostic name.** `record_injection!` writes `:inj_E_kin` and `:inj_Z` (WP-11: `:inj_Omega`); a name not yet stored reads as `NaN`. The store needs no per-model layout.
- **The fluxes use the advection also when `nonlinear = false`.** The flux then describes the nonlinear transfer of the current field, not a term of the evolved equations.
