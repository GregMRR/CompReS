# Output: NetCDF files, output classes, diagnostics registry

Date: 2026-11-19. Status: accepted.

Output is NetCDF (NCDatasets.jl). Each run directory holds:

| File | Content |
|---|---|
| `global.nc` | GLOBAL class, unlimited time dimension |
| `spectral.nc` | SPECTRAL class, unlimited time dimension, shell centres and mode counts |
| `fields/snap_NNNNNN.nc` | one FIELDS snapshot per file |
| `run.log` | runtime log (ADR-0010) |
| `params.toml` | human-readable copy of all parameters |

Each class has its own sampling period, stored as a file attribute. All parameters are global attributes of every `.nc` file. A snapshot holds the physical-space fields (`ux`, `uy` and, compressible only, `rho`), the forcing state, the seed and the step, and it is the only restart input. A restart reproduces the saved state to round-off (`‖û_restart − û‖₂/‖û‖₂ ≤ 10·eps(T)`), because the physical-to-spectral round trip is not exact in floating point, and it reproduces the forcing sequence exactly (amended by ADR-0014, ADR-0015).

A **diagnostic** is a named function `(state, workspace) → value` registered under one output class. Injection rates need the state at the start of a step, so `step!` calls the hook `record_injection!`, which stores them in the workspace, and their diagnostics read them (amended by ADR-0016, ADR-0017). A run selects diagnostics by name; a new diagnostic is one function and one registration line.

Naming: a quantity on the solenoidal part has suffix `_inc`, one on the compressive part `_comp` (`E_kin_inc`, `Pi_kin_comp`). Vorticity quantities (`Z`, `Omega`, `P`) have no split. Rates are `inj_<q>`, `dissS_<q>` (small scale) and `dissL_<q>` (large scale). The incompressible model writes the same names with the `_comp` ones absent, so analysis scripts work on both models.

Compressible GLOBAL quantities are grid means of the physical fields; compressible spectra sum over all grid modes, because `√ρ u` and `ω/√ρ` are not band-limited (amended by ADR-0020).

Default GLOBAL: incompressible `E_kin`, `Z`, `P` and the injection and dissipation rates of `E_kin` and `Z`. Compressible: mass, `E_kin`, `E_kin_inc`, `E_kin_comp`, `E_int`, `E_tot`, `Z`, `Omega`, `u_rms`, `Ma`, `p_theta` (`⟨pθ⟩`), and the rates of `E_kin` and `Omega`.

Default SPECTRAL: incompressible `E_kin_inc_k`, `Z_k`, `Pi_kin`, `Pi_Z`. Compressible: `E_kin_inc_k`, `E_kin_comp_k`, `Omega_k`, `rho_k`, `Pi_kin`, `Pi_kin_inc`, `Pi_kin_comp`, `Pi_Omega`.

A run is a short Julia script against the package API; there is no separate input-file format.
