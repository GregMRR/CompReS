# Output format

Output is NetCDF, written with NCDatasets.jl (ADR-0009, ADR-0024, ADR-0027).
`outputs` of [`Simulation`](@ref) selects the output classes; a class that is absent writes no file.

## Run directory

| File | Content |
|---|---|
| `global.nc` | GLOBAL class: one value per diagnostic and sample time |
| `spectral.nc` | SPECTRAL class: one value per shell, diagnostic and sample time |
| `fields/snap_NNNNNN.nc` | one snapshot of the FIELDS class per file |
| `run.log` | runtime log: header, one progress line per snapshot, stop line, timing table |
| `params.toml` | all parameters, as nested tables |

## Sampling

Each class has its own period.
It samples at the start time `t₀` and at the step nearest to each multiple `jP` of its period `P`.
The multiples refer to `t = 0`, so a restart keeps the sampling times.

## `global.nc`

- Unlimited dimension `time`.
- Variables `time` and one `(time)` variable per diagnostic, all of the float type `T` of the grid.
- Variable `step` (`Int64`): the step number of each record.
- The writer flushes each record, so a run that stops at its walltime leaves complete records.
- The injection rates `inj_*` are NaN at `t₀`, because no step has run. The variables have no `_FillValue` attribute.

## `spectral.nc`

- Unlimited dimension `time` and fixed dimension `shell`, with one entry per grid shell.
- Coordinate variables `shell_k` (shell centres) and `shell_edge` (upper edges) of type `T`, and `shell_count` (retained modes) and `shell_count_all` (all grid modes) of type `Int64`.
- Each diagnostic is a `(shell, time)` variable. Its attribute `coordinates` names its shell coordinate, `shell_k` or `shell_edge`.
  Spectra use `shell_k`; fluxes `Pi_*` use `shell_edge`, because a flux belongs to a cutoff.

## Snapshots

- The file name is `snap_NNNNNN.nc`, where `NNNNNN` is the snapshot index. The index counts the FIELDS samples, walltime snapshots and failure snapshots of a run, starting at 0 at `t₀`. It is not the step number. The attribute `snapshot_index` stores it.
- Dimensions `x` and `y`; variables `x` and `y` (grid-point coordinates).
- Physical fields `ux`, `uy` and, in the compressible model, `rho`.
- Scalars `t` (type `T`) and `n` (`Int64`, the step).
- One variable `forcing_<name>` for each entry of the forcing state. A complex vector has the dimension `(complex = 2, index)`, with real parts and imaginary parts. The seed is the attribute `forcing_seed`.
- A snapshot is also the restart input. [`restart!`](@ref) rejects a file with another grid size, forcing type, forcing seed or field set.
  A restart reproduces the state to round-off and the forcing sequence exactly.

A restart writes the restored state again as the next snapshot index, because every run segment samples FIELDS at `t₀`.
`global.nc` and `spectral.nc` keep the records of steps before the restart step and drop the rest, so the time axis stays monotonic.

## Parameters

`params.toml` holds one nested table with the keys `float_type`, `model` (with `large_scale` and `eos`), `grid`, `forcing` and `run`.
The table `run` holds the driver settings: `scheme`, `integrating_factor`, `Δt`, `n_check`, `C_stop`, `max_walltime`, `profile` and `<CLASS>_period`.
It holds scalar values only; arrays and the forcing state are not in it.

Each `.nc` file carries the same table as global attributes.
The attribute name joins the keys with `_` (`model_eos_cs`, `run_Δt`, `forcing_seed`).
Booleans are stored as 0 or 1.
Each file also has the attributes `output_class` and `period` (not for FIELDS).

## Variable names

A quantity on the solenoidal part has the suffix `_inc`, on the compressive part `_comp` (`E_kin_inc`, `Pi_kin_comp`).
Vorticity quantities (`Z`, `Omega`, `P`) have no split.
Rates are `inj_<q>` (injection), `dissS_<q>` (small-scale dissipation) and `dissL_<q>` (large-scale dissipation); the compressible model also writes `kappa_<q>` for the density diffusion.
The incompressible model writes the same names without the `_comp` ones, so one analysis script reads both models.

Default GLOBAL variables:

- Incompressible: `E_kin`, `Z`, `P`, and the injection and dissipation rates of `E_kin` and `Z`.
- Compressible: `mass`, `E_kin`, `E_kin_inc`, `E_kin_comp`, `E_int`, `E_tot`, `Z`, `Omega`, `u_rms`, `Ma`, `p_theta`, and the rates of `E_kin` and `Omega`.

Default SPECTRAL variables:

- Incompressible: `E_kin_inc_k`, `Z_k`, `Pi_kin`, `Pi_Z`.
- Compressible: `E_kin_inc_k`, `E_kin_comp_k`, `Omega_k`, `rho_k`, `Pi_kin`, `Pi_kin_inc`, `Pi_kin_comp`, `Pi_Omega`.

`Output(period, diagnostics = nothing)` selects other names; `nothing` selects the defaults.
Compressible GLOBAL quantities are grid means of the physical fields.

## Reading the files

Julia, with NCDatasets.jl:

```julia
using NCDatasets
ds = NCDataset("runs/example/global.nc")
t, Ma = ds["time"][:], ds["Ma"][:]
close(ds)

NCDataset("runs/example/fields/snap_000010.nc") do ds
    ux = ds["ux"][:, :]          # size (Nx, Ny)
    ds.attrib["model_eos_cs"]    # a parameter
end
```

Python, with netCDF4:

```python
from netCDF4 import Dataset
with Dataset("runs/example/spectral.nc") as ds:
    k = ds["shell_k"][:]
    E = ds["E_kin_inc_k"][:]      # shape (time, shell) in Python order
    flux = ds["Pi_kin"][:]
```

Python, with xarray:

```python
import xarray as xr
g = xr.open_dataset("runs/example/global.nc")
g["E_kin"].plot()
snaps = xr.open_mfdataset("runs/example/fields/snap_*.nc", combine="nested", concat_dim="snapshot")
```

Julia reads a `(shell, time)` variable with `shell` as the first index.
Python reverses the dimension order, so the same variable has shape `(time, shell)`.
