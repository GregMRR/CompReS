# NetCDF output of the classes GLOBAL, SPECTRAL and FIELDS, parameter files and restart (WP-12, ADR-0009, ADR-0024).

"""
    Output(; period, diagnostics = nothing)

The sampling setting of one output class: the sampling `period` in time units and, for the classes GLOBAL
and SPECTRAL, the names of the diagnostics to write. `diagnostics = nothing` selects the defaults of the
model ([`default_diagnostics`](@ref)); [`selected_diagnostics`](@ref) resolves it. The class FIELDS takes
no diagnostics. The constructor throws an `ArgumentError` for `period ≤ 0`.
"""
struct Output
    period::Float64
    diagnostics::Union{Nothing,Vector{Symbol}}
end

function Output(; period::Real, diagnostics = nothing)
    period > 0 || throw(ArgumentError("the output period must be positive, got $period"))
    return Output(Float64(period), isnothing(diagnostics) ? nothing : Symbol.(collect(diagnostics)))
end

"""
    selected_diagnostics(output::Output, class::Symbol, model) -> Vector{Pair{Symbol,Diagnostic}}

The diagnostics that `output` selects for the output class `class` (`:GLOBAL` or `:SPECTRAL`) of `model`,
as [`select_diagnostics`](@ref) returns them. The writers [`open_global`](@ref) and
[`open_spectral`](@ref) take this vector.
"""
selected_diagnostics(output::Output, class::Symbol, model) =
    select_diagnostics(class, something(output.diagnostics, default_diagnostics(class, model)))

# --- Parameters -----------------------------------------------------------------------------------------

# Nested table of the scalar parameters of the parameter struct `x`: numbers, booleans and strings, and the
# tables of nested parameter structs of this package. Arrays, the backend and the names in `skip` are left out.
function parameter_table(x; skip = ())
    table = Dict{String,Any}("type" => string(nameof(typeof(x))))
    for name in fieldnames(typeof(x))
        name in skip && continue
        value = getfield(x, name)
        if value isa Union{Number,Bool,AbstractString}
            table[string(name)] = value
        elseif value isa Symbol
            table[string(name)] = string(value)
        elseif parentmodule(typeof(value)) === CompReS && !(value isa AbstractArray)
            table[string(name)] = parameter_table(value)
        end
    end
    return table
end

"""
    parameters(problem; run = (;)) -> Dict{String,Any}

The parameters of `problem` as a nested table with the sections `float_type`, `model` (including its
equation of state), `grid` and `forcing`, and `run` when the `NamedTuple` `run` is not empty (the settings of
the driver, for example the scheme and the time step; its values are numbers or strings). Each parameter
struct gives its scalar fields and its name under `type`; arrays and the state of the forcing are left out.
[`write_parameters`](@ref) and the NetCDF writers write this table.
"""
function parameters(problem; run::NamedTuple = (;))
    table = Dict{String,Any}("float_type" => string(typeof(problem.grid.Lx)),
                             "model" => parameter_table(problem.model),
                             "grid" => parameter_table(problem.grid),
                             "forcing" => parameter_table(problem.forcing; skip = keys(forcing_state(problem.forcing))))
    isempty(run) || (table["run"] = Dict{String,Any}(string(k) => v for (k, v) in pairs(run)))
    return table
end

"""
    write_parameters(path, problem; run = (;)) -> path

Write [`parameters`](@ref)`(problem; run)` to the TOML file `path`, the human-readable copy of ADR-0009.
"""
function write_parameters(path, problem; run::NamedTuple = (;))
    open(path, "w") do io
        TOML.print(io, parameters(problem; run); sorted = true)
    end
    return path
end

# NetCDF attribute value of a parameter: booleans become 0 or 1.
attribute_value(v::Bool) = Int8(v)
attribute_value(v) = v

# Write the nested table as global attributes; the attribute name joins the keys with "_".
function set_attributes!(ds, table, prefix = "")
    for (key, value) in table
        name = isempty(prefix) ? key : string(prefix, "_", key)
        value isa AbstractDict ? set_attributes!(ds, value, name) : (ds.attrib[name] = attribute_value(value))
    end
    return ds
end

# New NetCDF file at `path` with the parameters of `problem` as global attributes, plus the class and period.
function create_file(path, problem, class, period, run)
    ds = NCDataset(path, "c")
    set_attributes!(ds, parameters(problem; run))
    ds.attrib["output_class"] = String(class)
    isnothing(period) || (ds.attrib["period"] = Float64(period))
    return ds
end

# --- GLOBAL and SPECTRAL files --------------------------------------------------------------------------

"""
    open_global(path, problem, diagnostics; period, run = (;)) -> NCDataset

Create the NetCDF file `path` of the output class GLOBAL and return it open. `diagnostics` is the vector
that [`selected_diagnostics`](@ref) returns. The file has the unlimited dimension `time` and, along it, the
variables `time` (the time `t`), `step` (the step number `n`) and one variable per diagnostic, in the float
type `T` of the grid. The global attributes hold the [`parameters`](@ref) with the keys of each table joined
by `_` (for example `model_ν`, `forcing_seed`), `output_class` and `period`.

[`write_global!`](@ref) appends one record per call; the caller closes the file with `close`.
"""
function open_global(path, problem, diagnostics; period::Real, run::NamedTuple = (;))
    T = typeof(problem.grid.Lx)
    ds = create_file(path, problem, :GLOBAL, period, run)
    defDim(ds, "time", Inf)
    defVar(ds, "time", T, ("time",))
    defVar(ds, "step", Int64, ("time",))
    for (name, _) in diagnostics
        defVar(ds, String(name), T, ("time",))
    end
    return ds
end

"""
    write_global!(ds, state, problem, workspace, diagnostics) -> ds

Append one record of the GLOBAL file `ds` of [`open_global`](@ref): the time and step of `state` and the
value of every diagnostic in `diagnostics`, called as `d(state, problem, workspace)`. The values
that a diagnostic reads from the workspace (the injection rates, ADR-0016) are `NaN` before the first step;
the file stores them as plain values. The function flushes the file to disk.
"""
function write_global!(ds, state, problem, workspace, diagnostics)
    i = ds.dim["time"] + 1
    ds["time"][i] = state.t
    ds["step"][i] = state.n
    for (name, d) in diagnostics
        ds[String(name)][i] = d(state, problem, workspace)
    end
    sync(ds)
    return ds
end

"""
    open_spectral(path, problem, diagnostics; period, run = (;)) -> NCDataset

Create the NetCDF file `path` of the output class SPECTRAL and return it open. The file has the unlimited
dimension `time` and the dimension `shell`, with one entry per shell of the grid. The variables are `time` and
`step` as in [`open_global`](@ref); the shell coordinates `shell_k` (centres) and `shell_edge` (upper edges) in
the type `T` of the grid; the mode counts `shell_count` (retained modes) and `shell_count_all` (all grid
modes); and one variable `(shell, time)` per diagnostic. The attribute `coordinates` of a diagnostic variable
names its coordinate variable, as the [`Diagnostic`](@ref) fixes it. The global attributes are those of
[`open_global`](@ref).

[`write_spectral!`](@ref) appends one record per call; the caller closes the file with `close`.
"""
function open_spectral(path, problem, diagnostics; period::Real, run::NamedTuple = (;))
    grid = problem.grid
    T = typeof(grid.Lx)
    ds = create_file(path, problem, :SPECTRAL, period, run)
    defDim(ds, "time", Inf)
    defDim(ds, "shell", length(grid.shell_k))
    defVar(ds, "time", T, ("time",))
    defVar(ds, "step", Int64, ("time",))
    defVar(ds, "shell_k", T, ("shell",))[:] = grid.shell_k
    defVar(ds, "shell_edge", T, ("shell",))[:] = grid.shell_edge
    defVar(ds, "shell_count", Int64, ("shell",))[:] = grid.shell_count
    defVar(ds, "shell_count_all", Int64, ("shell",))[:] = grid.shell_count_all
    for (name, d) in diagnostics
        v = defVar(ds, String(name), T, ("shell", "time"))
        v.attrib["coordinates"] = String(d.coordinate)
    end
    return ds
end

"""
    write_spectral!(ds, state, problem, workspace, diagnostics) -> ds

Append one record to the SPECTRAL file `ds` of [`open_spectral`](@ref): the time and step of `state` and the
shell vector of every diagnostic in `diagnostics`. The function flushes the file to disk.
"""
function write_spectral!(ds, state, problem, workspace, diagnostics)
    i = ds.dim["time"] + 1
    ds["time"][i] = state.t
    ds["step"][i] = state.n
    for (name, d) in diagnostics
        ds[String(name)][:, i] = d(state, problem, workspace)
    end
    sync(ds)
    return ds
end

# --- FIELDS snapshots and restart -----------------------------------------------------------------------

# Forcing state in the NetCDF file: a complex vector becomes the variable (complex, index) of real parts and
# imaginary parts; a real scalar becomes a scalar variable.
function write_forcing_value!(ds, name, z::AbstractVector{<:Complex{T}}) where {T}
    haskey(ds.dim, "complex") || defDim(ds, "complex", 2)
    defDim(ds, name * "_index", length(z))
    defVar(ds, name, T, ("complex", name * "_index"))[:, :] = Array(reinterpret(reshape, T, z))
end

write_forcing_value!(ds, name, x::AbstractFloat) = (defVar(ds, name, typeof(x), ())[] = x)

read_forcing_value(ds, name, ::AbstractVector{Complex{T}}) where {T} =
    Vector(reinterpret(reshape, Complex{T}, Array{T}(ds[name][:, :])))

read_forcing_value(ds, name, ::T) where {T<:AbstractFloat} = T(ds[name][])

"""
    write_snapshot(dir, index, state, problem, workspace; run = (;)) -> path

Write the FIELDS snapshot number `index` of `state` to `dir/snap_NNNNNN.nc` (`NNNNNN` is `index`, zero-padded
to six digits) and return the path. `dir` must exist. The snapshot is the only restart input (ADR-0009).

The dimensions are `x` (`Nx`) and `y` (`Ny`). The variables, in the float type `T` of the grid, are the
coordinates `x` and `y`; the physical-space fields of the state, `ux`, `uy` and, compressible model only,
`rho`, of size `(x, y)`; the time `t`; the step `n` (`Int64`); and the [`forcing_state`](@ref) of
`problem.forcing` as the variables `forcing_<name>`. A complex vector is stored with the dimension
`complex` of size 2 (real part, imaginary part). The attributes are those of [`open_global`](@ref) with
`output_class = "FIELDS"` and no `period`, plus `snapshot_index = index` (`Int64`), which [`restart!`](@ref)
reads to continue the numbering. The function overwrites `workspace.physical`.
"""
function write_snapshot(dir, index::Integer, state, problem, workspace; run::NamedTuple = (;))
    grid = problem.grid
    T = typeof(grid.Lx)
    path = joinpath(dir, "snap_$(lpad(index, 6, '0')).nc")
    ds = create_file(path, problem, :FIELDS, nothing, run)
    ds.attrib["snapshot_index"] = Int64(index)
    try
        defDim(ds, "x", grid.Nx)
        defDim(ds, "y", grid.Ny)
        defVar(ds, "x", T, ("x",))[:] = T(grid.Lx) .* (0:(grid.Nx - 1)) ./ grid.Nx
        defVar(ds, "y", T, ("y",))[:] = T(grid.Ly) .* (0:(grid.Ny - 1)) ./ grid.Ny
        for name in keys(state.fields)
            f = workspace.physical[name]
            to_physical!(f, state.fields[name], workspace.transforms)
            defVar(ds, String(name), T, ("x", "y"))[:, :] = Array(f)
        end
        defVar(ds, "t", T, ())[] = state.t
        defVar(ds, "n", Int64, ())[] = state.n
        for (name, value) in pairs(forcing_state(problem.forcing))
            write_forcing_value!(ds, "forcing_$name", value)
        end
    finally
        close(ds)
    end
    return path
end

"""
    read_snapshot!(state, problem, workspace, path; upscale = false) -> state

Restore `state` and the forcing state of `problem.forcing` from the snapshot `path` of [`write_snapshot`](@ref).
The function transforms the physical fields to spectral space and truncates them to the Galerkin truncation,
sets `state.t` and `state.n`, and calls [`set_forcing_state!`](@ref). The restored fields equal the saved ones
to `‖û_restart − û‖₂/‖û‖₂ ≤ 10·eps(T)` (ADR-0014, ADR-0015); with the forcing state, the seed and the step,
the forcing sequence continues exactly (ADR-0007). `problem` must have the grid size, the fields and the
forcing type (and seed) of the run that wrote the file. The function overwrites `workspace.physical`, and
throws an `ArgumentError` otherwise.

**Upscaling** (`upscale = true`, ADR-0031). The snapshot may come from another grid size of the same box
(`Lx`, `Ly`) and another float type. The function copies each Fourier mode
``\\hat f(\\mathbf{k}) = \\langle f\\, e^{-i\\mathbf{k}\\cdot\\mathbf{x}} \\rangle`` of the snapshot to the same
wavevector of `problem.grid`, sets every other mode to zero (Fourier padding), and applies [`truncate!`](@ref)
and [`enforce_hermitian!`](@ref). A refinement keeps every mode; a coarsening removes the modes outside the
new Galerkin truncation. The forcing rules change:

- another seed: a warning; the function restores the forcing state, and the run continues with the new seed;
- another forcing type or another number of forced modes: a warning; the forcing keeps its current state.

A support of the same type and mode count but with other wavevectors is not detected, because the
snapshot stores no wavevectors. The function throws an `ArgumentError` for another box or a missing field.
"""
function read_snapshot!(state, problem, workspace, path; upscale::Bool = false)
    grid, forcing = problem.grid, problem.forcing
    NCDataset(path) do ds
        if !upscale
            saved = (Nx = ds.attrib["grid_Nx"], Ny = ds.attrib["grid_Ny"], forcing = ds.attrib["forcing_type"])
            saved == (Nx = grid.Nx, Ny = grid.Ny, forcing = string(nameof(typeof(forcing)))) ||
                throw(ArgumentError("the snapshot $path was written for another grid size or forcing type: $saved"))
            hasproperty(forcing, :seed) && ds.attrib["forcing_seed"] != forcing.seed &&
                throw(ArgumentError("the snapshot $path has the forcing seed $(ds.attrib["forcing_seed"]), the forcing $(forcing.seed)"))
        else
            box = (Lx = Float64(ds.attrib["grid_Lx"]), Ly = Float64(ds.attrib["grid_Ly"]))
            isapprox(box.Lx, grid.Lx; rtol = 1e-6) && isapprox(box.Ly, grid.Ly; rtol = 1e-6) ||
                throw(ArgumentError("the snapshot $path has the box $box; upscaling keeps the box ($(grid.Lx), $(grid.Ly))"))
        end
        all(name -> haskey(ds, String(name)), keys(state.fields)) ||
            throw(ArgumentError("the snapshot $path lacks a field of $(keys(state.fields))"))
        for name in keys(state.fields)
            f̂ = state.fields[name]
            if upscale
                copyto!(f̂, padded_spectrum(Array(ds[String(name)][:, :]), grid))
                truncate!(f̂, grid)
                enforce_hermitian!(f̂, grid)
            else
                copyto!(workspace.physical[name], Array(ds[String(name)][:, :]))
                to_spectral!(f̂, workspace.physical[name], workspace.transforms)
                truncate!(f̂, grid)
            end
        end
        state.t = ds["t"][]
        state.n = ds["n"][]
        upscale ? restore_forcing_upscaled!(forcing, ds, path) : set_forcing_state!(forcing, read_forcing_state(ds, forcing))
    end
    return state
end

# The forcing state that the snapshot ds holds for the entries of forcing_state(forcing).
function read_forcing_state(ds, forcing)
    template = forcing_state(forcing)
    return NamedTuple{keys(template)}(map((k, v) -> read_forcing_value(ds, "forcing_$k", v), keys(template), values(template)))
end

# Forcing state after an upscaling read (ADR-0031): restored when the type and the number of forced modes
# match, with a warning for another seed; otherwise a warning, and the forcing keeps its state.
function restore_forcing_upscaled!(forcing, ds, path)
    saved_type, type = ds.attrib["forcing_type"], string(nameof(typeof(forcing)))
    if saved_type != type
        @warn "the snapshot $path has the forcing type $saved_type, the run $type; the forcing keeps its initial state"
        return forcing
    end
    hasproperty(forcing, :seed) && ds.attrib["forcing_seed"] != forcing.seed &&
        @warn "the snapshot $path has the forcing seed $(ds.attrib["forcing_seed"]), the run $(forcing.seed); " *
              "the forcing state is restored and the run continues with the new seed"
    saved = read_forcing_state(ds, forcing)
    if map(length, saved) != map(length, forcing_state(forcing))
        @warn "the snapshot $path has another number of forced modes than the run; the forcing keeps its initial state"
        return forcing
    end
    return set_forcing_state!(forcing, saved)
end

# Host spectral field on `grid` of the physical field f of another grid size of the same box: the modes of
# f that both grids store, without the Nyquist column and row of either grid; zero elsewhere (Fourier padding).
function padded_spectrum(f::AbstractMatrix, grid::Grid{T}) where {T}
    Mx, My = size(f)
    source = FFTW.rfft(T.(f)) ./ T(Mx * My)
    f̂ = zeros(Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny)
    i_max, j_max = min(Mx, grid.Nx) ÷ 2 - 1, min(My, grid.Ny) ÷ 2 - 1
    for j in -j_max:j_max, i in 0:i_max
        f̂[i + 1, mod(j, grid.Ny) + 1] = source[i + 1, mod(j, My) + 1]
    end
    return f̂
end
