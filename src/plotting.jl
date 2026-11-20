# Plotting interface (WP-15, ADR-0012, ADR-0029): the generic plotting functions, whose methods the extensions
# CompReSMakieExt and CompReSPlotsExt add, and the readers of the output files that both extensions share.
# This file loads no plotting library.

# The plotting libraries and the extensions that implement the plotting functions for them.
const PLOTTING_EXTENSIONS = (Makie = :CompReSMakieExt, Plots = :CompReSPlotsExt)

# The library of the only loaded plotting extension.
function loaded_library()
    loaded = [lib for (lib, ext) in pairs(PLOTTING_EXTENSIONS) if !isnothing(Base.get_extension(@__MODULE__, ext))]
    isempty(loaded) && throw(ArgumentError("no plotting extension is loaded: load CairoMakie, GLMakie or Plots"))
    length(loaded) == 1 ||
        throw(ArgumentError("both plotting extensions are loaded: select one with library = :Makie or :Plots"))
    return only(loaded)
end

"""
    plot_field(path; field = :omega, library) -> figure

Plot one physical-space field of a FIELDS snapshot as a colour map over the grid points, with the time of the
snapshot in the title. `path` is a snapshot file `snap_NNNNNN.nc` or an output directory; for an output
directory the function takes the snapshot of the highest index in `path/fields`. `field` is a field of the
snapshot (`:ux`, `:uy`, and `:rho` in the compressible model) or `:omega`, the vorticity

```math
\\omega = \\partial_x u_y - \\partial_y u_x,
```

which the function computes spectrally from `ux` and `uy` on the CPU (specification, Conventions). A field with
values of both signs gets a colour range symmetric about zero.

The keyword `library` selects the plotting extension: `:Makie` (loaded by `using CairoMakie` or `using
GLMakie`) returns a `Makie.Figure`, `:Plots` (loaded by `using Plots`) returns a `Plots.Plot`. The default is
the only loaded extension; the function throws an `ArgumentError` when none or both are loaded.
"""
plot_field(path; library::Symbol = loaded_library(), kwargs...) = plot_field(Val(library), path; kwargs...)

"""
    plot_spectrum(path; names = nothing, time = nothing, library) -> figure

Plot shell spectra of a SPECTRAL file on log–log axes against the shell centres `shell_k`, with a vertical line
at the largest retained wavenumber ``k_{\\max}``. `path` is a file `spectral.nc` or an output directory. `names`
are the spectra to plot; `nothing` selects every variable of the file whose coordinate is `shell_k`. `time`
selects the record nearest to that time; `nothing` selects the last record. The axes leave out the shell
``k = 0`` and every value that is not positive.

The function throws an `ArgumentError` for a name that is not in the file or whose coordinate is `shell_edge`
(a flux, which takes both signs).

The keyword `library` selects the plotting extension: `:Makie` (loaded by `using CairoMakie` or `using
GLMakie`) returns a `Makie.Figure`, `:Plots` (loaded by `using Plots`) returns a `Plots.Plot`. The default is
the only loaded extension; the function throws an `ArgumentError` when none or both are loaded.
"""
plot_spectrum(path; library::Symbol = loaded_library(), kwargs...) = plot_spectrum(Val(library), path; kwargs...)

"""
    plot_global(path; names = nothing, library) -> figure

Plot GLOBAL diagnostics against time, one panel per diagnostic. `path` is a file `global.nc` or an output
directory. `names` are the diagnostics to plot; `nothing` selects every diagnostic of the file. The injection
rates are `NaN` at the start time of each run segment (ADR-0024); the plot leaves those points out.

The function throws an `ArgumentError` for a name that is not in the file.

The keyword `library` selects the plotting extension: `:Makie` (loaded by `using CairoMakie` or `using
GLMakie`) returns a `Makie.Figure`, `:Plots` (loaded by `using Plots`) returns a `Plots.Plot`. The default is
the only loaded extension; the function throws an `ArgumentError` when none or both are loaded.
"""
plot_global(path; library::Symbol = loaded_library(), kwargs...) = plot_global(Val(library), path; kwargs...)

"""
    animate_snapshots(outdir, file; field = :omega, spectrum = :E_kin_inc_k, framerate = 10, decades = 12,
                      library) -> file

Write an animation of the FIELDS snapshots of the output directory `outdir` to the file `file` and return
`file`. The file extension of `file` selects the format, for example `.mp4` or `.gif`. Each frame shows
`field` of one snapshot, as [`plot_field`](@ref) draws it, and, unless `spectrum` is `nothing`, the shell
spectrum `spectrum` of `outdir/spectral.nc` at the record nearest to the snapshot time, as
[`plot_spectrum`](@ref) draws it. The frames follow the snapshot index.

All frames share one colour range (the extrema of `field` over all snapshots) and one spectrum axis, which
spans `decades` decades below the largest spectrum value over all records, so the frames show the evolution of
the spectrum tail. `framerate` is in frames per second.

The keyword `library` selects the plotting extension, `:Makie` or `:Plots`, as in [`plot_field`](@ref).
"""
animate_snapshots(outdir, file; library::Symbol = loaded_library(), kwargs...) =
    animate_snapshots(Val(library), outdir, file; kwargs...)

# --- Readers shared by the extensions -------------------------------------------------------------------

# The file `name` of the output directory `path`, or `path` itself when it is a file.
output_file(path, name) = isdir(path) ? joinpath(path, name) : path

# The snapshot files of the output directory `dir`, in index order.
function snapshot_files(dir)
    fields = joinpath(dir, "fields")
    files = isdir(fields) ? sort!(filter(f -> startswith(f, "snap_") && endswith(f, ".nc"), readdir(fields))) : String[]
    isempty(files) && throw(ArgumentError("the output directory $dir holds no snapshot in $fields"))
    return joinpath.(fields, files)
end

# The snapshot file `path`, or the snapshot of the highest index when `path` is an output directory.
snapshot_file(path) = isdir(path) ? last(snapshot_files(path)) : path

# Vorticity of the physical velocity `ux`, `uy` of the snapshot `ds`, computed spectrally on the CPU.
function snapshot_vorticity(ds)
    ux, uy = Array(ds["ux"][:, :]), Array(ds["uy"][:, :])
    T = eltype(ux)
    g = Grid(; Lx = ds.attrib["grid_Lx"], Ly = ds.attrib["grid_Ly"], Nx = size(ux, 1), Ny = size(ux, 2), T)
    tr = Transforms(g)
    û, ω̂, ∂ = (zeros(Complex{T}, g.Nx ÷ 2 + 1, g.Ny) for _ in 1:3)
    ∂x!(ω̂, to_spectral!(û, uy, tr), tr)
    ω̂ .-= ∂y!(∂, to_spectral!(û, ux, tr), tr)
    return to_physical!(similar(ux), ω̂, tr)
end

# The field `field` of the snapshot `path` (a name of the file or :omega), its coordinates and time.
function snapshot_field(path, field::Symbol)
    return NCDataset(path) do ds
        if field !== :omega && !haskey(ds, String(field))
            throw(ArgumentError("the snapshot $path has no field $field; use :omega or one of $(keys(ds))"))
        end
        values = field === :omega ? snapshot_vorticity(ds) : Array(ds[String(field)][:, :])
        (x = Array(ds["x"][:]), y = Array(ds["y"][:]), values, t = ds["t"][])
    end
end

# Colour range of the extrema `(lo, hi)`: symmetric about zero when the values take both signs.
# The second value tells whether the range is symmetric.
function color_range(lo, hi)
    lo < 0 < hi && return (-max(-lo, hi), max(-lo, hi)), true
    return (lo, hi), false
end

# Spectra `names` of the SPECTRAL file `path` (nothing: every variable with the coordinate shell_k).
# Returns the record times, the shell centres without the shell k = 0, k_max, and name => (shell, time) matrix.
function read_spectra(path, names)
    return NCDataset(path) do ds
        spectra = [name for name in keys(ds) if get(ds[name].attrib, "coordinates", "") == "shell_k"]
        selected = isnothing(names) ? spectra : String.(collect(names))
        for name in selected
            haskey(ds, name) || throw(ArgumentError("$path has no variable $name; its spectra are $spectra"))
            name in spectra || throw(ArgumentError("$name in $path is not a shell spectrum (coordinate shell_k)"))
        end
        (time = Array(ds["time"][:]), k = ds["shell_k"][2:end], k_max = ds.attrib["grid_k_max"],
         spectra = [name => Array(ds[name][:, :])[2:end, :] for name in selected])
    end
end

# Index of the record nearest to `time` in `times`; the last record for `time = nothing`.
record_index(times, ::Nothing) = lastindex(times)
record_index(times, time::Real) = argmin(abs.(times .- time))

# Values for log axes: NaN replaces every value that is not positive.
log_values(v) = map(x -> x > 0 ? x : oftype(x, NaN), v)

# Diagnostics `names` of the GLOBAL file `path` (nothing: every diagnostic), as the times and name => values.
function read_global(path, names)
    return NCDataset(path) do ds
        diagnostics = [name for name in keys(ds) if name ∉ ("time", "step") && dimnames(ds[name]) == ("time",)]
        selected = isnothing(names) ? diagnostics : String.(collect(names))
        for name in selected
            name in diagnostics || throw(ArgumentError("$path has no diagnostic $name; it has $diagnostics"))
        end
        (time = Array(ds["time"][:]), values = [name => Array(ds[name][:]) for name in selected])
    end
end

# Rows and columns of a near-square layout of `n` panels.
function panel_layout(n)
    columns = ceil(Int, sqrt(n))
    return cld(n, columns), columns
end

# The data of the frames of animate_snapshots: the snapshot files, the shared colour range of `field` and,
# unless `spectrum` is nothing, the spectra and the shared spectrum-axis limits over `decades` decades.
function animation_data(outdir, field, spectrum, decades)
    files = snapshot_files(outdir)
    lo, hi = mapreduce(f -> extrema(snapshot_field(f, field).values), ((a, b), (c, d)) -> (min(a, c), max(b, d)),
                       files)
    colors = color_range(lo, hi)
    isnothing(spectrum) && return (; files, colors, spectra = nothing)
    spectra = read_spectra(output_file(outdir, "spectral.nc"), (spectrum,))
    values = last(only(spectra.spectra))
    top = maximum(values)
    limits = (top * 10.0^(-decades), 10 * top)
    return (; files, colors, spectra = (; spectra.time, spectra.k, spectra.k_max, values, limits))
end
