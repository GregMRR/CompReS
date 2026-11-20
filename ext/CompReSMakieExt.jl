# Package extension of CompReS, loaded when Makie is loaded (CairoMakie or GLMakie): the Makie methods of the
# plotting functions (WP-15, ADR-0012, ADR-0029). It uses only the Makie API.
module CompReSMakieExt

using Makie
using CompReS
using CompReS: output_file, snapshot_file, snapshot_field, color_range, read_spectra, record_index, log_values,
    read_global, panel_layout, animation_data

colormap(symmetric) = symmetric ? :balance : :viridis

# Colour map of `values` over the grid points `x`, `y` in the grid position `position`, with a colour bar.
# `values` and `title` may be observables.
function field_panel!(position, x, y, values, (range, symmetric), name, title)
    ax = Axis(position[1, 1]; xlabel = "x", ylabel = "y", aspect = DataAspect(), title)
    hm = heatmap!(ax, x, y, values; colorrange = range, colormap = colormap(symmetric))
    Colorbar(position[1, 2], hm; label = name)
    return ax
end

# Log–log axis of shell spectra in the grid position `position`, with a vertical line at k_max.
function spectrum_axis(position, k_max, title)
    ax = Axis(position; xlabel = "k", ylabel = "spectrum", xscale = log10, yscale = log10, title)
    vlines!(ax, [k_max]; color = :gray, linestyle = :dash, label = "k_max")
    return ax
end

function CompReS.plot_field(::Val{:Makie}, path; field::Symbol = :omega)
    (; x, y, values, t) = snapshot_field(snapshot_file(path), field)
    fig = Figure()
    field_panel!(fig[1, 1], x, y, values, color_range(extrema(values)...), string(field),
                 "$field, t = $(round(t; sigdigits = 4))")
    return fig
end

function CompReS.plot_spectrum(::Val{:Makie}, path; names = nothing, time::Union{Nothing,Real} = nothing)
    data = read_spectra(output_file(path, "spectral.nc"), names)
    i = record_index(data.time, time)
    fig = Figure()
    ax = spectrum_axis(fig[1, 1], data.k_max, "t = $(round(data.time[i]; sigdigits = 4))")
    for (name, values) in data.spectra
        lines!(ax, data.k, log_values(values[:, i]); label = name)
    end
    axislegend(ax)
    return fig
end

function CompReS.plot_global(::Val{:Makie}, path; names = nothing)
    data = read_global(output_file(path, "global.nc"), names)
    rows, columns = panel_layout(length(data.values))
    fig = Figure(; size = (400 * columns, 250 * rows))
    for (j, (name, values)) in enumerate(data.values)
        ax = Axis(fig[fldmod1(j, columns)...]; xlabel = "t", title = name)
        lines!(ax, data.time, values)
    end
    return fig
end

function CompReS.animate_snapshots(::Val{:Makie}, outdir, file; field::Symbol = :omega, spectrum = :E_kin_inc_k,
                                   framerate::Integer = 10, decades::Real = 12)
    data = animation_data(outdir, field, spectrum, decades)
    frame = Observable(snapshot_field(first(data.files), field))
    fig = Figure(; size = isnothing(spectrum) ? (600, 500) : (1100, 500))
    field_panel!(fig[1, 1], frame[].x, frame[].y, lift(f -> f.values, frame), data.colors, string(field),
                 lift(f -> "$field, t = $(round(f.t; sigdigits = 4))", frame))
    if !isnothing(spectrum)
        s = data.spectra
        i = lift(f -> record_index(s.time, f.t), frame)
        ax = spectrum_axis(fig[1, 2], s.k_max, lift(i -> "$spectrum, t = $(round(s.time[i]; sigdigits = 4))", i))
        limits!(ax, s.k[1], 1.2 * s.k_max, s.limits...)
        lines!(ax, s.k, lift(i -> log_values(s.values[:, i]), i))
    end
    record(fig, file, data.files; framerate) do snapshot
        frame[] = snapshot_field(snapshot, field)
    end
    return file
end

end
