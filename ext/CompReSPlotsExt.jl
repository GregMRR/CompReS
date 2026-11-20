# Package extension of CompReS, loaded when Plots is loaded: the Plots methods of the plotting functions
# (WP-15, ADR-0012, ADR-0029), the light option.
module CompReSPlotsExt

using Plots
using CompReS
using CompReS: output_file, snapshot_file, snapshot_field, color_range, read_spectra, record_index, log_values,
    read_global, panel_layout, animation_data

colormap(symmetric) = symmetric ? :balance : :viridis

# Colour map of `values` (size (x, y)) over the grid points `x`, `y`; Plots takes the matrix as (y, x).
function field_plot(x, y, values, (range, symmetric), title)
    return heatmap(x, y, permutedims(values); clims = range, c = colormap(symmetric), aspect_ratio = :equal,
                   xlims = (x[1], x[end]), ylims = (y[1], y[end]), xlabel = "x", ylabel = "y", title)
end

# Log–log plot of shell spectra with a vertical line at k_max.
function spectrum_plot(k_max, title)
    p = plot(; xscale = :log10, yscale = :log10, xlabel = "k", ylabel = "spectrum", title)
    vline!(p, [k_max]; color = :gray, linestyle = :dash, label = "k_max")
    return p
end

function CompReS.plot_field(::Val{:Plots}, path; field::Symbol = :omega)
    (; x, y, values, t) = snapshot_field(snapshot_file(path), field)
    return field_plot(x, y, values, color_range(extrema(values)...), "$field, t = $(round(t; sigdigits = 4))")
end

function CompReS.plot_spectrum(::Val{:Plots}, path; names = nothing, time::Union{Nothing,Real} = nothing)
    data = read_spectra(output_file(path, "spectral.nc"), names)
    i = record_index(data.time, time)
    p = spectrum_plot(data.k_max, "t = $(round(data.time[i]; sigdigits = 4))")
    for (name, values) in data.spectra
        plot!(p, data.k, log_values(values[:, i]); label = name)
    end
    return p
end

function CompReS.plot_global(::Val{:Plots}, path; names = nothing)
    data = read_global(output_file(path, "global.nc"), names)
    rows, columns = panel_layout(length(data.values))
    panels = [plot(data.time, values; xlabel = "t", title = name, legend = false) for (name, values) in data.values]
    return plot(panels...; layout = (rows, columns), size = (400 * columns, 250 * rows))
end

function CompReS.animate_snapshots(::Val{:Plots}, outdir, file; field::Symbol = :omega, spectrum = :E_kin_inc_k,
                                   framerate::Integer = 10, decades::Real = 12)
    data = animation_data(outdir, field, spectrum, decades)
    animation = Animation()
    for snapshot in data.files
        (; x, y, values, t) = snapshot_field(snapshot, field)
        p = field_plot(x, y, values, data.colors, "$field, t = $(round(t; sigdigits = 4))")
        if !isnothing(spectrum)
            s = data.spectra
            i = record_index(s.time, t)
            sp = spectrum_plot(s.k_max, "$spectrum, t = $(round(s.time[i]; sigdigits = 4))")
            plot!(sp, s.k, log_values(s.values[:, i]); label = string(spectrum), xlims = (s.k[1], 1.2 * s.k_max),
                  ylims = s.limits)
            p = plot(p, sp; layout = (1, 2), size = (1100, 500))
        end
        frame(animation, p)
    end
    Plots.buildanimation(animation, file, endswith(file, ".gif"); fps = framerate, show_msg = false)
    return file
end

end
