# Draw the mean time per step against the grid size, one figure per model, CPU and GPU lines (ADR-0032).
#
#   julia --project=benchmark benchmark/plot_scaling.jl [--machine=label] [--file=...] [--outdir=docs/src/assets]
#
# Reads `benchmark/results/<machine>-scaling.toml` (written by `benchmark/scaling.jl`) and writes
# `<outdir>/scaling_<model>.png` for each model in the file.
using CairoMakie, TOML

function parse_arguments(args)
    options = Dict("machine" => get(ENV, "COMPRES_MACHINE", gethostname()), "file" => "",
                   "outdir" => joinpath(@__DIR__, "..", "docs", "src", "assets"))
    for arg in args
        m = match(r"^--(\w+)=(.*)$", arg)
        (m !== nothing && haskey(options, m[1])) || error("unknown argument $arg; keys: $(join(keys(options), ", "))")
        options[m[1]] = m[2]
    end
    isempty(options["file"]) && (options["file"] = joinpath(@__DIR__, "results", options["machine"] * "-scaling.toml"))
    return options
end

const DEVICE_LABEL = Dict("cpu" => "CPU", "cuda" => "GPU (CUDA)")

function main(args)
    o = parse_arguments(args)
    results = TOML.parsefile(o["file"])
    meta, configs = results["meta"], results["config"]
    mkpath(o["outdir"])
    for model in sort(unique(c["model"] for c in values(configs)))
        Ns = sort(unique(n for c in values(configs) if c["model"] == model for n in c["N"]))
        fig = Figure(size = (720, 500))
        ax = Axis(fig[1, 1]; xscale = log2, yscale = log10, xticks = (Ns, string.(Ns)),
                  xlabel = "resolution N (grid N × N)", ylabel = "mean time per step (s)",
                  title = "$model, RK4 with integrating factor",
                  subtitle = "$(meta["machine"]): $(meta["cpu"]) with $(meta["threads"]) threads; $(meta["gpu"])")
        for device in ("cpu", "cuda")
            matching = [c for c in values(configs) if c["model"] == model && c["device"] == device]
            isempty(matching) && continue
            c = only(matching)
            scatterlines!(ax, c["N"], c["t_step"]; label = "$(DEVICE_LABEL[device]), $(c["T"]), mean of $(c["steps"]) steps")
        end
        axislegend(ax; position = :lt)
        file = joinpath(o["outdir"], "scaling_$model.png")
        save(file, fig)
        println("wrote ", file)
    end
end

main(ARGS)
