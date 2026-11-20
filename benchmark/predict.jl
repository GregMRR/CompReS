# Predict the cost of a target run from a benchmark result file (WP-16, ADR-0010, ADR-0030).
#
#   julia --project=benchmark benchmark/predict.jl --N=2048 --model=compressible --device=cuda --T=Float32 \
#       --t_end=100 --dt=1e-4 [--scheme=RK4] [--integrating_factor=true] [--file=benchmark/results/<hostname>.toml]
#
# The prediction covers `step!` only: the output, the diagnostics and the forcing of a real run add to it.
# A fit made on one machine does not transfer to another one (ADR-0010).
using TOML

function parse_arguments(args)
    options = Dict("N" => nothing, "model" => nothing, "device" => nothing, "T" => nothing, "t_end" => nothing,
                   "dt" => nothing, "scheme" => "RK4", "integrating_factor" => "true",
                   "file" => joinpath(@__DIR__, "results", gethostname() * ".toml"))
    for arg in args
        m = match(r"^--(\w+)=(.*)$", arg)
        (m !== nothing && haskey(options, m[1])) || error("unknown argument $arg; keys: $(join(keys(options), ", "))")
        options[m[1]] = m[2]
    end
    missing_keys = [k for (k, v) in options if v === nothing]
    isempty(missing_keys) || error("missing arguments: ", join(sort(missing_keys), ", "))
    return options
end

function main(args)
    o = parse_arguments(args)
    isfile(o["file"]) || error("no result file $(o["file"]); run benchmark/run.jl first")
    results = TOML.parsefile(o["file"])
    key = join((o["model"], o["scheme"], parse(Bool, o["integrating_factor"]) ? "IF" : "noIF", o["device"], o["T"]), "-")
    haskey(results["config"], key) || error("no configuration $key in $(o["file"]); it has: ",
                                           join(sort(collect(keys(results["config"]))), ", "))
    c = results["config"][key]
    N, t_end, Δt = parse(Int, o["N"]), parse(Float64, o["t_end"]), parse(Float64, o["dt"])
    lo, hi = extrema(c["N"])
    lo <= N <= hi || @warn "N = $N lies outside the measured range $lo to $hi; the prediction extrapolates"
    t_step = c["fit_a"] * N^2 * log(N) + c["fit_b"] * N^2 + c["fit_c"]
    memory = c["memory_fit_m"] * N^2 + c["memory_fit_c"]
    steps = ceil(Int, t_end / Δt)
    println("configuration  ", key, "  (measured on ", results["meta"]["hostname"], ", ", results["meta"]["date"], ")")
    println("time per step  ", round(1e3 * t_step; sigdigits = 4), " ms")
    println("steps          ", steps)
    println("walltime       ", round(steps * t_step / 3600; sigdigits = 4), " h  (", round(steps * t_step; sigdigits = 4), " s)")
    println("memory         ", round(memory / 2^30; sigdigits = 4), " GiB on the ", c["device"] == "cuda" ? "device" : "host")
end

main(ARGS)
