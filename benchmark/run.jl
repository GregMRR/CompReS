# Benchmark of the time per step and the memory of CompReS (WP-16, ADR-0010, ADR-0030).
#
#   julia --project=benchmark -t auto benchmark/run.jl [--key=a,b,...] ...
#
# Keys (defaults in brackets):
#   --targets  device:Type pairs [cpu:Float64,cuda:Float64,cuda:Float32]
#   --models   incompressible, compressible [both]
#   --schemes  Euler, RK2, SSPRK3, RK4 [RK4]
#   --ifs      integrating factor, true or false [true]
#   --N        grid sizes [256,512,1024,2048,4096]
#   --machine  machine label [$COMPRES_MACHINE, else the hostname]
#   --out      result file [benchmark/results/<machine>.toml]
# The run adds its configurations to the result file and keeps the other configurations in it.
using CompReS, CUDA, KernelAbstractions, TOML
using Dates: now

const SCHEMES = Dict("Euler" => Euler(), "RK2" => RK2(), "SSPRK3" => SSPRK3(), "RK4" => RK4())
const FLOATS = Dict("Float64" => Float64, "Float32" => Float32)
const WARMUP_STEPS = 2
const TIME_BUDGET = 3.0   # seconds of measured steps per point
const MIN_REPS, MAX_REPS = 5, 100
const Δt = 1e-4

"Key of one configuration in the result file."
config_key(model, scheme, integrating_factor, device, T) =
    join((model, scheme, integrating_factor ? "IF" : "noIF", device, T), "-")

"""
    fit_step_time(N, t) -> (coefficients, residual)

Fit `t = a·N² log N + b·N² + c` (natural logarithm) to the points `(N, t)`. The fit minimises the
relative error: each row of the least-squares system is divided by `t`, so the small grids count as
much as the large ones. The residual is the maximum over the points of `|fit - t| / t`.
"""
function fit_step_time(N, t)
    A = [N .^ 2 .* log.(N) N .^ 2 ones(length(N))]
    coefficients = (A ./ t) \ ones(length(t))
    return coefficients, maximum(abs.(A * coefficients .- t) ./ t)
end

"Fit `bytes = m·N² + c` by least squares on the relative error; return `(m, c)` and the residual."
function fit_memory(N, bytes)
    A = [N .^ 2 ones(length(N))]
    coefficients = (A ./ bytes) \ ones(length(bytes))
    return coefficients, maximum(abs.(A * coefficients .- bytes) ./ bytes)
end

# Bytes in use by the run: free device memory on a GPU (this includes the FFT plan workspaces), the Julia heap on the CPU.
function used_memory(backend)
    GC.gc(true)
    if backend isa CUDABackend
        CUDA.reclaim()
        return Int(CUDA.total_memory() - CUDA.free_memory())
    end
    return Int(Base.gc_live_bytes())
end

function make_model(model, T)
    model == "compressible" && return Compressible(; eos = Isothermal(cs = 10), μ = 0.01, ζ = 0.005)
    return Incompressible(; ν = 0.01)
end

# Median time of `step!` in seconds and memory in bytes; `nothing` when the grid does not fit the device.
function measure(model, scheme, integrating_factor, backend, T, N)
    before = used_memory(backend)
    grid = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N, T, backend)
    sim = try
        Simulation(grid, make_model(model, T); scheme, integrating_factor, Δt, outdir = mktempname())
    catch err
        err isa OutOfGPUMemoryError || rethrow()
        return nothing
    end
    set_initial_condition!(sim, random_velocity(grid; spectrum = k -> k^4 * exp(-k^2 / 8), E_kin = 0.5, seed = 1))
    memory = used_memory(backend) - before
    step_once!() = (step!(sim.state, sim.problem, sim.workspace, sim.Δt); KernelAbstractions.synchronize(backend))
    for _ in 1:WARMUP_STEPS
        step_once!()
    end
    first_time = @elapsed step_once!()
    reps = clamp(round(Int, TIME_BUDGET / first_time), MIN_REPS, MAX_REPS)
    times = [(@elapsed step_once!()) for _ in 1:reps]
    sim = grid = nothing
    return (t_step = sort(times)[(reps + 1) ÷ 2], memory = memory)
end

mktempname() = joinpath(tempdir(), "compres_benchmark_unused")

function parse_arguments(args)
    options = Dict("targets" => "cpu:Float64,cuda:Float64,cuda:Float32", "models" => "incompressible,compressible",
                   "schemes" => "RK4", "ifs" => "true", "N" => "256,512,1024,2048,4096",
                   "machine" => get(ENV, "COMPRES_MACHINE", gethostname()), "out" => "")
    for arg in args
        m = match(r"^--(\w+)=(.*)$", arg)
        (m !== nothing && haskey(options, m[1])) || error("unknown argument $arg; keys: $(join(keys(options), ", "))")
        options[m[1]] = m[2]
    end
    isempty(options["out"]) && (options["out"] = joinpath(@__DIR__, "results", options["machine"] * ".toml"))
    return options
end

function main(args)
    options = parse_arguments(args)
    split_list(key) = split(options[key], ',')
    targets = [Tuple(split(s, ':')) for s in split_list("targets")]
    models, schemes = split_list("models"), split_list("schemes")
    ifs = parse.(Bool, split_list("ifs"))
    Ns = parse.(Int, split_list("N"))
    results = isfile(options["out"]) ? TOML.parsefile(options["out"]) : Dict{String,Any}()
    results["meta"] = Dict("machine" => options["machine"], "date" => string(now()), "julia" => string(VERSION),
                           "threads" => Threads.nthreads(), "cpu" => Sys.cpu_info()[1].model,
                           "gpu" => CUDA.functional() ? CUDA.name(CUDA.device()) : "none",
                           "gpu_memory_bytes" => CUDA.functional() ? Int(CUDA.total_memory()) : 0,
                           "step_time_unit" => "s", "N_log" => "natural")
    configs = get!(results, "config", Dict{String,Any}())
    for (device, Tname) in targets, model in models, scheme in schemes, integrating_factor in ifs
        device == "cuda" && !CUDA.functional() && (@warn "CUDA is not functional; skip $device"; continue)
        backend = device == "cuda" ? CUDABackend() : CPU()
        T = FLOATS[Tname]
        key = config_key(model, scheme, integrating_factor, device, Tname)
        points = Tuple{Int,Float64,Int}[]
        for N in Ns
            point = measure(model, SCHEMES[scheme], integrating_factor, backend, T, N)
            if point === nothing
                @warn "$key: N = $N does not fit the device memory; point dropped"
                continue
            end
            push!(points, (N, point.t_step, point.memory))
            println(rpad(key, 44), " N = ", lpad(N, 5), "  t_step = ", round(1e3 * point.t_step; sigdigits = 4),
                    " ms  memory = ", round(point.memory / 2^20; digits = 1), " MiB")
        end
        length(points) >= 3 || (@warn "$key: fewer than 3 points, no fit"; continue)
        N, t, bytes = (float.(getindex.(points, i)) for i in 1:3)
        (a, b, c), residual = fit_step_time(N, t)
        (m, m0), memory_residual = fit_memory(N, bytes)
        configs[key] = Dict("model" => model, "scheme" => scheme, "integrating_factor" => integrating_factor,
                            "device" => device, "T" => Tname, "N" => Int.(N), "t_step" => t,
                            "memory_bytes" => Int.(bytes), "fit_a" => a, "fit_b" => b, "fit_c" => c,
                            "fit_residual" => residual, "memory_fit_m" => m, "memory_fit_c" => m0,
                            "memory_fit_residual" => memory_residual)
        println(rpad(key, 44), " fit: a = ", a, ", b = ", b, ", c = ", c, ", residual = ", round(100residual; digits = 1), " %")
        mkpath(dirname(options["out"]))
        open(io -> TOML.print(io, results), options["out"], "w")
    end
    println("results in ", options["out"])
end

main(ARGS)
