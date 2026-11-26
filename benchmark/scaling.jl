# Mean time per step against the grid size, CPU and GPU, for the scaling plots of the docs (ADR-0032).
#
#   nice -n 10 julia --project=benchmark -t 32 benchmark/scaling.jl [--key=a,b,...] ...
#
# Keys (defaults in brackets):
#   --devices  cpu, cuda [cpu,cuda]
#   --models   incompressible, compressible [both]
#   --N        grid sizes [128,256,512,1024,2048,4096,6144]
#   --steps    measured steps per point [300]
#   --T        float type [Float64]
#   --machine  machine label [$COMPRES_MACHINE, else the hostname]
#   --out      result file [benchmark/results/<machine>-scaling.toml]
# A point is the wall time of `steps` consecutive calls of `step!` divided by `steps`, after 5 warm-up calls,
# with one device synchronisation at the end. An unrecorded point at N = 64 runs first in each configuration.
# The run saves the result file after each point.
# `benchmark/plot_scaling.jl` draws the result file.
using CompReS, CUDA, KernelAbstractions, TOML
using Dates: now

const FLOATS = Dict("Float64" => Float64, "Float32" => Float32)
const WARMUP_STEPS = 5
const Δt = 1e-4

function make_model(model)
    model == "compressible" && return Compressible(; eos = Isothermal(cs = 10), μ = 0.01, ζ = 0.005)
    return Incompressible(; ν = 0.01)
end

# Mean time of `step!` in seconds over `steps` calls; `nothing` when the grid does not fit the device.
function measure(model, backend, T, N, steps)
    GC.gc(true)
    backend isa CUDABackend && CUDA.reclaim()
    grid = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N, T, backend)
    sim = try
        sim = Simulation(grid, make_model(model); scheme = RK4(), integrating_factor = true, Δt,
                         outdir = joinpath(tempdir(), "compres_benchmark_unused"))
        set_initial_condition!(sim, random_velocity(grid; spectrum = k -> k^4 * exp(-k^2 / 8), E_kin = 0.5, seed = 1))
        sim
    catch err
        (err isa OutOfGPUMemoryError || err isa CUDA.CUFFT.CUFFTError) || rethrow()
        return nothing
    end
    run_steps(n) = (for _ in 1:n; step!(sim.state, sim.problem, sim.workspace, sim.Δt); end;
                    KernelAbstractions.synchronize(backend))
    run_steps(WARMUP_STEPS)
    return (@elapsed run_steps(steps)) / steps
end

function parse_arguments(args)
    options = Dict("devices" => "cpu,cuda", "models" => "incompressible,compressible",
                   "N" => "128,256,512,1024,2048,4096,6144", "steps" => "300", "T" => "Float64",
                   "machine" => get(ENV, "COMPRES_MACHINE", gethostname()), "out" => "")
    for arg in args
        m = match(r"^--(\w+)=(.*)$", arg)
        (m !== nothing && haskey(options, m[1])) || error("unknown argument $arg; keys: $(join(keys(options), ", "))")
        options[m[1]] = m[2]
    end
    isempty(options["out"]) && (options["out"] = joinpath(@__DIR__, "results", options["machine"] * "-scaling.toml"))
    return options
end

function main(args)
    options = parse_arguments(args)
    split_list(key) = split(options[key], ',')
    Ns, steps, Tname = parse.(Int, split_list("N")), parse(Int, options["steps"]), options["T"]
    results = isfile(options["out"]) ? TOML.parsefile(options["out"]) : Dict{String,Any}()
    configs = get!(results, "config", Dict{String,Any}())
    for device in split_list("devices"), model in split_list("models")
        device == "cuda" && !CUDA.functional() && (@warn "CUDA is not functional; skip $device"; continue)
        backend = device == "cuda" ? CUDABackend() : CPU()
        key = join((model, device, Tname), "-")
        configs[key] = Dict{String,Any}("model" => model, "device" => device, "T" => Tname, "steps" => steps,
                                        "N" => Int[], "t_step" => Float64[])
        measure(model, backend, FLOATS[Tname], 64, 20)   # unrecorded: compiles the kernels and initialises the device
        for N in Ns
            t = measure(model, backend, FLOATS[Tname], N, steps)
            if t === nothing
                @warn "$key: N = $N does not fit the device memory; point dropped"
                continue
            end
            push!(configs[key]["N"], N); push!(configs[key]["t_step"], t)
            println(rpad(key, 32), " N = ", lpad(N, 5), "  mean t_step = ", round(1e3t; sigdigits = 4), " ms")
            results["meta"] = Dict("machine" => options["machine"], "date" => string(now()), "julia" => string(VERSION),
                                   "threads" => Threads.nthreads(), "cpu" => Sys.cpu_info()[1].model,
                                   "cpu_logical_cores" => Sys.CPU_THREADS,
                                   "gpu" => CUDA.functional() ? CUDA.name(CUDA.device()) : "none",
                                   "scheme" => "RK4", "integrating_factor" => true, "step_time_unit" => "s")
            mkpath(dirname(options["out"]))
            open(io -> TOML.print(io, results), options["out"], "w")
        end
    end
    println("results in ", options["out"])
end

main(ARGS)
