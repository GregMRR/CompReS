# Simulation driver: set-up, initial condition, time loop with sampling, safeguards, runtime log and restart
# (WP-13, ADR-0009, ADR-0010, ADR-0027).

const OUTPUT_CLASSES = (:GLOBAL, :SPECTRAL, :FIELDS)

"""
    Simulation(grid::Grid{T}, model; forcing = NoForcing(), scheme = RK4(), integrating_factor = true, Δt,
               outputs = (;), outdir, max_walltime = Inf, n_check = 10, C_stop = 1, profile = false,
               overwrite = false)

A run of `model` ([`Incompressible`](@ref) or [`Compressible`](@ref)) on `grid`, driven by `forcing`, with the
time step `Δt`, the [`RungeKuttaScheme`](@ref) `scheme` and the integrating-factor switch
`integrating_factor`. [`set_initial_condition!`](@ref) or [`restart!`](@ref) sets the state, and
[`run!`](@ref) advances it and writes the output to the directory `outdir`.

The constructor converts the model (with its equation of state), the forcing, `Δt` and `C_stop` to the float
type `T` of the grid ([`convert_float`](@ref), ADR-0006). Without an initial condition, the state is at rest,
with the density ``\\rho = 1`` in the compressible model.

Keywords:

- `outputs`: a `NamedTuple` with one [`Output`](@ref) per sampled output class, under the keys `GLOBAL`,
  `SPECTRAL` and `FIELDS`. A class that is absent writes no file.
- `max_walltime`: the wall-clock time in seconds after which [`run!`](@ref) writes a snapshot and stops.
- `n_check`: the number of steps between two checks of the state (see [`run!`](@ref)).
- `C_stop`: the largest CFL number that the checks accept. The CFL number is
  ```math
  C = \\Delta t\\, \\max_{\\mathbf{x}} \\Big[\\frac{|u_x| + c}{\\Delta x} + \\frac{|u_y| + c}{\\Delta y}\\Big],
  ```
  with ``c = 0`` in the incompressible model and ``c = \\max_{\\mathbf{x}} c_s(\\rho)`` in the compressible
  model ([`sound_speed²`](@ref)), over the grid points.
- `profile`: `true` adds the times of the tasks of [`step!`](@ref) to the runtime log; on a GPU, the device
  synchronisations that this needs slow the run (ADR-0010).
- `overwrite`: `true` lets a new run delete the output files of another run in `outdir`.

The constructor throws an `ArgumentError` for `Δt ≤ 0`, `n_check < 1`, `C_stop ≤ 0`, `max_walltime < 0`, an
unknown output class or diagnostic name, a compressive forcing with the incompressible model
([`Problem`](@ref)), and for [`SSPRK3`](@ref) with the integrating factor when its growth factor
``e^{\\sigma\\Delta t/2}`` overflows the float type `T` on the largest linear rate ``\\sigma``. It warns when a
forced run uses [`Euler`](@ref) or [`RK2`](@ref): both are unstable on undamped oscillating modes.

The fields that a user reads are `grid`, `problem` (the [`Problem`](@ref)), `state` (the [`State`](@ref)),
`workspace` (the [`Workspace`](@ref)), `Δt`, `C_stop` and `cfl`, the CFL number of the last check.
"""
mutable struct Simulation{T<:AbstractFloat,G<:Grid,P<:Problem,S<:State,W<:Workspace,O<:NamedTuple,D<:NamedTuple}
    const grid::G
    const problem::P
    const state::S
    const workspace::W
    const Δt::T
    const outputs::O
    const diagnostics::D
    const outdir::String
    const max_walltime::Float64
    const n_check::Int
    const C_stop::T
    const profile::Bool
    const overwrite::Bool
    const timer::TimerOutput
    cfl::T
    new_run::Bool
    next_snapshot::Int
    snapshot_step::Int
    log_position::Int
end

# Spectral fields of the state at rest: zero velocity and, compressible model only, ρ = 1.
zero_velocity(grid::Grid{T}) where {T} =
    (ux = KernelAbstractions.zeros(grid.backend, Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny),
     uy = KernelAbstractions.zeros(grid.backend, Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny))
rest_fields(grid::Grid, ::Incompressible) = zero_velocity(grid)
rest_fields(grid::Grid, ::Compressible) = merge(zero_velocity(grid), uniform_density(grid))

# Largest linear rate σ of the model, over all fields and Helmholtz parts.
largest_rate(r::IncompressibleRates) = maximum(r.velocity)
largest_rate(r::CompressibleRates) = max(maximum(r.solenoidal), maximum(r.compressive), maximum(r.density))

# Reject a time stepper whose integrating factor multiplies by a growth factor exp(στ) (τ > 0) that overflows T.
# With the Horner form only a backward propagation (negative coefficient) grows; SSPRK3 has one, −Δt/2 (ADR-0023).
function check_growth_factor(stepper::TimeStepper, rates, Δt::T) where {T}
    stepper.integrating_factor || return nothing
    ops = Iterators.flatten((stepper.stage_ops..., stepper.final_ops))
    τ = maximum(op -> op.term == 0 ? -op.coeff * Δt : zero(T), ops)
    σ = largest_rate(rates)
    isfinite(exp(σ * τ)) ||
        throw(ArgumentError("$(nameof(typeof(stepper.scheme))) with the integrating factor multiplies by " *
                            "exp(στ) = exp($(σ * τ)), which overflows $T; reduce Δt, use RK4 or " *
                            "integrating_factor = false"))
    return nothing
end

function Simulation(grid::Grid{T}, model; forcing = NoForcing(), scheme::RungeKuttaScheme = RK4(),
                    integrating_factor::Bool = true, Δt::Real, outputs::NamedTuple = (;),
                    outdir::AbstractString, max_walltime::Real = Inf, n_check::Integer = 10, C_stop::Real = 1,
                    profile::Bool = false, overwrite::Bool = false) where {T}
    Δt > 0 || throw(ArgumentError("Δt must be positive, got $Δt"))
    n_check >= 1 || throw(ArgumentError("n_check must be at least 1, got $n_check"))
    C_stop > 0 || throw(ArgumentError("C_stop must be positive, got $C_stop"))
    max_walltime >= 0 || throw(ArgumentError("max_walltime must be nonnegative, got $max_walltime"))
    all(in(OUTPUT_CLASSES), keys(outputs)) ||
        throw(ArgumentError("the output classes are $OUTPUT_CLASSES, got $(keys(outputs))"))
    all(o -> o isa Output, outputs) || throw(ArgumentError("each output class takes an Output"))
    problem = Problem(grid, model, convert_float(T, forcing))
    state = State(rest_fields(grid, problem.model))
    timer = TimerOutput()
    ws = Workspace(grid, state; scheme, integrating_factor, timer = profile ? timer : nothing)
    check_growth_factor(ws.stepper, problem.rates, T(Δt))
    forcing isa NoForcing || !(scheme isa Union{Euler,RK2}) ||
        @warn "$(nameof(typeof(scheme))) is unstable on undamped oscillating modes; use SSPRK3 or RK4 for a forced run"
    classes = filter(class -> haskey(outputs, class), (:GLOBAL, :SPECTRAL))
    diagnostics = NamedTuple{classes}(map(class -> selected_diagnostics(outputs[class], class, problem.model), classes))
    return Simulation(grid, problem, state, ws, T(Δt), outputs, diagnostics, String(outdir), Float64(max_walltime),
                      Int(n_check), T(C_stop), profile, overwrite, timer, T(NaN), true, 0, -1, 0)
end

"""
    set_initial_condition!(sim::Simulation, fields::NamedTuple) -> sim

Copy the spectral `fields` into the state of `sim` and set `t = 0` and `n = 0`; [`run!`](@ref) then starts a
new run. `fields` holds `ux` and `uy` and, compressible model only, optionally `rho`; without `rho` the
density is uniform, ``\\rho = 1``. The initial-condition functions return such fields, for example
[`random_velocity`](@ref); `merge` combines a velocity with a density. The function throws an
`ArgumentError` for a missing velocity component, a field that the model does not evolve, or a field of
another size.
"""
function set_initial_condition!(sim::Simulation, fields::NamedTuple)
    names = keys(sim.state.fields)
    (haskey(fields, :ux) && haskey(fields, :uy)) ||
        throw(ArgumentError("the initial condition needs ux and uy, got $(keys(fields))"))
    all(in(names), keys(fields)) ||
        throw(ArgumentError("the model evolves the fields $names; the initial condition has $(keys(fields))"))
    all(f -> size(f) == size(sim.state.fields.ux), fields) ||
        throw(ArgumentError("the initial fields must have the spectral size $(size(sim.state.fields.ux))"))
    for name in names
        # Only rho can be absent.
        copyto!(sim.state.fields[name], haskey(fields, name) ? fields[name] : uniform_density(sim.grid).rho)
    end
    sim.state.t = 0
    sim.state.n = 0
    sim.new_run = true
    sim.next_snapshot = 0
    return sim
end

"""
    restart!(sim::Simulation, path) -> sim

Restore the state of `sim` and the state of its forcing from the snapshot `path` ([`read_snapshot!`](@ref));
[`run!`](@ref) then continues that run in `sim.outdir`. `sim` must be set up as the run that wrote the
snapshot (same grid size, model fields, forcing type and seed); `Δt` and the other driver settings may
differ. The snapshot numbering continues after the index of `path`. The function throws an
`ArgumentError` for a snapshot of another run.
"""
function restart!(sim::Simulation, path)
    read_snapshot!(sim.state, sim.problem, sim.workspace, path)
    sim.next_snapshot = NCDataset(ds -> Int(ds.attrib["snapshot_index"]) + 1, path)
    sim.new_run = false
    return sim
end

"""
    cfl_number(sim::Simulation) -> T

The CFL number of the state of `sim` (see [`Simulation`](@ref)), from the velocity and, compressible model,
the sound speed at the grid points. The density must be positive. The function overwrites
`sim.workspace.physical` and `sim.workspace.vorticity`.
"""
function cfl_number(sim::Simulation{T}) where {T}
    (; grid, state, workspace) = sim
    tr = workspace.transforms
    for name in keys(state.fields)
        to_physical!(workspace.physical[name], state.fields[name], tr)
    end
    Δx, Δy = T(grid.Lx) / grid.Nx, T(grid.Ly) / grid.Ny
    local_rate = workspace.vorticity
    local_rate .= abs.(workspace.physical.ux) ./ Δx .+ abs.(workspace.physical.uy) ./ Δy
    c = largest_sound_speed(sim.problem.model, workspace)
    return sim.Δt * (maximum(local_rate) + c * (1 / Δx + 1 / Δy))
end

largest_sound_speed(model::Incompressible{T}, ws) where {T} = zero(T)
largest_sound_speed(model::Compressible, ws) = sqrt(mapreduce(ρ -> sound_speed²(model.eos, ρ), max, ws.physical.rho))

# Check of the state every n_check steps: finite fields, positive density, CFL number at most C_stop. Return
# the failure message, or nothing.
function check_state(sim::Simulation)
    (; state, workspace) = sim
    all(f -> all(isfinite, f), state.fields) || return "the state holds NaN or Inf"
    if haskey(state.fields, :rho)
        ρ = workspace.physical.rho
        ρ_min = minimum(to_physical!(ρ, state.fields.rho, workspace.transforms))
        ρ_min > 0 || return "the density is not positive: min ρ = $ρ_min"
    end
    sim.cfl = cfl_number(sim)
    sim.cfl <= sim.C_stop || return "the CFL number $(sim.cfl) exceeds C_stop = $(sim.C_stop)"
    return nothing
end

# A diagnostic that times itself under its name in the timer of the run.
struct TimedDiagnostic{D}
    diagnostic::D
    timer::TimerOutput
    label::String
end

(d::TimedDiagnostic)(state, problem, ws) = @timeit d.timer d.label d.diagnostic(state, problem, ws)

# The driver settings written to params.toml and to the attributes of every file (ADR-0024).
run_parameters(sim::Simulation) =
    (; scheme = string(nameof(typeof(sim.workspace.stepper.scheme))),
     integrating_factor = sim.workspace.stepper.integrating_factor, Δt = sim.Δt, n_check = sim.n_check,
     C_stop = sim.C_stop, max_walltime = sim.max_walltime, profile = sim.profile,
     (Symbol(class, "_period") => output.period for (class, output) in pairs(sim.outputs))...)

fields_dir(sim) = joinpath(sim.outdir, "fields")
log_path(sim) = joinpath(sim.outdir, "run.log")
class_path(sim, class) = joinpath(sim.outdir, class === :GLOBAL ? "global.nc" : "spectral.nc")

# Output files of a run in outdir: the files that a new run would overwrite.
function run_files(sim)
    files = filter(isfile, [class_path(sim, :GLOBAL), class_path(sim, :SPECTRAL), log_path(sim)])
    isdir(fields_dir(sim)) &&
        append!(files, joinpath.(fields_dir(sim), filter(startswith("snap_"), readdir(fields_dir(sim)))))
    return files
end

# Records with step < n₀ of the time variables `names` of the file `path`, or nothing without the file.
function kept_records(path, names, n₀)
    isfile(path) || return nothing
    return NCDataset(path) do ds
        saved = Set(name for name in keys(ds) if "time" in dimnames(ds[name]))
        saved == Set(names) ||
            throw(ArgumentError("$path holds the variables $(sort!(collect(saved))); this run writes $(sort(names)). " *
                                "Use the same diagnostics or another outdir"))
        k = count(<(n₀), ds["step"][:])
        Dict(name => Array(ds[name][ntuple(_ -> Colon(), ndims(ds[name]) - 1)..., 1:k]) for name in names)
    end
end

# Open the file of the output class `class`. A continued run keeps the records before its first step n₀.
function open_class(sim, class, run)
    path = class_path(sim, class)
    diagnostics = sim.diagnostics[class]
    names = ["time"; "step"; [String(name) for (name, _) in diagnostics]]
    kept = sim.new_run ? nothing : kept_records(path, names, sim.state.n)
    open_file = class === :GLOBAL ? open_global : open_spectral
    ds = open_file(path, sim.problem, diagnostics; period = sim.outputs[class].period, run)
    if !isnothing(kept) && !isempty(kept["step"])
        for name in names
            ds[name][ntuple(_ -> Colon(), ndims(ds[name]) - 1)..., 1:length(kept["step"])] = kept[name]
        end
    end
    return ds
end

# Append text to run.log and rewrite the timer table after it; the next call overwrites the table.
function write_log!(sim::Simulation, text = "")
    isfile(log_path(sim)) || touch(log_path(sim))
    open(log_path(sim), "r+") do io
        seek(io, sim.log_position)
        print(io, text)
        sim.log_position = position(io)
        print_timer(IOContext(io, :displaysize => (24, 160)), sim.timer)
        println(io)
        truncate(io, position(io))
    end
    return nothing
end

function progress_line(sim, start)
    (; t, n) = sim.state
    return "t = $t  n = $n  C = $(round(sim.cfl; sigdigits = 3))  wall = $(round(time() - start; digits = 1)) s\n"
end

function write_header!(sim::Simulation, t_end)
    mkpath(sim.outdir)
    sim.log_position = isfile(log_path(sim)) ? filesize(log_path(sim)) : 0
    (; grid, problem, state) = sim
    header = """
        == run segment started $(Libc.strftime("%Y-%m-%d %H:%M:%S", time())): t = $(state.t), n = $(state.n), t_end = $t_end
        model $(nameof(typeof(problem.model))), grid $(grid.Nx) × $(grid.Ny), $(typeof(sim.Δt)), \
        backend $(nameof(typeof(grid.backend))), $(Threads.nthreads()) threads
        $(join(("$k = $v" for (k, v) in pairs(run_parameters(sim))), ", "))
        """
    write_log!(sim, header)
    return nothing
end

"""
    run!(sim::Simulation; t_end) -> sim

Advance `sim` from its current time ``t_0`` to `t_end`, write the output classes of `sim.outputs` to
`sim.outdir` and record the runtime log `run.log` (ADR-0009, ADR-0010).

**Steps and times.** The time of step ``n`` is ``t_n = t_0 + (n - n_0)\\Delta t``, evaluated in Float64 from the
step count and assigned to `sim.state.t`, so that the time does not drift by accumulated rounding. The run
stops at the first step with ``t_n \\ge t_\\mathrm{end} - \\Delta t/2``, the step nearest to `t_end`.

**Sampling.** Each output class samples at ``t_0``, then at the step nearest to each multiple ``jP`` of its
period ``P``: step ``n`` samples when ``\\lfloor (t_n + \\Delta t/2)/P \\rfloor`` exceeds its value at step
``n - 1``. The multiples refer to ``t = 0``, so a restart keeps the sampling times. At ``t_0`` the injection
rates `inj_*` are NaN, because no step has run (ADR-0019). GLOBAL and SPECTRAL append to `global.nc` and
`spectral.nc`; FIELDS writes the snapshot `fields/snap_NNNNNN.nc`, numbered by the snapshot index, and
`params.toml` holds the [`parameters`](@ref) with the driver settings under `run`.

**Checks.** At ``t_0`` and every `n_check` steps, before the sampling, the run checks the state fields: every
value finite (not NaN or Inf), the density positive at every grid point (compressible model), and the CFL
number of [`Simulation`](@ref) at most `C_stop`. On a failure, `run!` writes a snapshot of the failed state
and throws an `ErrorException` that names the failure.

**Walltime.** After each step, and at ``t_0``, a wall-clock time since the call of at least `max_walltime`
writes a snapshot and ends the call normally. Neither stop writes a second snapshot of a step that FIELDS has
just sampled.

**Files of a run.** After [`set_initial_condition!`](@ref) (or at rest), `run!` starts a new run: it creates
the files, and it throws an `ArgumentError` when `outdir` holds the files of another run, unless
`overwrite = true`, which deletes them first. After [`restart!`](@ref) or a previous `run!`, the run continues:
`global.nc` and `spectral.nc` keep their records of the steps before ``n_0`` and drop the rest; `run.log`
gets a new segment.

**Runtime log.** `run.log` holds a header per call, one progress line per FIELDS sample, the stop reason,
and the TimerOutputs table, rewritten at each FIELDS sample and at the stop. The table times the steps, the
checks, and each output class, with one entry per diagnostic; the remainder of a class is its file output.
With `profile = true`, it also times the tasks of [`step!`](@ref).

The function throws an `ArgumentError` for `t_end` before the current time.
"""
function run!(sim::Simulation{T}; t_end::Real) where {T}
    start = time()
    (; state, problem, workspace) = sim
    t_end >= state.t || throw(ArgumentError("t_end = $t_end lies before the current time $(state.t)"))
    if sim.new_run
        files = run_files(sim)
        isempty(files) || sim.overwrite ||
            throw(ArgumentError("$(sim.outdir) holds the output of another run; choose another outdir or set overwrite = true"))
        foreach(rm, files)
    end
    reset_timer!(sim.timer)
    write_header!(sim, t_end)
    mkpath(fields_dir(sim))
    run = run_parameters(sim)
    write_parameters(joinpath(sim.outdir, "params.toml"), problem; run)
    empty!(workspace.diagnostics.injection)
    classes = keys(sim.diagnostics)
    datasets = NamedTuple{classes}(map(class -> open_class(sim, class, run), classes))
    sim.new_run = false
    sim.snapshot_step = -1
    segment = (; datasets, run, start,
               diagnostics = map(d -> [name => TimedDiagnostic(f, sim.timer, String(name)) for (name, f) in d],
                                 sim.diagnostics))
    t₀, n₀, Δt = Float64(state.t), state.n, Float64(sim.Δt)
    s_end = max(0, ceil(Int, (t_end - t₀) / Δt - 1 / 2))
    stop = "t_end = $t_end reached"
    try
        check!(sim, segment, 0)
        foreach(class -> sample!(sim, segment, class), keys(sim.outputs))
        for s in 1:s_end
            time() - start >= sim.max_walltime && break
            @timeit sim.timer "step" begin
                step!(state, problem, workspace, sim.Δt)
                KernelAbstractions.synchronize(sim.grid.backend)
            end
            t, t_previous = t₀ + s * Δt, t₀ + (s - 1) * Δt
            state.t = T(t)
            check!(sim, segment, s)
            for (class, output) in pairs(sim.outputs)
                period_index(t, output.period, Δt) > period_index(t_previous, output.period, Δt) &&
                    sample!(sim, segment, class)
            end
        end
        if state.n - n₀ < s_end
            stop = "max_walltime = $(sim.max_walltime) s reached"
            sim.snapshot_step == state.n || write_snapshot!(sim, segment)
            @info "run! stopped at the walltime limit at t = $(state.t), n = $(state.n)"
        end
    catch e
        stop = "error: $(sprint(showerror, e))"
        rethrow()
    finally
        foreach(close, datasets)
        write_log!(sim, "stop: $stop; " * progress_line(sim, start))
    end
    return sim
end

# Index of the period multiple jP that the step at time t has reached, with the step nearest to jP counted as
# reaching it.
period_index(t, period, Δt) = floor(Int, (t + Δt / 2) / period)

function write_snapshot!(sim::Simulation, segment)
    (; state, problem, workspace) = sim
    @timeit sim.timer "FIELDS" write_snapshot(fields_dir(sim), sim.next_snapshot, state, problem, workspace;
                                              segment.run)
    sim.next_snapshot += 1
    sim.snapshot_step = state.n
    return nothing
end

function sample!(sim::Simulation, segment, class)
    if class === :FIELDS
        write_snapshot!(sim, segment)
        write_log!(sim, progress_line(sim, segment.start))
    else
        write_class! = class === :GLOBAL ? write_global! : write_spectral!
        @timeit sim.timer String(class) write_class!(segment.datasets[class], sim.state, sim.problem, sim.workspace,
                                                     segment.diagnostics[class])
    end
    return nothing
end

# Check the state at the steps s that are multiples of n_check; on a failure, write a snapshot and throw.
function check!(sim::Simulation, segment, s)
    s % sim.n_check == 0 || return nothing
    failure = @timeit sim.timer "check" check_state(sim)
    isnothing(failure) && return nothing
    sim.snapshot_step == sim.state.n || write_snapshot!(sim, segment)
    error("$failure at t = $(sim.state.t), n = $(sim.state.n); snapshot $(sim.next_snapshot - 1) holds this state")
end
