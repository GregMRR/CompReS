# Explicit Runge–Kutta schemes from Butcher tables, with and without the integrating factor (WP-7, ADR-0005, ADR-0023).

"""
    RungeKuttaScheme

Supertype of the explicit Runge–Kutta schemes [`Euler`](@ref), [`RK2`](@ref), [`SSPRK3`](@ref) and
[`RK4`](@ref). Each scheme is a Butcher table ``(a_{ij}, b_i, c_i)`` with ``s`` stages, read by
[`butcher_table`](@ref). For ``\\partial_t \\hat{\\mathbf{u}} = F(\\hat{\\mathbf{u}})`` one step reads

```math
\\hat{\\mathbf{u}}^{(i)} = \\hat{\\mathbf{u}}_n + \\Delta t \\sum_{j<i} a_{ij} \\mathbf{k}_j, \\qquad
\\mathbf{k}_i = F(\\hat{\\mathbf{u}}^{(i)}), \\qquad
\\hat{\\mathbf{u}}_{n+1} = \\hat{\\mathbf{u}}_n + \\Delta t \\sum_i b_i \\mathbf{k}_i .
```
"""
abstract type RungeKuttaScheme end

"""
    Euler()

The forward Euler scheme, order 1, one stage:

```math
\\begin{array}{c|c} 0 & \\\\ \\hline & 1 \\end{array}
```

It is unstable on undamped oscillating modes; use it for tests only (ADR-0005).
"""
struct Euler <: RungeKuttaScheme end

"""
    RK2()

Heun's scheme, order 2, two stages:

```math
\\begin{array}{c|cc} 0 & & \\\\ 1 & 1 & \\\\ \\hline & \\tfrac12 & \\tfrac12 \\end{array}
```

It is unstable on undamped oscillating modes; use it for tests only (ADR-0005).
"""
struct RK2 <: RungeKuttaScheme end

"""
    SSPRK3()

The strong-stability-preserving scheme of Shu and Osher, order 3, three stages:

```math
\\begin{array}{c|ccc} 0 & & & \\\\ 1 & 1 & & \\\\ \\tfrac12 & \\tfrac14 & \\tfrac14 & \\\\ \\hline
& \\tfrac16 & \\tfrac16 & \\tfrac23 \\end{array}
```
"""
struct SSPRK3 <: RungeKuttaScheme end

"""
    RK4()

The classical Runge–Kutta scheme, order 4, four stages:

```math
\\begin{array}{c|cccc} 0 & & & & \\\\ \\tfrac12 & \\tfrac12 & & & \\\\ \\tfrac12 & 0 & \\tfrac12 & & \\\\
1 & 0 & 0 & 1 & \\\\ \\hline & \\tfrac16 & \\tfrac13 & \\tfrac13 & \\tfrac16 \\end{array}
```
"""
struct RK4 <: RungeKuttaScheme end

"""
    butcher_table(scheme::RungeKuttaScheme) -> (a, b, c)

The Butcher table of `scheme` with exact `Rational{Int}` entries: the strictly lower-triangular
``s \\times s`` matrix `a`, the weights `b` and the nodes `c`, with ``c_1 = 0``.
"""
butcher_table(::Euler) = (a = zeros(Rational{Int}, 1, 1), b = [1 // 1], c = [0 // 1])
butcher_table(::RK2) = (a = Rational{Int}[0 0; 1 0], b = [1 // 2, 1 // 2], c = [0 // 1, 1 // 1])
butcher_table(::SSPRK3) =
    (a = Rational{Int}[0 0 0; 1 0 0; 1//4 1//4 0], b = [1 // 6, 1 // 6, 2 // 3], c = [0 // 1, 1 // 1, 1 // 2])
butcher_table(::RK4) =
    (a = Rational{Int}[0 0 0 0; 1//2 0 0 0; 0 1//2 0 0; 0 0 1 0], b = [1 // 6, 1 // 3, 1 // 3, 1 // 6],
     c = [0 // 1, 1 // 2, 1 // 2, 1 // 1])

# One operation of a stage combination, applied to an accumulator that starts as û_n. term = 0 propagates
# the accumulator over coeff·Δt; term = j ≥ 1 adds coeff·Δt·k_j.
struct CombineOp{T<:AbstractFloat}
    term::Int
    coeff::T
end

# Operations that build û_n + Δt Σ_j w_j k_j at node `target`, for the weights `w` of the tendencies k_j at
# the nodes `c`. Without the integrating factor they are plain additions. With it, the tendency k_j lives at
# time c_j, and the Lawson form propagates each term to the target: e^{L c Δt} û_n + Δt Σ_j w_j
# e^{L (c - c_j) Δt} k_j. The operations add the terms in increasing c_j and propagate the accumulator
# between two successive nodes (Horner scheme), so the propagation count equals the number of distinct nodes.
function combine_ops(::Type{T}, w, c, target, integrating_factor) where {T}
    terms = sort!([j for j in eachindex(w) if !iszero(w[j])]; by = j -> c[j])
    ops = CombineOp{T}[]
    node = 0 // 1
    for j in terms
        if integrating_factor && c[j] != node
            push!(ops, CombineOp(0, T(c[j] - node)))
            node = c[j]
        end
        push!(ops, CombineOp(j, T(w[j])))
    end
    integrating_factor && target != node && push!(ops, CombineOp(0, T(target - node)))
    return ops
end

"""
    TimeStepper(scheme::RungeKuttaScheme, integrating_factor::Bool, fields::NamedTuple; timer = nothing)

The time-stepping part of the [`Workspace`](@ref): the scheme, the integrating-factor switch and the stage
buffers for the spectral `fields` of a [`State`](@ref), allocated once. [`step!`](@ref) reads it as
`workspace.stepper`.

The fields are `scheme`, `integrating_factor`, `nodes` (the ``c_i`` in the float type `T` of `fields`),
`tendencies` (one set of spectral fields per stage, for ``\\mathbf{k}_i``), `stage` (a [`State`](@ref)
whose fields hold the stage value ``\\hat{\\mathbf{u}}^{(i)}``) and `timer`. The stage combinations of the
Butcher table are precomputed as lists of additions and propagations.

`timer` is `nothing` or a `TimerOutput` of TimerOutputs.jl. With a `TimerOutput`, [`step!`](@ref) times each
of its tasks into it, with a device synchronisation before and after each task (the profile mode of the
runtime log, ADR-0010). The synchronisations slow a GPU run. With `nothing`, `step!` times nothing.
"""
struct TimeStepper{S<:RungeKuttaScheme,T<:AbstractFloat,F<:NamedTuple,ST<:State,TM<:Union{Nothing,TimerOutput}}
    scheme::S
    integrating_factor::Bool
    nodes::Vector{T}
    stage_ops::Vector{Vector{CombineOp{T}}}
    final_ops::Vector{CombineOp{T}}
    tendencies::Vector{F}
    stage::ST
    timer::TM
end

function TimeStepper(scheme::RungeKuttaScheme, integrating_factor::Bool, fields::NamedTuple; timer = nothing)
    T = real(eltype(first(fields)))
    (; a, b, c) = butcher_table(scheme)
    s = length(b)
    stage_ops = [combine_ops(T, a[i, :], c, c[i], integrating_factor) for i in 1:s]
    final_ops = combine_ops(T, b, c, 1 // 1, integrating_factor)
    tendencies = [map(similar, fields) for _ in 1:s]
    return TimeStepper(scheme, integrating_factor, T.(c), stage_ops, final_ops, tendencies,
                       State(map(similar, fields)), timer)
end

# Evaluate `expr`. When `timer` is a TimerOutput (profile mode), time `expr` under `label`, with a device
# synchronisation before and after it, so that the time covers the kernels of `expr` only.
macro profiled(timer, backend, label, expr)
    return quote
        local to = $(esc(timer))
        if to === nothing
            $(esc(expr))
        else
            KernelAbstractions.synchronize($(esc(backend)))
            @timeit to $(esc(label)) begin
                local value = $(esc(expr))
                KernelAbstractions.synchronize($(esc(backend)))
                value
            end
        end
    end
end

@kernel function combine_kernel!(out, x, w, k)
    I = @index(Global)
    @inbounds out[I] = x[I] + w * k[I]
end

# out = x + w k for every field; out may be x.
function combine!(out, x, w, k)
    for (o, xf, kf) in zip(values(out), values(x), values(k))
        combine_kernel!(get_backend(o))(o, xf, w, kf; ndrange = length(o))
    end
    return out
end

# Apply the operations to the accumulator that starts as `src`, writing into `dst` (`dst` may be `src`). Every
# stage after the first and the final combination have at least one nonzero weight, so `ops` is not empty.
function apply_ops!(dst, src, ops, tendencies, rates, Δt)
    from = src
    for op in ops
        if op.term == 0
            propagate!(dst, from, rates, op.coeff * Δt)
        else
            combine!(dst, from, op.coeff * Δt, tendencies[op.term])
        end
        from = dst
    end
    return nothing
end

# Write the tendency k = F(û) of the stage state `s` into k: the right-hand side, plus the linear terms when
# the integrating factor does not apply them.
function tendency!(k, s, problem, ws, st, backend)
    @profiled st.timer backend "rhs!" rhs!(k, s, problem, ws)
    st.integrating_factor || @profiled st.timer backend "linear terms" linear_tendency!(k, s.fields, problem.rates)
    return k
end

"""
    step!(state, problem, workspace, Δt) -> state

Advance `state` by one time step `Δt` with the explicit Runge–Kutta scheme and the integrating-factor switch
of `workspace.stepper` (a [`TimeStepper`](@ref)), then add `Δt` to `state.t` and 1 to `state.n`.

The order of the calls in one step is:

1. [`update_forcing!`](@ref)`(problem.forcing, state, workspace, Δt)`, once: the forcing
   ``\\hat{\\mathbf{f}}_n`` stays constant over the stages (ADR-0005).
2. [`record_injection!`](@ref)`(problem, state, workspace, Δt)`, with ``\\hat{\\mathbf{u}}_n`` and
   ``\\hat{\\mathbf{f}}_n``.
3. The stages: [`rhs!`](@ref)`(k, stage, problem, workspace)` gives the nonlinear term and the forcing,
   ``N(\\hat{\\mathbf{u}})``.

Without the integrating factor, each stage adds the linear terms with [`linear_tendency!`](@ref), so
``F = L\\hat{\\mathbf{u}} + N`` and the step is the scheme of [`RungeKuttaScheme`](@ref). With the integrating
factor (Lawson form), [`propagate!`](@ref) integrates the linear part exactly, ``E(\\tau) = e^{L\\tau}``:

```math
\\hat{\\mathbf{u}}^{(i)} = E(c_i\\Delta t)\\,\\hat{\\mathbf{u}}_n
+ \\Delta t \\sum_{j<i} a_{ij}\\, E\\big((c_i - c_j)\\Delta t\\big)\\, N(\\hat{\\mathbf{u}}^{(j)}), \\qquad
\\hat{\\mathbf{u}}_{n+1} = E(\\Delta t)\\,\\hat{\\mathbf{u}}_n
+ \\Delta t \\sum_i b_i\\, E\\big((1 - c_i)\\Delta t\\big)\\, N(\\hat{\\mathbf{u}}^{(i)}) .
```

Each scheme keeps its order with the integrating factor. In [`SSPRK3`](@ref), ``c_3 < c_2``, so the third
stage applies the growth factor ``E(-\\Delta t/2)`` to ``N(\\hat{\\mathbf{u}}^{(2)})``; it overflows when
``\\sigma\\Delta t/2`` exceeds `log(floatmax(T))` (88 for Float32) on a mode with rate ``\\sigma``.

When `workspace.stepper.timer` is a `TimerOutput` (profile mode, ADR-0010), the function times the tasks
`update_forcing!`, `record_injection!`, `rhs!`, `linear terms` and `stage combination` into it.

`problem` provides the fields `forcing` and `rates` and a method of [`rhs!`](@ref); `rates` takes
[`linear_tendency!`](@ref) and [`propagate!`](@ref) methods. `workspace` provides the field `stepper`, and
whatever `rhs!`, `update_forcing!` and `record_injection!` read. A test problem that defines these
hooks runs through the same `step!` as the solver. The function allocates no array: it writes the
stages into the buffers of `workspace.stepper` and the result into `state.fields`.
"""
function step!(state::State{T}, problem, ws, Δt) where {T}
    st = ws.stepper
    h = T(Δt)
    backend = get_backend(first(state.fields))
    @profiled st.timer backend "update_forcing!" update_forcing!(problem.forcing, state, ws, h)
    @profiled st.timer backend "record_injection!" record_injection!(problem, state, ws, h)
    ks = st.tendencies
    tendency!(ks[1], state, problem, ws, st, backend)
    stage = st.stage
    stage.n = state.n
    for i in 2:length(ks)
        @profiled st.timer backend "stage combination" apply_ops!(stage.fields, state.fields, st.stage_ops[i], ks,
                                                                  problem.rates, h)
        stage.t = state.t + st.nodes[i] * h
        tendency!(ks[i], stage, problem, ws, st, backend)
    end
    @profiled st.timer backend "stage combination" apply_ops!(state.fields, state.fields, st.final_ops, ks,
                                                              problem.rates, h)
    state.t += h
    state.n += 1
    return state
end
