# Registry of the named diagnostics of the output classes GLOBAL and SPECTRAL (WP-8, ADR-0009, ADR-0022).

"""
    Diagnostic

A registered diagnostic: the function `f(state, problem, workspace)` and the `coordinate` of its values.
A GLOBAL diagnostic returns one number of the float type `T` of the state and has the coordinate
`:none`. A SPECTRAL diagnostic returns a `Vector{T}` with one value per shell of the grid, and its
coordinate names the [`Grid`](@ref) field of the wavenumbers where the values are reported:
`:shell_k` (shell centres) or `:shell_edge` (upper shell edges). Calling a `Diagnostic` with
`(state, problem, workspace)` calls `f`.

A diagnostic is a generic function: one registration serves every model, and each model adds its
method (ADR-0022).
"""
struct Diagnostic
    f::Function
    coordinate::Symbol
end

(d::Diagnostic)(state, problem, workspace) = d.f(state, problem, workspace)

# One registry per output class; FIELDS holds snapshots, not diagnostics.
const DIAGNOSTICS = (GLOBAL = Dict{Symbol,Diagnostic}(), SPECTRAL = Dict{Symbol,Diagnostic}())

function class_registry(class::Symbol)
    haskey(DIAGNOSTICS, class) ||
        throw(ArgumentError("output class $class has no diagnostics; use $(join(keys(DIAGNOSTICS), " or "))"))
    return DIAGNOSTICS[class]
end

"""
    register_diagnostic!(class::Symbol, name::Symbol, f; coordinate = :shell_k) -> Diagnostic

Register the function `f(state, problem, workspace)` as the diagnostic `name` of the output class
`class`, `:GLOBAL` or `:SPECTRAL` (ADR-0009). `coordinate` applies to `:SPECTRAL` only: `:shell_k` for
values at the shell centres ``k_n = n\\Delta k``, `:shell_edge` for values at the upper shell edges
``K_n = (n + \\tfrac12)\\Delta k``. A new registration under an existing name replaces the old one.

The function throws an `ArgumentError` for another class or coordinate.
"""
function register_diagnostic!(class::Symbol, name::Symbol, f; coordinate::Symbol = :shell_k)
    registry = class_registry(class)
    class === :SPECTRAL && !(coordinate in (:shell_k, :shell_edge)) &&
        throw(ArgumentError("a SPECTRAL coordinate is :shell_k or :shell_edge, not :$coordinate"))
    return registry[name] = Diagnostic(f, class === :GLOBAL ? :none : coordinate)
end

"""
    select_diagnostics(class::Symbol, names) -> Vector{Pair{Symbol,Diagnostic}}

The registered diagnostics `names` of the output class `class`, as `name => diagnostic` pairs in the
order of `names`. The function throws an `ArgumentError` for an unknown class or name.
"""
function select_diagnostics(class::Symbol, names)
    registry = class_registry(class)
    for name in names
        haskey(registry, name) ||
            throw(ArgumentError("no $class diagnostic :$name; registered: $(join(sort!(collect(keys(registry))), ", "))"))
    end
    return [name => registry[name] for name in names]
end

"""
    default_diagnostics(class::Symbol, model) -> Tuple of Symbol

The names of the default diagnostics of the output class `class` for `model` (ADR-0009).
"""
function default_diagnostics end

"""
    record_injection!(problem, state, workspace, Δt) -> workspace

Store the injection rates of the step `Δt` that starts at `state` in `workspace.diagnostics.injection`,
under the names of their GLOBAL diagnostics. The time stepping calls it once per step, right after
[`update_forcing!`](@ref), so `state` holds ``\\hat{\\mathbf{u}}_n`` and `workspace.forcing` holds
``\\hat{\\mathbf{f}}_n`` (ADR-0016, ADR-0017). A diagnostic sees the state after the step, so the
`inj_*` diagnostics read the stored values; a rate not yet stored reads as `NaN`.

This method, for any problem without its own method, stores nothing.
"""
record_injection!(problem, state, workspace, Δt) = workspace

# GLOBAL diagnostic function that reads the injection rate `name` stored by record_injection!.
stored_injection(name::Symbol) =
    (state, problem, workspace) -> (rates = workspace.diagnostics.injection; get(rates, name, valtype(rates)(NaN)))

# Shell and mode reductions shared by the diagnostics of every model. A per-mode function
# f(kx, ky, ax, ay, bx, by) receives the derivative wavenumbers of the mode and the values of two spectral
# vector fields a and b there; the reductions weight it by the half-plane weight w (1 on the kx = 0
# column and the Nyquist column, 2 elsewhere), so a sum over the stored modes is a sum over the full plane.

@kernel function mode_density_kernel!(density, f, a, b, kx, ky, nkx)
    c, r = @index(Global, NTuple)
    @inbounds density[c, r] = ifelse(c == 1 || c == nkx, 1, 2) *
                              f(kx[c], ky[r], a.ux[c, r], a.uy[c, r], b.ux[c, r], b.uy[c, r])
end

@kernel function shell_sum_kernel!(shells, density, order, start)
    n = @index(Global)
    s = zero(eltype(shells))
    @inbounds for p in start[n]:(start[n + 1] - 1)
        s += density[order[p]]
    end
    @inbounds shells[n] = s
end

# Write w·f of every mode into workspace.diagnostics.density.
function mode_density!(f, a, b, ws)
    tr = ws.transforms
    density = ws.diagnostics.density
    mode_density_kernel!(tr.grid.backend)(density, f, a, b, tr.kx_deriv, tr.ky_deriv, size(density, 1);
                                          ndrange = size(density))
    return density
end

# Σ w f over all modes.
mode_sum(f, a, b, ws) = sum(mode_density!(f, a, b, ws))

# Σ w f over the modes of each shell, as a host vector with one entry per shell of the grid.
function shell_sums(f, a, b, ws)
    d = ws.diagnostics
    mode_density!(f, a, b, ws)
    shell_sum_kernel!(ws.transforms.grid.backend)(d.shells, d.density, d.shell_order, d.shell_start;
                                                  ndrange = length(d.shells))
    return Array(d.shells)
end
