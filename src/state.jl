# State of a run and the preallocated work arrays of the time loop (WP-4, ADR-0008).

# Field names of the two models, in storage order (ADR-0008).
const MODEL_FIELDS = ((:ux, :uy), (:ux, :uy, :rho))

"""
    State(fields::NamedTuple; t = 0, n = 0)

The evolved fields of a run at time `t` and step `n`.

`fields` holds spectral fields in the real-FFT layout of the [`Grid`](@ref), with element type
`Complex{T}`: `(ux, uy)` for the incompressible model, or `(ux, uy, rho)` for the compressible
model, whose state is `(\\mathbf{u}, \\rho)` (specification, compressible model; ADR-0002). The
`State` keeps the given arrays, not copies. The time `t` takes the float type `T` of the fields.

The time loop changes `t`, `n` and the values in the arrays; the set of arrays stays fixed.
The constructor throws an `ArgumentError` for other field names, and for fields of unequal size or
element type.
"""
mutable struct State{T<:AbstractFloat,F<:NamedTuple}
    const fields::F
    t::T
    n::Int
end

function State(fields::NamedTuple; t = 0, n = 0)
    keys(fields) in MODEL_FIELDS ||
        throw(ArgumentError("state fields must be $(join(MODEL_FIELDS, " or ")), not $(keys(fields))"))
    first_field = first(fields)
    all(f -> size(f) == size(first_field) && eltype(f) == eltype(first_field), fields) ||
        throw(ArgumentError("state fields must share one size and one element type"))
    T = real(eltype(first_field))
    return State{T,typeof(fields)}(fields, T(t), Int(n))
end

"""
    Workspace(grid::Grid, state::State)

The [`Transforms`](@ref) of `grid` and the work arrays of the time loop, allocated once on the
backend of the grid. The time loop allocates nothing else.

The fields are:

- `transforms`: the [`Transforms`](@ref) of `grid`.
- `physical`: a `NamedTuple` of real `Nx × Ny` arrays, one per state field and with the same name,
  for the physical-space values of that field (`ux`, `uy` and, compressible model only, `rho`).
- `vorticity`: a real `Nx × Ny` array for the physical-space vorticity `\\omega` (WP-6).

A work package adds its buffers as one more field of this struct, built in this constructor.
A buffer per state field is a `NamedTuple` built with `map` over `state.fields`, as `physical` is;
it then follows the model without a branch. Code reads a buffer by name, for example
`ws.physical.ux`, and the typed fields keep every access type-stable.
"""
struct Workspace{TR<:Transforms,P<:NamedTuple,A<:AbstractMatrix}
    transforms::TR
    physical::P
    vorticity::A
end

function Workspace(grid::Grid{T}, state::State{T}) where {T}
    physical = map(_ -> KernelAbstractions.zeros(grid.backend, T, grid.Nx, grid.Ny), state.fields)
    vorticity = KernelAbstractions.zeros(grid.backend, T, grid.Nx, grid.Ny)
    return Workspace(Transforms(grid), physical, vorticity)
end
