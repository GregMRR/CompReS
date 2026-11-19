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
- `compressible`: `nothing` for the incompressible model. For the compressible model, the buffers of
  its [`rhs!`](@ref) (WP-10): `viscous_x`, `viscous_y`, real `Nx × Ny` arrays for the physical-space
  `\\mu\\nabla^2\\mathbf{u} + \\zeta\\nabla\\theta`, and `spectral`, one spectral field.
- `forcing`: a `NamedTuple` `(ux, uy)` of spectral fields for the forcing `\\hat{\\mathbf{f}}` of the
  current step, which [`update_forcing!`](@ref) writes (WP-9). They are zero outside the forcing support.
- `diagnostics`: the buffers of the diagnostics (WP-8): `density`, a real array of the spectral size for
  one value per mode; `shells`, a real vector with one value per shell of the grid; `shell_order`
  and `shell_start`, the integer vectors of the shell reduction (the linear indices of all modes sorted
  by shell, and the position of the first mode of each shell in that order, plus one end position);
  `nonlinear`, a `NamedTuple` `(ux, uy)` of spectral fields for the nonlinear term of the fluxes; and
  `injection`, a `Dict{Symbol,T}` of the injection rates that [`record_injection!`](@ref) stores.

A work package adds its buffers as one more field of this struct, built in this constructor.
A buffer per state field is a `NamedTuple` built with `map` over `state.fields`, as `physical` is;
it then follows the model without a branch. Code reads a buffer by name, for example
`ws.physical.ux`, and the typed fields keep every access type-stable.
"""
struct Workspace{TR<:Transforms,P<:NamedTuple,A<:AbstractMatrix,C,F<:NamedTuple,D<:NamedTuple}
    transforms::TR
    physical::P
    vorticity::A
    compressible::C
    forcing::F
    diagnostics::D
end

# Buffers of the compressible rhs! (WP-10); the incompressible model needs none.
compressible_buffers(grid, ::NamedTuple{(:ux, :uy)}) = nothing
compressible_buffers(grid::Grid{T}, ::NamedTuple{(:ux, :uy, :rho)}) where {T} =
    (viscous_x = KernelAbstractions.zeros(grid.backend, T, grid.Nx, grid.Ny),
     viscous_y = KernelAbstractions.zeros(grid.backend, T, grid.Nx, grid.Ny),
     spectral = KernelAbstractions.zeros(grid.backend, Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny))

# Buffers of the diagnostics (WP-8). The shell reduction reads the modes of shell n at the positions
# shell_start[n]:(shell_start[n + 1] - 1) of shell_order.
function diagnostic_buffers(grid::Grid{T}) where {T}
    shell = vec(Array(grid.shell_index))
    nshell = length(grid.shell_k)
    order = Int32.(sortperm(shell))
    start = Int32.(searchsortedfirst.(Ref(shell[order]), 1:(nshell + 1)))
    on_device(a) = copyto!(KernelAbstractions.allocate(grid.backend, eltype(a), size(a)), a)
    spectral() = KernelAbstractions.zeros(grid.backend, Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny)
    return (density = KernelAbstractions.zeros(grid.backend, T, grid.Nx ÷ 2 + 1, grid.Ny),
            shells = KernelAbstractions.zeros(grid.backend, T, nshell),
            shell_order = on_device(order), shell_start = on_device(start),
            nonlinear = (ux = spectral(), uy = spectral()), injection = Dict{Symbol,T}())
end

function Workspace(grid::Grid{T}, state::State{T}) where {T}
    physical = map(_ -> KernelAbstractions.zeros(grid.backend, T, grid.Nx, grid.Ny), state.fields)
    vorticity = KernelAbstractions.zeros(grid.backend, T, grid.Nx, grid.Ny)
    forcing = (ux = KernelAbstractions.zeros(grid.backend, Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny),
               uy = KernelAbstractions.zeros(grid.backend, Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny))
    return Workspace(Transforms(grid), physical, vorticity, compressible_buffers(grid, state.fields), forcing,
                     diagnostic_buffers(grid))
end
