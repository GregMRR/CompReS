# Model types, linear rates, explicit linear tendency and integrating-factor propagator (WP-5, ADR-0004, ADR-0005).

"""
    LargeScaleDissipation(; α = 0, m = 0)

Large-scale dissipation ``-\\alpha(-\\nabla^2)^{-m}\\mathbf{u}``, with Fourier rate
``\\alpha k^{-2m}`` on every mode except ``\\mathbf{k} = 0`` (ADR-0004). `m = 0` is linear friction;
`m ≥ 1` is hypoviscosity. The default `α = 0` switches the term off.
The constructor throws an `ArgumentError` for `α < 0` or `m < 0`.
"""
struct LargeScaleDissipation{T<:AbstractFloat}
    α::T
    m::Int
end

function LargeScaleDissipation(; α::Real = 0, m::Integer = 0)
    α >= 0 || throw(ArgumentError("α must be nonnegative, got $α"))
    m >= 0 || throw(ArgumentError("m must be nonnegative, got $m"))
    return LargeScaleDissipation(float(α), Int(m))
end

"""
    Incompressible(; ν, large_scale = LargeScaleDissipation(), nonlinear = true)

Incompressible model with kinematic viscosity `ν` and the [`LargeScaleDissipation`](@ref)
`large_scale`. Its linear part acts on each velocity mode as (specification, incompressible model,
Fourier form)

```math
\\partial_t \\hat u_i = -\\big(\\nu k^2 + \\alpha k^{-2m}\\big)\\hat u_i .
```

`nonlinear = false` is a test switch: the right-hand side then drops the advection and keeps
the linear terms and the forcing. The constructor throws an `ArgumentError` for `ν < 0`.
"""
struct Incompressible{T<:AbstractFloat}
    ν::T
    large_scale::LargeScaleDissipation{T}
    nonlinear::Bool
end

function Incompressible(; ν::Real, large_scale::LargeScaleDissipation = LargeScaleDissipation(),
                        nonlinear::Bool = true)
    ν >= 0 || throw(ArgumentError("ν must be nonnegative, got $ν"))
    T = float(promote_type(typeof(ν), typeof(large_scale.α)))
    return Incompressible{T}(T(ν), convert_float(T, large_scale), nonlinear)
end

"""
    Compressible(; eos, μ, ζ, κ = 0, large_scale = LargeScaleDissipation(), nonlinear = true)

Compressible model with the equation of state `eos`, shear viscosity `μ`, bulk viscosity `ζ`,
density diffusivity `κ` and the [`LargeScaleDissipation`](@ref) `large_scale`.

The small-scale dissipation ``\\rho^{-1}(\\mu\\nabla^2\\mathbf{u} + \\zeta\\nabla\\theta)`` splits into
a linear part with ``\\rho_0 = 1`` and an explicit remainder (ADR-0004). The linear part is diagonal in
the Helmholtz basis (specification, compressible model, Fourier form):

```math
\\partial_t \\hat{\\mathbf{u}}^s = -\\big(\\mu k^2 + \\alpha k^{-2m}\\big)\\hat{\\mathbf{u}}^s, \\qquad
\\partial_t \\hat{\\mathbf{u}}^c = -\\big((\\mu + \\zeta) k^2 + \\alpha k^{-2m}\\big)\\hat{\\mathbf{u}}^c, \\qquad
\\partial_t \\hat\\rho = -\\kappa k^2 \\hat\\rho .
```

`nonlinear = false` is a test switch: the right-hand side then drops the advection, the pressure,
the viscous remainder and the continuity flux, so `ρ` stays fixed for `κ = 0`. The constructor
throws an `ArgumentError` for a negative `μ`, `ζ` or `κ`.
"""
struct Compressible{T<:AbstractFloat,E}
    eos::E
    μ::T
    ζ::T
    κ::T
    large_scale::LargeScaleDissipation{T}
    nonlinear::Bool
end

function Compressible(; eos, μ::Real, ζ::Real, κ::Real = 0,
                      large_scale::LargeScaleDissipation = LargeScaleDissipation(), nonlinear::Bool = true)
    for (name, value) in ((:μ, μ), (:ζ, ζ), (:κ, κ))
        value >= 0 || throw(ArgumentError("$name must be nonnegative, got $value"))
    end
    T = float(promote_type(typeof(μ), typeof(ζ), typeof(κ), typeof(large_scale.α)))
    return Compressible{T,typeof(eos)}(eos, T(μ), T(ζ), T(κ), convert_float(T, large_scale), nonlinear)
end

convert_float(::Type{T}, ls::LargeScaleDissipation) where {T<:AbstractFloat} = LargeScaleDissipation{T}(T(ls.α), ls.m)
convert_float(::Type{T}, m::Incompressible) where {T<:AbstractFloat} =
    Incompressible{T}(T(m.ν), convert_float(T, m.large_scale), m.nonlinear)
function convert_float(::Type{T}, m::Compressible) where {T<:AbstractFloat}
    eos = convert_float(T, m.eos)
    return Compressible{T,typeof(eos)}(eos, T(m.μ), T(m.ζ), T(m.κ), convert_float(T, m.large_scale), m.nonlinear)
end

"""
    IncompressibleRates

Linear rates of the [`Incompressible`](@ref) model, from [`linear_rates`](@ref).
`velocity` holds ``\\sigma(\\mathbf{k}) = \\nu k^2 + \\alpha k^{-2m}`` (with ``\\sigma(0) = 0``) on the
real-FFT layout of the grid; both velocity components decay at this rate.
"""
struct IncompressibleRates{A<:AbstractMatrix}
    velocity::A
end

"""
    CompressibleRates

Linear rates of the [`Compressible`](@ref) model, from [`linear_rates`](@ref), on the real-FFT layout
of the grid (``\\sigma(0) = 0`` in each array):

- `solenoidal`: ``\\mu k^2 + \\alpha k^{-2m}`` for the solenoidal part of the velocity.
- `compressive`: ``(\\mu + \\zeta) k^2 + \\alpha k^{-2m}`` for the compressive part of the velocity.
- `density`: ``\\kappa k^2`` for the density.
- `kx`, `ky`: the grid wavenumbers, for the Helmholtz split inside the kernels. On the retained modes
  this split equals [`helmholtz_split!`](@ref); the two differ only on the Nyquist column and row,
  where [`helmholtz_split!`](@ref) uses zero derivative wavenumbers.
"""
struct CompressibleRates{A<:AbstractMatrix}
    solenoidal::A
    compressive::A
    density::A
    kx::A
    ky::A
end

# Large-scale rate α k^{-2m} for the squared wavenumber k2; zero at k = 0 (ADR-0004).
large_scale_rate(ls, k2) = iszero(k2) ? zero(k2) : ls.α / k2^ls.m

# Device array of the rate function `σ(k2)` on the grid.
function rate_array(σ, grid::Grid)
    host = map(σ, Array(grid.k2))
    return copyto!(KernelAbstractions.allocate(grid.backend, eltype(host), size(host)), host)
end

"""
    linear_rates(grid::Grid{T}, model) -> IncompressibleRates or CompressibleRates

Build the linear rates of `model` on `grid`, once per run: one device array of float type `T` per
field and Helmholtz component (see [`IncompressibleRates`](@ref) and [`CompressibleRates`](@ref)).
The function converts the model parameters to `T`. [`linear_tendency!`](@ref) and
[`propagate!`](@ref) read the arrays.
"""
function linear_rates(grid::Grid{T}, model::Incompressible) where {T}
    m = convert_float(T, model)
    return IncompressibleRates(rate_array(k2 -> m.ν * k2 + large_scale_rate(m.large_scale, k2), grid))
end

function linear_rates(grid::Grid{T}, model::Compressible) where {T}
    m = convert_float(T, model)
    ls = m.large_scale
    return CompressibleRates(rate_array(k2 -> m.μ * k2 + large_scale_rate(ls, k2), grid),
                             rate_array(k2 -> (m.μ + m.ζ) * k2 + large_scale_rate(ls, k2), grid),
                             rate_array(k2 -> m.κ * k2, grid), grid.kx, grid.ky)
end

# Map fs·u^s + fc·u^c of the velocity mode (x, y), diagonal in the Helmholtz basis. The propagator
# passes decay factors as fs and fc; the tendency passes negated rates.
@inline function helmholtz_diagonal(fs, fc, kx, ky, x, y)
    cx, cy = compressive_part(kx, ky, x, y)
    return fs * x + (fc - fs) * cx, fs * y + (fc - fs) * cy
end

@kernel function scalar_propagate_kernel!(ox, oy, ux, uy, σ, τ)
    I = @index(Global)
    @inbounds begin
        e = exp(-σ[I] * τ)
        ox[I] = e * ux[I]
        oy[I] = e * uy[I]
    end
end

@kernel function scalar_tendency_kernel!(dx, dy, ux, uy, σ)
    I = @index(Global)
    @inbounds begin
        dx[I] -= σ[I] * ux[I]
        dy[I] -= σ[I] * uy[I]
    end
end

@kernel function compressible_propagate_kernel!(ox, oy, oρ, ux, uy, ρ, σs, σc, σρ, kx, ky, τ)
    I = @index(Global)
    @inbounds begin
        ox[I], oy[I] = helmholtz_diagonal(exp(-σs[I] * τ), exp(-σc[I] * τ), kx[I], ky[I], ux[I], uy[I])
        oρ[I] = exp(-σρ[I] * τ) * ρ[I]
    end
end

@kernel function compressible_tendency_kernel!(dx, dy, dρ, ux, uy, ρ, σs, σc, σρ, kx, ky)
    I = @index(Global)
    @inbounds begin
        tx, ty = helmholtz_diagonal(-σs[I], -σc[I], kx[I], ky[I], ux[I], uy[I])
        dx[I] += tx
        dy[I] += ty
        dρ[I] -= σρ[I] * ρ[I]
    end
end

"""
    linear_tendency!(dû, û, rates) -> dû

Add the linear terms ``L\\hat{\\mathbf{u}}`` of the model to the tendency `dû`, for explicit time
stepping. `û` and `dû` are spectral fields named `(:ux, :uy)` for [`IncompressibleRates`](@ref) and
`(:ux, :uy, :rho)` for [`CompressibleRates`](@ref). Per mode,
``L\\hat u_i = -\\sigma\\hat u_i`` in the incompressible model, and

```math
L\\hat{\\mathbf{u}} = -\\sigma_s\\hat{\\mathbf{u}}^s - \\sigma_c\\hat{\\mathbf{u}}^c, \\qquad
L\\hat\\rho = -\\sigma_\\rho\\hat\\rho
```

in the compressible model, with the rates of [`linear_rates`](@ref). The function allocates nothing.
"""
function linear_tendency!(dû, û, r::IncompressibleRates)
    backend = get_backend(r.velocity)
    scalar_tendency_kernel!(backend)(dû.ux, dû.uy, û.ux, û.uy, r.velocity; ndrange = length(r.velocity))
    return dû
end

function linear_tendency!(dû, û, r::CompressibleRates)
    backend = get_backend(r.solenoidal)
    compressible_tendency_kernel!(backend)(dû.ux, dû.uy, dû.rho, û.ux, û.uy, û.rho, r.solenoidal,
                                           r.compressive, r.density, r.kx, r.ky; ndrange = length(r.solenoidal))
    return dû
end

"""
    propagate!(out, û, rates, τ) -> out

Write ``\\exp(L\\tau)\\hat{\\mathbf{u}}``, the integrating-factor propagator of the linear terms
(ADR-0005), into `out`. The fields are named as in [`linear_tendency!`](@ref); `out` may be `û`.
Per mode, the propagator is ``e^{-\\sigma\\tau}`` in the incompressible model, and

```math
e^{L\\tau}\\hat{\\mathbf{u}} = e^{-\\sigma_s\\tau} P\\hat{\\mathbf{u}}
+ e^{-\\sigma_c\\tau}(I - P)\\hat{\\mathbf{u}}, \\qquad
e^{L\\tau}\\hat\\rho = e^{-\\sigma_\\rho\\tau}\\hat\\rho
```

in the compressible model, with the solenoidal projection ``P_{ij} = \\delta_{ij} - k_ik_j/k^2``.
The kernel computes the exponentials from the rates, so any `τ` (for example ``c_i\\Delta t`` of a
Butcher table) costs no storage and no allocation. `τ` is converted to the float type of the rates.
"""
function propagate!(out, û, r::IncompressibleRates, τ)
    σ = r.velocity
    scalar_propagate_kernel!(get_backend(σ))(out.ux, out.uy, û.ux, û.uy, σ, eltype(σ)(τ); ndrange = length(σ))
    return out
end

function propagate!(out, û, r::CompressibleRates, τ)
    σ = r.solenoidal
    compressible_propagate_kernel!(get_backend(σ))(out.ux, out.uy, out.rho, û.ux, û.uy, û.rho, σ,
                                                   r.compressive, r.density, r.kx, r.ky, eltype(σ)(τ);
                                                   ndrange = length(σ))
    return out
end
