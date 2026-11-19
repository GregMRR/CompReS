# Equations of state: the barotropic law `p(ρ)` and the functions derived from it (ADR-0008).

"""
    Isothermal(; cs)

Isothermal equation of state ``p = c_s^2 \\rho`` with constant sound speed `cs`
(specification, "Equation of state"). The enthalpy and the internal energy use the reference
density ``\\rho_0 = 1``.
"""
struct Isothermal{T<:AbstractFloat}
    cs::T
end

Isothermal(; cs::Real) = Isothermal(float(cs))

"""
    Polytropic(; K, γ)

Polytropic equation of state ``p = K \\rho^\\gamma`` with ``\\gamma \\neq 1``
(specification, "Equation of state").
"""
struct Polytropic{T<:AbstractFloat}
    K::T
    γ::T
end

function Polytropic(; K::Real, γ::Real)
    γ == 1 && throw(ArgumentError("Polytropic needs γ ≠ 1; use Isothermal for γ = 1"))
    T = float(promote_type(typeof(K), typeof(γ)))
    return Polytropic{T}(T(K), T(γ))
end

"""
    convert_float(T, x)

Return the parameter struct `x` with every floating-point parameter converted to the float type `T`
(ADR-0006). Each parameter struct of the package has one method.
"""
function convert_float end

convert_float(::Type{T}, eos::Isothermal) where {T<:AbstractFloat} = Isothermal{T}(T(eos.cs))
convert_float(::Type{T}, eos::Polytropic) where {T<:AbstractFloat} = Polytropic{T}(T(eos.K), T(eos.γ))

"""
    pressure(eos, ρ)

Pressure ``p(\\rho)`` of the equation of state `eos`: ``c_s^2 \\rho`` for [`Isothermal`](@ref),
``K \\rho^\\gamma`` for [`Polytropic`](@ref).
"""
pressure(eos::Isothermal, ρ) = eos.cs^2 * ρ
pressure(eos::Polytropic, ρ) = eos.K * ρ^eos.γ

"""
    enthalpy(eos, ρ)

Enthalpy ``h(\\rho)`` with ``h' = p'/\\rho`` (specification, "Equation of state"):
``c_s^2 \\ln\\rho`` for [`Isothermal`](@ref), ``K\\gamma\\rho^{\\gamma-1}/(\\gamma-1)`` for
[`Polytropic`](@ref).
"""
enthalpy(eos::Isothermal, ρ) = eos.cs^2 * log(ρ)
enthalpy(eos::Polytropic, ρ) = eos.K * eos.γ / (eos.γ - 1) * ρ^(eos.γ - 1)

"""
    internal_energy(eos, ρ)

Internal energy per unit mass ``e(\\rho)`` with ``e' = p/\\rho^2`` (specification, "Equation of
state"): ``c_s^2 \\ln\\rho`` for [`Isothermal`](@ref), ``K\\rho^{\\gamma-1}/(\\gamma-1)`` for
[`Polytropic`](@ref).
"""
internal_energy(eos::Isothermal, ρ) = eos.cs^2 * log(ρ)
internal_energy(eos::Polytropic, ρ) = eos.K / (eos.γ - 1) * ρ^(eos.γ - 1)

"""
    sound_speed²(eos, ρ)

Squared sound speed ``c_s^2 = p'(\\rho)`` (specification, "Equation of state"): ``c_s^2`` for
[`Isothermal`](@ref), ``K\\gamma\\rho^{\\gamma-1}`` for [`Polytropic`](@ref).
"""
sound_speed²(eos::Isothermal, ρ) = eos.cs^2 + zero(ρ)
sound_speed²(eos::Polytropic, ρ) = eos.K * eos.γ * ρ^(eos.γ - 1)

"""
    density_from_enthalpy(eos, h)

Density ``\\rho`` with ``h(\\rho) = h``, the inverse of [`enthalpy`](@ref):
``\\exp(h/c_s^2)`` for [`Isothermal`](@ref),
``\\big((\\gamma-1)h/(K\\gamma)\\big)^{1/(\\gamma-1)}`` for [`Polytropic`](@ref).
"""
density_from_enthalpy(eos::Isothermal, h) = exp(h / eos.cs^2)
density_from_enthalpy(eos::Polytropic, h) = ((eos.γ - 1) * h / (eos.K * eos.γ))^(1 / (eos.γ - 1))
