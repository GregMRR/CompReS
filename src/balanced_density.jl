# Balanced density of a solenoidal velocity field: the initial ρ that makes θ = 0 and ∂_tθ = 0 at t = 0
# (WP-14, ADR-0013).

# ĥ ← −k_ik_j(u_iu_j)^/k² on the disc, zero at k = 0 and outside the disc; ĥ holds (u_xu_x)^ on entry.
@kernel function enthalpy_poisson_kernel!(ĥ, xy, yy, keep, kx, ky, k2)
    I = @index(Global)
    @inbounds begin
        a, b, q = kx[I], ky[I], k2[I]
        h = -(a * a * ĥ[I] + 2 * a * b * xy[I] + b * b * yy[I]) / ifelse(q == 0, one(q), q)
        ĥ[I] = ifelse(keep[I] & (q != 0), h, zero(h))
    end
end

# The constant h₀ with ⟨ρ(h + h₀)⟩ = 1 for the physical field h of zero mean, ρ = density_from_enthalpy.
# Newton steps on the increasing function f(h₀) = ⟨ρ(h + h₀)⟩ − 1, with f' = ⟨ρ/c_s²(ρ)⟩ because dh/dρ = c_s²/ρ.
# A step that leaves the bracket [h(1) − max h, h(1) − min h], where f changes sign, becomes a bisection step.
function enthalpy_constant(eos, h)
    T = eltype(h)
    mean_density(h₀) = sum(x -> density_from_enthalpy(eos, x + h₀), h) / length(h)
    slope(h₀) = sum(x -> (ρ = density_from_enthalpy(eos, x + h₀); ρ / sound_speed²(eos, ρ)), h) / length(h)
    h₀ = enthalpy(eos, one(T))
    lo, hi = h₀ - maximum(h), h₀ - minimum(h)
    # A change of h₀ below this bound changes ⟨ρ⟩ by less than eps(T).
    tolerance = eps(T) * sound_speed²(eos, one(T))
    for _ in 1:100
        f = mean_density(h₀) - 1
        f < 0 ? (lo = h₀) : (hi = h₀)
        next = h₀ - f / slope(h₀)
        lo < next < hi || (next = (lo + hi) / 2)
        converged = abs(next - h₀) <= tolerance
        h₀ = next
        converged && break
    end
    return h₀
end

"""
    balanced_density(u, eos, grid::Grid{T}) -> (; rho)

The density that balances the solenoidal velocity `u` at ``t = 0``, as a spectral field of float type `T`
(ADR-0013). `u` is a `NamedTuple` with the spectral fields `ux` and `uy`, for example the output of
[`random_velocity`](@ref) or [`kolmogorov`](@ref); `eos` is the equation of state of the run.
Merge the result with `u` to get the fields of a compressible [`State`](@ref):
`merge(u, balanced_density(u, eos, grid))`.

A uniform ``\\rho`` with a solenoidal ``\\mathbf{u}`` emits acoustic waves of amplitude ``O(Ma^2)``. For
``\\mu = \\zeta = \\kappa = 0`` and ``\\theta = \\nabla\\cdot\\mathbf{u} = 0``, the compressible momentum equation
(specification, compressible model) gives ``\\partial_t\\theta = -\\partial_i\\partial_j(u_iu_j) - \\nabla^2 h(\\rho)``.
The function makes this tendency zero in three steps:

1. Solve the Poisson equation for the enthalpy on the Galerkin truncation, with zero mean:

```math
\\nabla^2 h = -\\partial_i\\partial_j(u_iu_j), \\qquad
\\hat h(\\mathbf{k}) = -\\frac{k_ik_j\\,\\widehat{(u_iu_j)}(\\mathbf{k})}{k^2} \\quad (0 < |\\mathbf{k}| < k_c),
\\qquad \\hat h(\\mathbf{0}) = 0 .
```

2. Find the constant ``h_0`` by a scalar root find (Newton steps with a bisection safeguard) such that

```math
\\big\\langle \\rho(h + h_0) \\big\\rangle = 1, \\qquad \\rho(h) = \\texttt{density\\_from\\_enthalpy}(h) .
```

3. Return ``\\hat\\rho``, the Fourier coefficients of ``\\rho(h + h_0)`` on the grid, truncated to the disc.

The constant goes on ``h``, not on ``\\rho``: a shift of ``\\rho`` after the inversion changes ``\\nabla h`` for
a nonlinear ``h(\\rho)``. The truncation of ``\\hat\\rho`` keeps the mean ``\\hat\\rho(\\mathbf{0}) = 1`` and leaves
a residual of ``\\partial_t\\theta(0)`` that falls spectrally with ``N``. With viscosity,
``\\partial_t\\theta(0) = O(\\mu Ma^2)``, because ``\\partial_i(1/\\rho)\\nabla^2u_i \\neq 0``. The products
``u_iu_j`` are evaluated on the grid; for a velocity inside the truncation their aliases fall outside the disc,
so ``\\hat h`` is exact on the disc.

The function converts `eos` to `T` (ADR-0006) and allocates its own [`Transforms`](@ref) and work arrays: call it
once per run, before the time loop. It leaves `u` unchanged. The routine assumes a low Mach number: for
[`Polytropic`](@ref), ``h + h_0`` must stay positive.
"""
function balanced_density(u, eos, grid::Grid{T}) where {T}
    eos = convert_float(T, eos)
    tr = Transforms(grid)
    ux = to_physical!(KernelAbstractions.allocate(grid.backend, T, grid.Nx, grid.Ny), u.ux, tr)
    uy = to_physical!(similar(ux), u.uy, tr)
    spectral() = KernelAbstractions.allocate(grid.backend, Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny)
    ĥ, xy, yy = to_spectral!(spectral(), ux .* ux, tr), to_spectral!(spectral(), ux .* uy, tr),
                to_spectral!(spectral(), uy .* uy, tr)
    enthalpy_poisson_kernel!(grid.backend)(ĥ, xy, yy, grid.truncation, grid.kx, grid.ky, grid.k2;
                                           ndrange = size(ĥ))
    h = to_physical!(ux, ĥ, tr)
    h₀ = enthalpy_constant(eos, h)
    ρ = map!(x -> density_from_enthalpy(eos, x + h₀), uy, h)
    return (rho = truncate!(to_spectral!(ĥ, ρ, tr), grid),)
end
