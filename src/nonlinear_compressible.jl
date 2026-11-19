# Right-hand side of the compressible model: continuity flux, advection, enthalpy gradient and the explicit
# remainder of the small-scale dissipation (WP-10, ADR-0002, ADR-0004).

# From the spectral state, write ω̂ into w and V̂ = (μ∇²u + ζ∇θ)^ into (vx, vy).
@kernel function vorticity_viscous_kernel!(w, vx, vy, ux, uy, kx, ky, μ, ζ)
    c, r = @index(Global, NTuple)
    @inbounds begin
        x, y, a, b = ux[c, r], uy[c, r], kx[c], ky[r]
        w[c, r] = i_times(a, y) - i_times(b, x)
        d = a * x + b * y
        k2 = a^2 + b^2
        vx[c, r] = -μ * k2 * x - ζ * a * d
        vy[c, r] = -μ * k2 * y - ζ * b * d
    end
end

# Overwrite the physical buffers, point by point:
#   (vx, vy) ← −ω ẑ × u + (1/ρ − 1)(vx, vy),   w ← ½|u|² + h(ρ),   (ux, uy) ← ρu.
@kernel function compressible_product_kernel!(ux, uy, w, vx, vy, ρ, eos)
    I = @index(Global)
    @inbounds begin
        x, y, ω, s = ux[I], uy[I], w[I], ρ[I]
        a = inv(s) - one(s)
        vx[I] = ω * y + a * vx[I]
        vy[I] = -ω * x + a * vy[I]
        w[I] = (x * x + y * y) / 2 + enthalpy(eos, s)
        ux[I] = s * x
        uy[I] = s * y
    end
end

# dρ ← −ik·(mx, my) s, truncated; dρ holds the unnormalised forward FFT of ρu_x on entry.
@kernel function continuity_kernel!(dρ, my, keep, kx, ky, s)
    c, r = @index(Global, NTuple)
    @inbounds dρ[c, r] = ifelse(keep[c, r], -s * (i_times(kx[c], dρ[c, r]) + i_times(ky[r], my[c, r])),
                                zero(eltype(dρ)))
end

# (dx, dy) ← s((dx, dy) − ik B), truncated; (dx, dy) and B hold unnormalised forward FFTs on entry.
@kernel function velocity_kernel!(dx, dy, B, keep, kx, ky, s)
    c, r = @index(Global, NTuple)
    @inbounds begin
        b = B[c, r]
        dx[c, r] = ifelse(keep[c, r], s * (dx[c, r] - i_times(kx[c], b)), zero(b))
        dy[c, r] = ifelse(keep[c, r], s * (dy[c, r] - i_times(ky[r], b)), zero(b))
    end
end

"""
    rhs!(dû, state, problem::Problem{<:Grid,<:Compressible}, workspace) -> dû

Write the tendency of the compressible model without its linear terms into `dû` (specification,
compressible model, Fourier space):

```math
\\partial_t \\hat\\rho = -ik_j\\, \\widehat{\\rho u_j}, \\qquad
\\partial_t \\hat u_i = -\\widehat{(\\omega\\, \\hat{\\mathbf{z}} \\times \\mathbf{u})}_i
- ik_i\\, \\widehat{\\big(\\tfrac12|\\mathbf{u}|^2 + h(\\rho)\\big)}
+ \\widehat{\\Big[\\Big(\\frac1\\rho - 1\\Big)\\big(\\mu\\nabla^2 u_i + \\zeta\\,\\partial_i\\theta\\big)\\Big]}
+ \\hat f_i ,
```

restricted to the Galerkin truncation. The Lamb form
``(\\mathbf{u}\\cdot\\nabla)\\mathbf{u} = \\nabla(\\tfrac12|\\mathbf{u}|^2) + \\omega\\,\\hat{\\mathbf{z}}\\times\\mathbf{u}``
gives the advection (specification, notes for the numerical tests, note 2), and the enthalpy
``h(\\rho)`` of the equation of state gives the pressure term. The last term is the explicit remainder
of the small-scale dissipation ``\\rho^{-1}(\\mu\\nabla^2\\mathbf{u} + \\zeta\\nabla\\theta)`` after its
linear part with ``\\rho_0 = 1`` (ADR-0004). All products are evaluated in physical space. One
evaluation takes 11 FFTs: 6 inverse (``u_x``, ``u_y``, ``\\rho``, ``\\omega`` and the two components of
``\\mu\\nabla^2\\mathbf{u} + \\zeta\\nabla\\theta``) and 5 forward (the two components of
``-\\omega\\,\\hat{\\mathbf{z}}\\times\\mathbf{u}`` plus the remainder, ``\\tfrac12|\\mathbf{u}|^2 + h``,
and the two components of ``\\rho\\mathbf{u}``). The forcing term ``\\hat f_i`` comes from
[`add_forcing!`](@ref).

The linear terms of [`linear_rates`](@ref) (the small-scale dissipation with ``\\rho_0 = 1``, the
large-scale dissipation and the density diffusivity ``-\\kappa k^2\\hat\\rho``) are not part of `dû`.
The time stepping adds them with [`linear_tendency!`](@ref), or applies them with [`propagate!`](@ref)
when it uses the integrating factor.

`dû` holds spectral fields named `(:ux, :uy, :rho)`, of the size and element type of `state.fields`;
the function overwrites every mode. The tendency of the mass ``\\langle\\rho\\rangle`` (the mode
``\\mathbf{k} = 0`` of `dû.rho`) is exactly zero. With `nonlinear = false` in the model, `dû` holds the
forcing term only, so ``\\rho`` stays fixed for ``\\kappa = 0``. The function leaves `state` unchanged
and overwrites the buffers `physical`, `vorticity` and `compressible` of `workspace`. It allocates no
array.
"""
function rhs!(dû, state, problem::Problem{<:Grid,<:Compressible}, ws)
    if problem.model.nonlinear
        compressible_terms!(dû, state.fields, problem.grid, problem.model, ws)
    else
        foreach(f -> fill!(f, zero(eltype(f))), dû)
    end
    return add_forcing!(dû, problem.forcing, state, ws)
end

# Write the tendency without the linear terms and the forcing into dû. Before the forward FFTs,
# dû holds ω̂ and the spectral μ∇²u + ζ∇θ.
function compressible_terms!(dû, û, grid::Grid{T}, model, ws) where {T}
    tr = ws.transforms
    backend = grid.backend
    (; viscous_x, viscous_y, spectral) = ws.compressible
    ux, uy, ρ = ws.physical.ux, ws.physical.uy, ws.physical.rho
    vorticity_viscous_kernel!(backend)(dû.rho, dû.ux, dû.uy, û.ux, û.uy, tr.kx_deriv, tr.ky_deriv,
                                       model.μ, model.ζ; ndrange = size(dû.ux))
    to_physical!(ws.vorticity, dû.rho, tr)
    to_physical!(viscous_x, dû.ux, tr)
    to_physical!(viscous_y, dû.uy, tr)
    to_physical!(ux, û.ux, tr)
    to_physical!(uy, û.uy, tr)
    to_physical!(ρ, û.rho, tr)
    compressible_product_kernel!(backend)(ux, uy, ws.vorticity, viscous_x, viscous_y, ρ, model.eos;
                                          ndrange = length(ux))
    s = inv(T(grid.Nx * grid.Ny))
    mul!(dû.rho, tr.forward, ux)
    mul!(spectral, tr.forward, uy)
    continuity_kernel!(backend)(dû.rho, spectral, grid.truncation, tr.kx_deriv, tr.ky_deriv, s;
                                ndrange = size(dû.rho))
    mul!(dû.ux, tr.forward, viscous_x)
    mul!(dû.uy, tr.forward, viscous_y)
    mul!(spectral, tr.forward, ws.vorticity)
    velocity_kernel!(backend)(dû.ux, dû.uy, spectral, grid.truncation, tr.kx_deriv, tr.ky_deriv, s;
                              ndrange = size(dû.ux))
    return dû
end
