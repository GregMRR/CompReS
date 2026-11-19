# Problem, forcing hook and right-hand side of the incompressible model in Lamb form (WP-6, ADR-0002, ADR-0003).

"""
    NoForcing()

The absence of forcing: [`add_forcing!`](@ref) adds nothing to the tendency.
"""
struct NoForcing end

"""
    add_forcing!(dû, forcing, state, workspace) -> dû

Add the spectral forcing ``\\hat{\\mathbf{f}}`` of `forcing` to the tendency `dû`. [`rhs!`](@ref) calls
this function after the nonlinear term. The method for [`NoForcing`](@ref) leaves `dû` unchanged.
"""
add_forcing!(dû, ::NoForcing, state, workspace) = dû

"""
    Problem(grid::Grid{T}, model, forcing)

Everything the right-hand side needs: `grid`, `model` (an [`Incompressible`](@ref) or
[`Compressible`](@ref) model), `forcing` and the linear rates `rates` of [`linear_rates`](@ref).
The constructor converts `model` to the float type `T` of the grid (ADR-0006).
[`rhs!`](@ref) evaluates the tendency of a `Problem`; the time stepping advances it.
"""
struct Problem{G<:Grid,M,F,R}
    grid::G
    model::M
    forcing::F
    rates::R
end

function Problem(grid::Grid{T}, model, forcing) where {T}
    m = convert_float(T, model)
    return Problem(grid, m, forcing, linear_rates(grid, m))
end

@kernel function vorticity_kernel!(ω̂, ux, uy, kx, ky)
    c, r = @index(Global, NTuple)
    @inbounds ω̂[c, r] = i_times(kx[c], uy[c, r]) - i_times(ky[r], ux[c, r])
end

# Overwrite (ux, uy) with (ω uy, −ω ux) = −ω ẑ × u, point by point.
@kernel function lamb_product_kernel!(ux, uy, ω)
    I = @index(Global)
    @inbounds begin
        w, x, y = ω[I], ux[I], uy[I]
        ux[I] = w * y
        uy[I] = -w * x
    end
end

# Scale the unnormalised forward FFTs (nx, ny) by s, project them on the solenoidal part and set the
# modes outside the Galerkin truncation to zero, in one pass.
@kernel function project_truncate_kernel!(nx, ny, keep, kx, ky, s)
    c, r = @index(Global, NTuple)
    @inbounds begin
        x, y = s * nx[c, r], s * ny[c, r]
        cx, cy = compressive_part(kx[c], ky[r], x, y)
        nx[c, r] = ifelse(keep[c, r], x - cx, zero(x))
        ny[c, r] = ifelse(keep[c, r], y - cy, zero(y))
    end
end

"""
    rhs!(dû, state, problem::Problem{<:Grid,<:Incompressible}, workspace) -> dû

Write the tendency of the incompressible model without its linear terms into `dû`:

```math
\\partial_t \\hat u_i = -P_{il}(\\mathbf{k})\\, \\widehat{(\\omega\\, \\hat{\\mathbf{z}} \\times \\mathbf{u})}_l
+ \\hat f_i ,
```

restricted to the Galerkin truncation (specification, incompressible model, Fourier space). The
Lamb form ``(\\mathbf{u}\\cdot\\nabla)\\mathbf{u} = \\nabla(\\tfrac12|\\mathbf{u}|^2) + \\omega\\,
\\hat{\\mathbf{z}} \\times \\mathbf{u}`` gives the advection: the projector ``P_{il}`` removes the
gradient. One evaluation takes 5 FFTs: 3 inverse (``u_x``, ``u_y``, ``\\omega``) and 2 forward.
The forcing term ``\\hat f_i`` comes from [`add_forcing!`](@ref).

The linear terms ``-\\sigma \\hat u_i`` of [`linear_rates`](@ref) are not part of `dû`. The time
stepping adds them with [`linear_tendency!`](@ref), or applies them with [`propagate!`](@ref) when it
uses the integrating factor.

`dû` holds spectral fields named `(:ux, :uy)`, of the size and element type of `state.fields`; the
function overwrites every mode. With `nonlinear = false` in the model, the advection is zero and
`dû` holds the forcing term only. The function leaves `state` unchanged and overwrites the buffers
`physical` and `vorticity` of `workspace`. It allocates no array.
"""
function rhs!(dû, state, problem::Problem{<:Grid,<:Incompressible}, ws)
    if problem.model.nonlinear
        advection!(dû, state.fields, problem.grid, ws)
    else
        fill!(dû.ux, zero(eltype(dû.ux)))
        fill!(dû.uy, zero(eltype(dû.uy)))
    end
    return add_forcing!(dû, problem.forcing, state, ws)
end

# Write −P(ω ẑ × u)^, truncated, into dû. dû.ux holds ω̂ until the forward FFTs overwrite it.
function advection!(dû, û, grid::Grid{T}, ws) where {T}
    tr = ws.transforms
    backend = grid.backend
    ux, uy = ws.physical.ux, ws.physical.uy
    vorticity_kernel!(backend)(dû.ux, û.ux, û.uy, tr.kx_deriv, tr.ky_deriv; ndrange = size(dû.ux))
    to_physical!(ws.vorticity, dû.ux, tr)
    to_physical!(ux, û.ux, tr)
    to_physical!(uy, û.uy, tr)
    lamb_product_kernel!(backend)(ux, uy, ws.vorticity; ndrange = length(ux))
    mul!(dû.ux, tr.forward, ux)
    mul!(dû.uy, tr.forward, uy)
    project_truncate_kernel!(backend)(dû.ux, dû.uy, grid.truncation, tr.kx_deriv, tr.ky_deriv,
                                      inv(T(grid.Nx * grid.Ny)); ndrange = size(dû.ux))
    return dû
end
