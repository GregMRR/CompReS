# GLOBAL and SPECTRAL diagnostics and injection rates of the compressible model (WP-11, ADR-0009, ADR-0020).
# The fields w = √ρ u and s = ω/√ρ are not band-limited. GLOBAL values are grid means of physical fields, and
# spectra and fluxes sum over all grid modes (shells up to the grid corner); by Parseval the two agree. Every
# function overwrites the buffers physical, vorticity, compressible and diagnostics of the workspace.

const CompressibleProblem = Problem{<:Grid,<:Compressible}

# Spectral operators of one mode for mode_map!: derivative wavenumbers (kx, ky), mode values a and b.
@inline curl_mode(kx, ky, a, b) = i_times(kx, b) - i_times(ky, a)
@inline divergence_mode(kx, ky, a, b) = i_times(kx, a) + i_times(ky, b)
@inline laplacian_mode(kx, ky, a, b) = -(kx^2 + ky^2) * a

# Re[a^s*·b] and Re[a^c*·b] of the Helmholtz parts of the mode a (helmholtz_split! convention) against the mode b.
@inline function solenoidal_product(kx, ky, ax, ay, bx, by)
    cx, cy = compressive_part(kx, ky, ax, ay)
    return velocity_product(ax - cx, ay - cy, bx, by)
end

@inline function compressive_product(kx, ky, ax, ay, bx, by)
    cx, cy = compressive_part(kx, ky, ax, ay)
    return velocity_product(cx, cy, bx, by)
end

# Re[a* b] of two scalar modes, passed to the reductions as the ux component of a pair.
@inline scalar_product(kx, ky, ax, ay, bx, by) = real(conj(ax) * bx)

# A scalar spectral field as the vector pair that mode_sum and shell_sums read.
scalar_pair(f̂) = (ux = f̂, uy = f̂)

@kernel function mode_map_kernel!(out, f, a, b, kx, ky)
    c, r = @index(Global, NTuple)
    @inbounds out[c, r] = f(kx[c], ky[r], a[c, r], b[c, r])
end

# out ← f(kx, ky, a, b) on every mode of the spectral fields a and b; out may be a or b.
function mode_map!(f, out, a, b, ws)
    tr = ws.transforms
    mode_map_kernel!(tr.grid.backend)(out, f, a, b, tr.kx_deriv, tr.ky_deriv; ndrange = size(out))
    return out
end

@kernel function point_map_kernel!(out, f, fields)
    I = @index(Global)
    @inbounds out[I] = f(map(a -> @inbounds(a[I]), fields)...)
end

# out ← f(fields...) at every grid point of the physical fields; out may be one of them.
function point_map!(f, out, fields...)
    point_map_kernel!(get_backend(out))(out, f, fields; ndrange = length(out))
    return out
end

# Grid mean ⟨f(fields...)⟩ of physical fields, through compressible.scratch. The sum runs over each column
# first, then over the column sums. With bounds checking on (as in Pkg.test), Base.sum of the whole array adds
# blocks of 1024 points in sequence, and its round-off breaks the 10·eps bound of the WP-11 injection test.
function point_mean(f, ws, fields...)
    a = point_map!(f, ws.compressible.scratch, fields...)
    return sum(sum(a; dims = 1)) / length(a)
end

# Physical field of f(kx, ky, a, b) on the spectral fields a and b, written into out through compressible.spectral.
physical_of!(out, f, a, b, ws) = to_physical!(out, mode_map!(f, ws.compressible.spectral, a, b, ws), ws.transforms)

# Physical u and ρ of the spectral state û, into workspace.physical.
function physical_velocity_density!(û, ws)
    tr = ws.transforms
    to_physical!(ws.physical.ux, û.ux, tr)
    to_physical!(ws.physical.uy, û.uy, tr)
    to_physical!(ws.physical.rho, û.rho, tr)
    return ws.physical
end

# Physical ω of the spectral state û, into workspace.vorticity.
physical_vorticity!(û, ws) = physical_of!(ws.vorticity, curl_mode, û.ux, û.uy, ws)

# Physical ∇·a into out, for the physical vector field a = (fx(fields...), fy(fields...)). The derivative is
# spectral on the grid field a, without truncation, so it is skew-adjoint over all grid modes.
function physical_divergence!(out, fx, fy, ws, fields...)
    (; scratch, pair) = ws.compressible
    tr = ws.transforms
    to_spectral!(pair.ux, point_map!(fx, scratch, fields...), tr)
    to_spectral!(pair.uy, point_map!(fy, scratch, fields...), tr)
    return physical_of!(out, divergence_mode, pair.ux, pair.uy, ws)
end

# ŵ = (√ρ u)^ into compressible.pair, from the physical u and ρ in workspace.physical.
function spectral_w!(ws)
    (; ux, uy, rho) = ws.physical
    (; scratch, pair) = ws.compressible
    to_spectral!(pair.ux, point_map!((x, s) -> sqrt(s) * x, scratch, ux, rho), ws.transforms)
    to_spectral!(pair.uy, point_map!((y, s) -> sqrt(s) * y, scratch, uy, rho), ws.transforms)
    return pair
end

# ŝ = (ω/√ρ)^ into compressible.pair.ux, from the physical ρ and ω in the workspace.
function spectral_s!(ws)
    (; scratch, pair) = ws.compressible
    return to_spectral!(pair.ux, point_map!((w, s) -> w / sqrt(s), scratch, ws.vorticity, ws.physical.rho),
                        ws.transforms)
end

"""
    mass(state, problem, workspace) -> T

The mean density, GLOBAL diagnostic `mass` of the compressible model, conserved by the continuity
equation also for ``\\kappa > 0`` (specification, compressible model, invariants):

```math
\\langle\\rho\\rangle = \\hat\\rho(\\mathbf{0}) .
```
"""
mass(state, ::CompressibleProblem, ws) = real(Array(view(state.fields.rho, 1:1, 1:1))[1])

"""
    kinetic_energy(state, problem::Problem{<:Grid,<:Compressible}, workspace) -> T

The kinetic energy, GLOBAL diagnostic `E_kin`, in the compressible model: the grid mean
(specification, compressible model, spectral variables)

```math
E_{\\rm kin} = \\tfrac12\\langle\\rho|\\mathbf{u}|^2\\rangle = \\tfrac12\\langle|\\mathbf{w}|^2\\rangle,
\\qquad \\mathbf{w} = \\sqrt\\rho\\,\\mathbf{u} .
```
"""
function kinetic_energy(state, ::CompressibleProblem, ws)
    p = physical_velocity_density!(state.fields, ws)
    return point_mean((x, y, s) -> s * (x * x + y * y) / 2, ws, p.ux, p.uy, p.rho)
end

"""
    kinetic_energy_inc(state, problem, workspace) -> T

The solenoidal kinetic energy, GLOBAL diagnostic `E_kin_inc` of the compressible model, a sum over all grid
modes of the solenoidal part ``\\hat{\\mathbf{w}}^s = P\\hat{\\mathbf{w}}`` of [`helmholtz_split!`](@ref)
(specification, compressible model, spectra):

```math
E_{\\rm kin}^{\\rm inc} = \\tfrac12 \\sum_{\\mathbf{k}} |\\hat{\\mathbf{w}}^s(\\mathbf{k})|^2 .
```

[`kinetic_energy_comp`](@ref) gives the compressive part; the two add up to ``E_{\\rm kin}``. On the Nyquist
column and row the derivative wavenumbers vanish, so these modes count as solenoidal in one direction.
"""
function kinetic_energy_inc(state, ::CompressibleProblem, ws)
    physical_velocity_density!(state.fields, ws)
    w = spectral_w!(ws)
    return mode_sum((kx, ky, ax, ay, bx, by) -> solenoidal_product(kx, ky, ax, ay, bx, by) / 2, w, w, ws)
end

"""
    kinetic_energy_comp(state, problem, workspace) -> T

The compressive kinetic energy, GLOBAL diagnostic `E_kin_comp` of the compressible model, with
``\\hat w^c_i = k_ik_j\\hat w_j/k^2`` (see [`kinetic_energy_inc`](@ref)):

```math
E_{\\rm kin}^{\\rm comp} = \\tfrac12 \\sum_{\\mathbf{k}} |\\hat{\\mathbf{w}}^c(\\mathbf{k})|^2 .
```
"""
function kinetic_energy_comp(state, ::CompressibleProblem, ws)
    physical_velocity_density!(state.fields, ws)
    w = spectral_w!(ws)
    return mode_sum((kx, ky, ax, ay, bx, by) -> compressive_product(kx, ky, ax, ay, bx, by) / 2, w, w, ws)
end

"""
    mean_internal_energy(state, problem, workspace) -> T

The internal energy, GLOBAL diagnostic `E_int` of the compressible model, the grid mean with the
[`internal_energy`](@ref) ``e(\\rho)`` of the equation of state (specification, compressible model,
invariants):

```math
E_{\\rm int} = \\langle\\rho\\, e(\\rho)\\rangle .
```
"""
function mean_internal_energy(state, problem::CompressibleProblem, ws)
    eos = problem.model.eos
    ρ = to_physical!(ws.physical.rho, state.fields.rho, ws.transforms)
    return point_mean(s -> s * internal_energy(eos, s), ws, ρ)
end

"""
    total_energy(state, problem, workspace) -> T

The total energy, GLOBAL diagnostic `E_tot` of the compressible model, conserved for
``\\mu = \\zeta = \\kappa = \\alpha = 0`` without forcing (specification, compressible model, invariants):

```math
E_{\\rm tot} = E_{\\rm kin} + E_{\\rm int} .
```
"""
total_energy(state, problem::CompressibleProblem, ws) =
    kinetic_energy(state, problem, ws) + mean_internal_energy(state, problem, ws)

"""
    potential_enstrophy(state, problem, workspace) -> T

The potential enstrophy, GLOBAL diagnostic `Omega` of the compressible model, the grid mean
(specification, compressible model, invariants):

```math
\\Omega = \\tfrac12\\langle\\omega^2/\\rho\\rangle = \\tfrac12\\langle s^2\\rangle, \\qquad s = \\omega/\\sqrt\\rho .
```
"""
function potential_enstrophy(state, ::CompressibleProblem, ws)
    ω = physical_vorticity!(state.fields, ws)
    ρ = to_physical!(ws.physical.rho, state.fields.rho, ws.transforms)
    return point_mean((w, s) -> w * w / (2 * s), ws, ω, ρ)
end

"""
    rms_velocity(state, problem, workspace) -> T

The root-mean-square velocity, GLOBAL diagnostic `u_rms` of the compressible model, without density weight
(ADR-0025):

```math
u_{\\rm rms} = \\langle|\\mathbf{u}|^2\\rangle^{1/2} = \\Big(\\sum_{\\mathbf{k}} |\\hat{\\mathbf{u}}(\\mathbf{k})|^2\\Big)^{1/2} .
```
"""
rms_velocity(state, ::CompressibleProblem, ws) =
    sqrt(mode_sum((kx, ky, ax, ay, bx, by) -> velocity_product(ax, ay, ax, ay), state.fields, state.fields, ws))

"""
    mach_number(state, problem, workspace) -> T

The Mach number, GLOBAL diagnostic `Ma` of the compressible model, with the [`sound_speed²`](@ref) of the
equation of state at the reference density ``\\rho_0 = 1`` (specification, equation of state; ADR-0025):

```math
\\mathrm{Ma} = u_{\\rm rms}/c_s(\\rho_0) .
```
"""
mach_number(state::State{T}, problem::CompressibleProblem, ws) where {T} =
    rms_velocity(state, problem, ws) / sqrt(sound_speed²(problem.model.eos, one(T)))

"""
    pressure_dilatation(state, problem, workspace) -> T

The pressure–dilatation term, GLOBAL diagnostic `p_theta` of the compressible model: the rate at which the
pressure moves energy from ``E_{\\rm int}`` to ``E_{\\rm kin}`` (specification, compressible model, global
budgets), with the [`pressure`](@ref) of the equation of state:

```math
\\langle p\\,\\theta\\rangle, \\qquad \\theta = \\partial_i u_i .
```
"""
function pressure_dilatation(state, problem::CompressibleProblem, ws)
    eos = problem.model.eos
    û = state.fields
    θ = physical_of!(ws.compressible.viscous_x, divergence_mode, û.ux, û.uy, ws)
    ρ = to_physical!(ws.physical.rho, û.rho, ws.transforms)
    return point_mean((s, d) -> pressure(eos, s) * d, ws, ρ, θ)
end

"""
    small_scale_energy_dissipation(state, problem::Problem{<:Grid,<:Compressible}, workspace) -> T

The small-scale dissipation rate of ``E_{\\rm kin}``, GLOBAL diagnostic `dissS_E_kin`, in the compressible
model: the mean viscous dissipation ``\\langle D\\rangle``, ``D = \\tau_{ij}\\partial_j u_i``, which in a
periodic domain reads (specification, compressible model, local balances)

```math
\\langle D\\rangle = \\mu\\langle\\partial_j u_i\\,\\partial_j u_i\\rangle + \\zeta\\langle\\theta^2\\rangle
= \\sum_{\\mathbf{k}} \\Big(\\mu k^2 |\\hat{\\mathbf{u}}|^2 + \\zeta |\\mathbf{k}\\cdot\\hat{\\mathbf{u}}|^2\\Big) .
```
"""
function small_scale_energy_dissipation(state, problem::CompressibleProblem, ws)
    (; μ, ζ) = problem.model
    return mode_sum((kx, ky, ax, ay, bx, by) -> μ * (kx^2 + ky^2) * velocity_product(ax, ay, ax, ay) +
                                               ζ * abs2(kx * ax + ky * ay),
                    state.fields, state.fields, ws)
end

"""
    small_scale_potential_enstrophy_dissipation(state, problem, workspace) -> T

The small-scale dissipation rate of ``\\Omega``, GLOBAL diagnostic `dissS_Omega` of the compressible model,
with ``q = \\omega/\\rho`` and ``\\mathbf{g} = \\rho^{-1}(\\mu\\nabla^2\\mathbf{u} + \\zeta\\nabla\\theta)``
(specification, compressible model, global budgets):

```math
-\\langle q \\operatorname{curl}\\mathbf{g}\\rangle .
```

The curl is spectral on the grid values of ``\\mathbf{g}``. The rate has no fixed sign; for uniform
``\\rho = \\rho_0`` it equals ``\\mu\\langle|\\nabla\\omega|^2\\rangle/\\rho_0^2``.
"""
function small_scale_potential_enstrophy_dissipation(state, problem::CompressibleProblem, ws)
    (; μ, ζ) = problem.model
    û = state.fields
    tr = ws.transforms
    (; viscous_x, viscous_y, spectral, pair) = ws.compressible
    ω, ρ = ws.vorticity, ws.physical.rho
    vorticity_viscous_kernel!(tr.grid.backend)(spectral, pair.ux, pair.uy, û.ux, û.uy, tr.kx_deriv, tr.ky_deriv,
                                               μ, ζ; ndrange = size(spectral))
    to_physical!(ω, spectral, tr)
    to_physical!(viscous_x, pair.ux, tr)
    to_physical!(viscous_y, pair.uy, tr)
    to_physical!(ρ, û.rho, tr)
    to_spectral!(pair.ux, point_map!(/, viscous_x, viscous_x, ρ), tr)
    to_spectral!(pair.uy, point_map!(/, viscous_y, viscous_y, ρ), tr)
    curl_g = physical_of!(viscous_x, curl_mode, pair.ux, pair.uy, ws)
    return -point_mean((w, s, c) -> w * c / s, ws, ω, ρ, curl_g)
end

"""
    large_scale_energy_dissipation(state, problem::Problem{<:Grid,<:Compressible}, workspace) -> T

The large-scale dissipation rate of ``E_{\\rm kin}``, GLOBAL diagnostic `dissL_E_kin`, in the compressible
model, for the [`LargeScaleDissipation`](@ref) ``-\\alpha(-\\nabla^2)^{-m}\\mathbf{u}`` of the velocity
equation (ADR-0004):

```math
\\alpha\\,\\big\\langle\\rho\\,\\mathbf{u}\\cdot(-\\nabla^2)^{-m}\\mathbf{u}\\big\\rangle .
```

The operator ``(-\\nabla^2)^{-m}`` skips the mode ``\\mathbf{k} = 0``.
"""
function large_scale_energy_dissipation(state, problem::CompressibleProblem, ws)
    ls = problem.model.large_scale
    û = state.fields
    (; viscous_x, viscous_y) = ws.compressible
    damp(kx, ky, a, b) = large_scale_rate(ls, kx^2 + ky^2) * a
    physical_of!(viscous_x, damp, û.ux, û.ux, ws)
    physical_of!(viscous_y, damp, û.uy, û.uy, ws)
    p = physical_velocity_density!(û, ws)
    return point_mean((x, y, s, a, b) -> s * (x * a + y * b), ws, p.ux, p.uy, p.rho, viscous_x, viscous_y)
end

"""
    large_scale_potential_enstrophy_dissipation(state, problem, workspace) -> T

The large-scale dissipation rate of ``\\Omega``, GLOBAL diagnostic `dissL_Omega` of the compressible model,
for the [`LargeScaleDissipation`](@ref) of the velocity equation (ADR-0004), with ``q = \\omega/\\rho``:

```math
\\alpha\\,\\big\\langle q\\,(-\\nabla^2)^{-m}\\omega\\big\\rangle .
```
"""
function large_scale_potential_enstrophy_dissipation(state, problem::CompressibleProblem, ws)
    ls = problem.model.large_scale
    û = state.fields
    ω = physical_vorticity!(û, ws)
    damped = physical_of!(ws.compressible.viscous_x,
                          (kx, ky, a, b) -> large_scale_rate(ls, kx^2 + ky^2) * curl_mode(kx, ky, a, b), û.ux, û.uy, ws)
    ρ = to_physical!(ws.physical.rho, û.rho, ws.transforms)
    return point_mean((w, s, l) -> w * l / s, ws, ω, ρ, damped)
end

"""
    density_diffusion_kinetic_energy(state, problem, workspace) -> T

The density-diffusion rate of ``E_{\\rm kin}``, GLOBAL diagnostic `kappa_E_kin` of the compressible model:
the change of ``E_{\\rm kin}`` per unit time by the density diffusivity ``\\kappa`` (specification,
compressible model, global budgets). It has no fixed sign:

```math
\\tfrac12\\kappa\\,\\langle|\\mathbf{u}|^2\\nabla^2\\rho\\rangle .
```
"""
function density_diffusion_kinetic_energy(state, problem::CompressibleProblem, ws)
    κ = problem.model.κ
    û = state.fields
    lap = physical_of!(ws.compressible.viscous_x, laplacian_mode, û.rho, û.rho, ws)
    p = physical_velocity_density!(û, ws)
    return point_mean((x, y, l) -> κ * (x * x + y * y) * l / 2, ws, p.ux, p.uy, lap)
end

"""
    density_diffusion_internal_energy(state, problem, workspace) -> T

The density-diffusion rate of ``E_{\\rm int}``, GLOBAL diagnostic `kappa_E_int` of the compressible model,
with the [`enthalpy`](@ref) ``h(\\rho)`` of the equation of state (specification, compressible model, local
balances):

```math
\\kappa\\,\\langle h(\\rho)\\,\\nabla^2\\rho\\rangle = -\\kappa\\,\\langle h'(\\rho)\\,|\\nabla\\rho|^2\\rangle \\le 0 .
```

The code evaluates the left form, the exact rate of the semi-discrete system (ADR-0025).
"""
function density_diffusion_internal_energy(state, problem::CompressibleProblem, ws)
    (; κ, eos) = problem.model
    û = state.fields
    lap = physical_of!(ws.compressible.viscous_x, laplacian_mode, û.rho, û.rho, ws)
    ρ = to_physical!(ws.physical.rho, û.rho, ws.transforms)
    return point_mean((s, l) -> κ * enthalpy(eos, s) * l, ws, ρ, lap)
end

"""
    density_diffusion_potential_enstrophy(state, problem, workspace) -> T

The density-diffusion rate of ``\\Omega``, GLOBAL diagnostic `kappa_Omega` of the compressible model, with
``q = \\omega/\\rho`` (specification, compressible model, global budgets). It has no fixed sign:

```math
-\\tfrac12\\kappa\\,\\langle q^2\\nabla^2\\rho\\rangle .
```
"""
function density_diffusion_potential_enstrophy(state, problem::CompressibleProblem, ws)
    κ = problem.model.κ
    û = state.fields
    lap = physical_of!(ws.compressible.viscous_x, laplacian_mode, û.rho, û.rho, ws)
    ω = physical_vorticity!(û, ws)
    ρ = to_physical!(ws.physical.rho, û.rho, ws.transforms)
    return point_mean((w, s, l) -> -κ * w * w * l / (2 * s * s), ws, ω, ρ, lap)
end

"""
    kinetic_energy_spectrum_inc(state, problem::Problem{<:Grid,<:Compressible}, workspace) -> Vector{T}

The shell spectrum of the solenoidal kinetic energy, SPECTRAL diagnostic `E_kin_inc_k`, in the compressible
model, at the shell centres ``k_n`` (`grid.shell_k`), over all grid modes (specification, compressible model,
spectra; see [`kinetic_energy_inc`](@ref)):

```math
E_{\\rm kin}^{\\rm inc}(k_n) = \\tfrac12 \\sum_{\\mathbf{q} \\in k_n} |\\hat{\\mathbf{w}}^s(\\mathbf{q})|^2 .
```
"""
function kinetic_energy_spectrum_inc(state, ::CompressibleProblem, ws)
    physical_velocity_density!(state.fields, ws)
    w = spectral_w!(ws)
    return shell_sums((kx, ky, ax, ay, bx, by) -> solenoidal_product(kx, ky, ax, ay, bx, by) / 2, w, w, ws)
end

"""
    kinetic_energy_spectrum_comp(state, problem, workspace) -> Vector{T}

The shell spectrum of the compressive kinetic energy, SPECTRAL diagnostic `E_kin_comp_k` of the
compressible model, at the shell centres ``k_n`` (`grid.shell_k`), over all grid modes:

```math
E_{\\rm kin}^{\\rm comp}(k_n) = \\tfrac12 \\sum_{\\mathbf{q} \\in k_n} |\\hat{\\mathbf{w}}^c(\\mathbf{q})|^2 .
```
"""
function kinetic_energy_spectrum_comp(state, ::CompressibleProblem, ws)
    physical_velocity_density!(state.fields, ws)
    w = spectral_w!(ws)
    return shell_sums((kx, ky, ax, ay, bx, by) -> compressive_product(kx, ky, ax, ay, bx, by) / 2, w, w, ws)
end

"""
    potential_enstrophy_spectrum(state, problem, workspace) -> Vector{T}

The shell spectrum of the potential enstrophy, SPECTRAL diagnostic `Omega_k` of the compressible model, at
the shell centres ``k_n`` (`grid.shell_k`), over all grid modes (specification, compressible model, spectra):

```math
\\Omega(k_n) = \\tfrac12 \\sum_{\\mathbf{q} \\in k_n} |\\hat s(\\mathbf{q})|^2, \\qquad s = \\omega/\\sqrt\\rho .
```
"""
function potential_enstrophy_spectrum(state, ::CompressibleProblem, ws)
    physical_vorticity!(state.fields, ws)
    to_physical!(ws.physical.rho, state.fields.rho, ws.transforms)
    ŝ = scalar_pair(spectral_s!(ws))
    return shell_sums((kx, ky, ax, ay, bx, by) -> scalar_product(kx, ky, ax, ay, bx, by) / 2, ŝ, ŝ, ws)
end

"""
    density_spectrum(state, problem, workspace) -> Vector{T}

The shell spectrum of the density, SPECTRAL diagnostic `rho_k` of the compressible model, at the shell
centres ``k_n`` (`grid.shell_k`) (specification, conventions). Shell 0 holds the mean,
``\\tfrac12\\langle\\rho\\rangle^2`` (ADR-0025):

```math
\\varrho(k_n) = \\tfrac12 \\sum_{\\mathbf{q} \\in k_n} |\\hat\\rho(\\mathbf{q})|^2 .
```
"""
function density_spectrum(state, ::CompressibleProblem, ws)
    ρ̂ = scalar_pair(state.fields.rho)
    return shell_sums((kx, ky, ax, ay, bx, by) -> scalar_product(kx, ky, ax, ay, bx, by) / 2, ρ̂, ρ̂, ws)
end

# ŵ into compressible.pair and N̂^w into diagnostics.nonlinear, with the Lamb form of the advection
# (specification, note 2): N^w = √ρ (∇(½|u|²) + ω ẑ × u) + u ∇·(ρu)/(2√ρ), ω ẑ × u = (−ω u_y, ω u_x).
function kinetic_energy_transfer!(state, ws)
    û = state.fields
    tr = ws.transforms
    (; ux, uy, rho) = ws.physical
    (; viscous_x, viscous_y, spectral, scratch, pair) = ws.compressible
    ω = physical_vorticity!(û, ws)
    physical_velocity_density!(û, ws)
    to_spectral!(spectral, point_map!((x, y) -> (x * x + y * y) / 2, scratch, ux, uy), tr)
    to_physical!(viscous_x, mode_map!((kx, ky, a, b) -> i_times(kx, a), pair.ux, spectral, spectral, ws), tr)
    to_physical!(viscous_y, mode_map!((kx, ky, a, b) -> i_times(ky, a), pair.uy, spectral, spectral, ws), tr)
    div_ρu = physical_divergence!(scratch, (x, y, s) -> s * x, (x, y, s) -> s * y, ws, ux, uy, rho)
    point_map!((gx, x, y, s, w, d) -> sqrt(s) * (gx - w * y) + x * d / (2 * sqrt(s)), viscous_x,
               viscous_x, ux, uy, rho, ω, div_ρu)
    point_map!((gy, x, y, s, w, d) -> sqrt(s) * (gy + w * x) + y * d / (2 * sqrt(s)), viscous_y,
               viscous_y, ux, uy, rho, ω, div_ρu)
    n̂ = ws.diagnostics.nonlinear
    to_spectral!(n̂.ux, viscous_x, tr)
    to_spectral!(n̂.uy, viscous_y, tr)
    return spectral_w!(ws), n̂
end

# Cumulative sums over shells of Σ w f(ŵ, N̂^w), over all grid modes.
function compressible_energy_flux(f, state, ws)
    ŵ, n̂ = kinetic_energy_transfer!(state, ws)
    shells = shell_sums(f, ŵ, n̂, ws)
    return cumsum!(shells, shells)
end

"""
    kinetic_energy_flux(state, problem::Problem{<:Grid,<:Compressible}, workspace) -> Vector{T}

The kinetic-energy flux, SPECTRAL diagnostic `Pi_kin`, in the compressible model, at the upper shell edges
``K_n = (n + \\tfrac12)\\Delta k`` (`grid.shell_edge`), over all grid modes (specification, compressible
model, fluxes):

```math
\\Pi_K(K_n) = \\sum_{|\\mathbf{q}| < K_n} \\mathrm{Re}\\big[\\hat{\\mathbf{w}}^* \\cdot \\hat{\\mathbf{N}}^w\\big],
\\qquad N^w_i = \\sqrt\\rho\\,u_j\\partial_j u_i + \\frac{u_i}{2\\sqrt\\rho}\\,\\partial_j(\\rho u_j) .
```

The advection takes the Lamb form ``u_j\\partial_j u_i = \\partial_i(\\tfrac12|\\mathbf{u}|^2) +
(\\omega\\,\\hat{\\mathbf{z}}\\times\\mathbf{u})_i``, with spectral derivatives of the grid fields
``\\tfrac12|\\mathbf{u}|^2`` and ``\\rho\\mathbf{u}`` and no truncation. The flux over all grid modes, the
last entry, is then zero to round-off (specification, note 2). The shells follow
[`kinetic_energy_flux`](@ref) of the incompressible model.
"""
kinetic_energy_flux(state, ::CompressibleProblem, ws) =
    compressible_energy_flux((kx, ky, ax, ay, bx, by) -> velocity_product(ax, ay, bx, by), state, ws)

"""
    kinetic_energy_flux_inc(state, problem, workspace) -> Vector{T}

The solenoidal part of the kinetic-energy flux, SPECTRAL diagnostic `Pi_kin_inc` of the compressible model,
at the upper shell edges ``K_n`` (`grid.shell_edge`), with ``\\hat{\\mathbf{w}}^s`` of
[`kinetic_energy_inc`](@ref) and ``\\hat{\\mathbf{N}}^w`` of [`kinetic_energy_flux`](@ref):

```math
\\Pi_K^{\\rm inc}(K_n) = \\sum_{|\\mathbf{q}| < K_n} \\mathrm{Re}\\big[\\hat{\\mathbf{w}}^{s*} \\cdot \\hat{\\mathbf{N}}^w\\big],
\\qquad \\Pi_K = \\Pi_K^{\\rm inc} + \\Pi_K^{\\rm comp} .
```
"""
kinetic_energy_flux_inc(state, ::CompressibleProblem, ws) = compressible_energy_flux(solenoidal_product, state, ws)

"""
    kinetic_energy_flux_comp(state, problem, workspace) -> Vector{T}

The compressive part of the kinetic-energy flux, SPECTRAL diagnostic `Pi_kin_comp` of the compressible
model, at the upper shell edges ``K_n`` (`grid.shell_edge`) (see [`kinetic_energy_flux_inc`](@ref)):

```math
\\Pi_K^{\\rm comp}(K_n) = \\sum_{|\\mathbf{q}| < K_n} \\mathrm{Re}\\big[\\hat{\\mathbf{w}}^{c*} \\cdot \\hat{\\mathbf{N}}^w\\big] .
```
"""
kinetic_energy_flux_comp(state, ::CompressibleProblem, ws) = compressible_energy_flux(compressive_product, state, ws)

"""
    potential_enstrophy_flux(state, problem, workspace) -> Vector{T}

The potential-enstrophy flux, SPECTRAL diagnostic `Pi_Omega` of the compressible model, at the upper shell
edges ``K_n = (n + \\tfrac12)\\Delta k`` (`grid.shell_edge`), over all grid modes (specification,
compressible model, fluxes):

```math
\\Pi_\\Omega(K_n) = \\sum_{|\\mathbf{q}| < K_n} \\mathrm{Re}\\big[\\hat s^*\\,\\hat N^s\\big], \\qquad
N^s = \\frac{\\partial_j(u_j\\omega)}{\\sqrt\\rho} - \\frac{\\omega}{2\\rho^{3/2}}\\,\\partial_j(\\rho u_j) .
```

The derivatives are spectral on the grid fields ``u_j\\omega`` and ``\\rho u_j``, without truncation. The
non-polynomial factors leave an aliasing floor at ``k_{\\max}`` (specification, note 2).
"""
function potential_enstrophy_flux(state, ::CompressibleProblem, ws)
    û = state.fields
    (; ux, uy, rho) = ws.physical
    (; viscous_x, viscous_y, scratch) = ws.compressible
    ω = physical_vorticity!(û, ws)
    physical_velocity_density!(û, ws)
    div_ρu = physical_divergence!(viscous_x, (x, y, s) -> s * x, (x, y, s) -> s * y, ws, ux, uy, rho)
    div_uω = physical_divergence!(viscous_y, (x, y, w) -> x * w, (x, y, w) -> y * w, ws, ux, uy, ω)
    n̂ = to_spectral!(ws.diagnostics.nonlinear.ux,
                     point_map!((s, w, d, e) -> e / sqrt(s) - w * d / (2 * s * sqrt(s)), scratch, rho, ω, div_ρu,
                                div_uω), ws.transforms)
    ŝ = spectral_s!(ws)
    shells = shell_sums(scalar_product, scalar_pair(ŝ), scalar_pair(n̂), ws)
    return cumsum!(shells, shells)
end

"""
    record_injection!(problem::Problem{<:Grid,<:Compressible}, state, workspace, Δt) -> workspace

Store the injection rates of ``E_{\\rm kin}`` and ``\\Omega`` for the step `Δt` as `inj_E_kin` and
`inj_Omega`. With ``\\mathbf{u}``, ``\\rho`` and the forcing ``\\mathbf{f}`` at the start of the step, the
forcing constant over the step and ``q = \\omega/\\rho`` (specification, compressible model, global budgets;
WP-11):

```math
\\mathrm{inj}_{E} = \\langle\\rho\\,\\mathbf{u}\\cdot\\mathbf{f}\\rangle + \\tfrac12\\Delta t\\,\\langle\\rho|\\mathbf{f}|^2\\rangle,
\\qquad
\\mathrm{inj}_{\\Omega} = \\langle q \\operatorname{curl}\\mathbf{f}\\rangle
+ \\tfrac12\\Delta t\\,\\langle(\\operatorname{curl}\\mathbf{f})^2/\\rho\\rangle .
```

The means are grid means of the physical fields. The ``\\tfrac12\\Delta t`` terms remove the
``O(\\Delta t/\\tau_f)`` bias of the plain products: with forcing only and fixed ``\\rho``, the rates equal
``(E_{\\rm kin}(t_{n+1}) - E_{\\rm kin}(t_n))/\\Delta t`` and the same for ``\\Omega``.
"""
function record_injection!(problem::CompressibleProblem, state::State{T}, ws, Δt) where {T}
    h = T(Δt) / 2
    û = state.fields
    tr = ws.transforms
    fx, fy = ws.compressible.viscous_x, ws.compressible.viscous_y
    p = physical_velocity_density!(û, ws)
    ω = physical_vorticity!(û, ws)
    to_physical!(fx, ws.forcing.ux, tr)
    to_physical!(fy, ws.forcing.uy, tr)
    rates = ws.diagnostics.injection
    rates[:inj_E_kin] = point_mean((x, y, s, a, b) -> s * (x * a + y * b + h * (a * a + b * b)), ws,
                                   p.ux, p.uy, p.rho, fx, fy)
    curl_f = physical_of!(fx, curl_mode, ws.forcing.ux, ws.forcing.uy, ws)
    rates[:inj_Omega] = point_mean((w, s, c) -> (w * c + h * c * c) / s, ws, ω, p.rho, curl_f)
    return ws
end

const COMPRESSIBLE_DEFAULTS =
    (GLOBAL = (:mass, :E_kin, :E_kin_inc, :E_kin_comp, :E_int, :E_tot, :Z, :Omega, :u_rms, :Ma, :p_theta,
               :inj_E_kin, :inj_Omega, :dissS_E_kin, :dissS_Omega, :dissL_E_kin, :dissL_Omega,
               :kappa_E_kin, :kappa_E_int, :kappa_Omega),
     SPECTRAL = (:E_kin_inc_k, :E_kin_comp_k, :Omega_k, :rho_k, :Pi_kin, :Pi_kin_inc, :Pi_kin_comp, :Pi_Omega))

default_diagnostics(class::Symbol, ::Compressible) = (class_registry(class); COMPRESSIBLE_DEFAULTS[class])

register_diagnostic!(:GLOBAL, :mass, mass)
register_diagnostic!(:GLOBAL, :E_kin_inc, kinetic_energy_inc)
register_diagnostic!(:GLOBAL, :E_kin_comp, kinetic_energy_comp)
register_diagnostic!(:GLOBAL, :E_int, mean_internal_energy)
register_diagnostic!(:GLOBAL, :E_tot, total_energy)
register_diagnostic!(:GLOBAL, :Omega, potential_enstrophy)
register_diagnostic!(:GLOBAL, :u_rms, rms_velocity)
register_diagnostic!(:GLOBAL, :Ma, mach_number)
register_diagnostic!(:GLOBAL, :p_theta, pressure_dilatation)
register_diagnostic!(:GLOBAL, :inj_Omega, stored_injection(:inj_Omega))
register_diagnostic!(:GLOBAL, :dissS_Omega, small_scale_potential_enstrophy_dissipation)
register_diagnostic!(:GLOBAL, :dissL_Omega, large_scale_potential_enstrophy_dissipation)
register_diagnostic!(:GLOBAL, :kappa_E_kin, density_diffusion_kinetic_energy)
register_diagnostic!(:GLOBAL, :kappa_E_int, density_diffusion_internal_energy)
register_diagnostic!(:GLOBAL, :kappa_Omega, density_diffusion_potential_enstrophy)
register_diagnostic!(:SPECTRAL, :E_kin_comp_k, kinetic_energy_spectrum_comp)
register_diagnostic!(:SPECTRAL, :Omega_k, potential_enstrophy_spectrum)
register_diagnostic!(:SPECTRAL, :rho_k, density_spectrum)
register_diagnostic!(:SPECTRAL, :Pi_kin_inc, kinetic_energy_flux_inc; coordinate = :shell_edge)
register_diagnostic!(:SPECTRAL, :Pi_kin_comp, kinetic_energy_flux_comp; coordinate = :shell_edge)
register_diagnostic!(:SPECTRAL, :Pi_Omega, potential_enstrophy_flux; coordinate = :shell_edge)
