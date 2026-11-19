# GLOBAL and SPECTRAL diagnostics and injection rates of the incompressible model (WP-8, ADR-0009).
# Sums run over the stored modes with the half-plane weight w (registry.jl); the state is truncated,
# so the derivative wavenumbers equal the grid wavenumbers on every nonzero mode.

const IncompressibleProblem = Problem{<:Grid,<:Incompressible}

# Re[a*·b] of two velocity modes.
@inline velocity_product(ax, ay, bx, by) = real(conj(ax) * bx + conj(ay) * by)

# Re[ω̂_a* ω̂_b] with the vorticity ω̂_a = i(kx ay − ky ax) of the velocity mode a.
@inline vorticity_product(kx, ky, ax, ay, bx, by) = real(conj(kx * ay - ky * ax) * (kx * by - ky * bx))

"""
    kinetic_energy(state, problem, workspace) -> T

The kinetic energy, GLOBAL diagnostic `E_kin`. In the incompressible model (specification,
incompressible model, invariants):

```math
E_{\\rm kin} = \\tfrac12\\langle|\\mathbf{u}|^2\\rangle = \\tfrac12 \\sum_{\\mathbf{k}} |\\hat{\\mathbf{u}}(\\mathbf{k})|^2 .
```
"""
kinetic_energy(state, ::IncompressibleProblem, ws) =
    mode_sum((kx, ky, ax, ay, bx, by) -> velocity_product(ax, ay, bx, by) / 2, state.fields, state.fields, ws)

"""
    enstrophy(state, problem, workspace) -> T

The enstrophy, GLOBAL diagnostic `Z`, in both models (specification, invariants). The vorticity is
band-limited, so the sum over the stored modes equals the grid mean:

```math
Z = \\tfrac12\\langle\\omega^2\\rangle = \\tfrac12 \\sum_{\\mathbf{k}} |\\hat\\omega(\\mathbf{k})|^2,
\\qquad \\hat\\omega = i(k_x \\hat u_y - k_y \\hat u_x) .
```
"""
enstrophy(state, ::Problem, ws) =
    mode_sum((kx, ky, ax, ay, bx, by) -> vorticity_product(kx, ky, ax, ay, bx, by) / 2, state.fields, state.fields, ws)

"""
    palinstrophy(state, problem, workspace) -> T

The palinstrophy, GLOBAL diagnostic `P`, in both models (specification, incompressible model, global
budgets):

```math
P = \\tfrac12\\langle|\\nabla\\omega|^2\\rangle = \\tfrac12 \\sum_{\\mathbf{k}} k^2 |\\hat\\omega(\\mathbf{k})|^2 .
```
"""
palinstrophy(state, ::Problem, ws) =
    mode_sum((kx, ky, ax, ay, bx, by) -> (kx^2 + ky^2) * vorticity_product(kx, ky, ax, ay, bx, by) / 2,
             state.fields, state.fields, ws)

"""
    small_scale_energy_dissipation(state, problem, workspace) -> T

The small-scale dissipation rate of ``E_{\\rm kin}``, GLOBAL diagnostic `dissS_E_kin`. In the
incompressible model (specification, incompressible model, global budgets):

```math
\\nu \\sum_{\\mathbf{k}} k^2 |\\hat{\\mathbf{u}}(\\mathbf{k})|^2 = 2\\nu Z .
```
"""
function small_scale_energy_dissipation(state, problem::IncompressibleProblem, ws)
    ν = problem.model.ν
    return mode_sum((kx, ky, ax, ay, bx, by) -> ν * (kx^2 + ky^2) * velocity_product(ax, ay, bx, by),
                    state.fields, state.fields, ws)
end

"""
    small_scale_enstrophy_dissipation(state, problem, workspace) -> T

The small-scale dissipation rate of ``Z``, GLOBAL diagnostic `dissS_Z`. In the incompressible model
(specification, incompressible model, global budgets):

```math
\\nu \\sum_{\\mathbf{k}} k^2 |\\hat\\omega(\\mathbf{k})|^2 = 2\\nu P .
```
"""
function small_scale_enstrophy_dissipation(state, problem::IncompressibleProblem, ws)
    ν = problem.model.ν
    return mode_sum((kx, ky, ax, ay, bx, by) -> ν * (kx^2 + ky^2) * vorticity_product(kx, ky, ax, ay, bx, by),
                    state.fields, state.fields, ws)
end

"""
    large_scale_energy_dissipation(state, problem, workspace) -> T

The large-scale dissipation rate of ``E_{\\rm kin}``, GLOBAL diagnostic `dissL_E_kin`, for the
[`LargeScaleDissipation`](@ref) of the model (ADR-0004):

```math
\\alpha \\sum_{\\mathbf{k} \\neq 0} k^{-2m} |\\hat{\\mathbf{u}}(\\mathbf{k})|^2 .
```
"""
function large_scale_energy_dissipation(state, problem::IncompressibleProblem, ws)
    ls = problem.model.large_scale
    return mode_sum((kx, ky, ax, ay, bx, by) -> large_scale_rate(ls, kx^2 + ky^2) * velocity_product(ax, ay, bx, by),
                    state.fields, state.fields, ws)
end

"""
    large_scale_enstrophy_dissipation(state, problem, workspace) -> T

The large-scale dissipation rate of ``Z``, GLOBAL diagnostic `dissL_Z`, for the
[`LargeScaleDissipation`](@ref) of the model (ADR-0004):

```math
\\alpha \\sum_{\\mathbf{k} \\neq 0} k^{-2m} |\\hat\\omega(\\mathbf{k})|^2 .
```
"""
function large_scale_enstrophy_dissipation(state, problem::IncompressibleProblem, ws)
    ls = problem.model.large_scale
    return mode_sum((kx, ky, ax, ay, bx, by) -> large_scale_rate(ls, kx^2 + ky^2) *
                                               vorticity_product(kx, ky, ax, ay, bx, by),
                    state.fields, state.fields, ws)
end

"""
    kinetic_energy_spectrum_inc(state, problem, workspace) -> Vector{T}

The shell spectrum of the solenoidal kinetic energy, SPECTRAL diagnostic `E_kin_inc_k`, at the shell
centres ``k_n`` (`grid.shell_k`). In the incompressible model the velocity is solenoidal
(specification, incompressible model, spectra):

```math
E(k_n) = \\tfrac12 \\sum_{\\mathbf{q} \\in k_n} |\\hat{\\mathbf{u}}(\\mathbf{q})|^2 .
```
"""
kinetic_energy_spectrum_inc(state, ::IncompressibleProblem, ws) =
    shell_sums((kx, ky, ax, ay, bx, by) -> velocity_product(ax, ay, bx, by) / 2, state.fields, state.fields, ws)

"""
    enstrophy_spectrum(state, problem, workspace) -> Vector{T}

The shell spectrum of the enstrophy, SPECTRAL diagnostic `Z_k`, at the shell centres ``k_n``
(`grid.shell_k`), in both models (specification, incompressible model, spectra):

```math
Z(k_n) = \\tfrac12 \\sum_{\\mathbf{q} \\in k_n} |\\hat\\omega(\\mathbf{q})|^2 .
```
"""
enstrophy_spectrum(state, ::Problem, ws) =
    shell_sums((kx, ky, ax, ay, bx, by) -> vorticity_product(kx, ky, ax, ay, bx, by) / 2, state.fields, state.fields, ws)

# Cumulative sums over shells of Σ w f, with b the tendency −P N̂ of the advection (WP-6).
function cumulative_flux(f, state, problem, ws)
    n̂ = advection!(ws.diagnostics.nonlinear, state.fields, problem.grid, ws)
    shells = shell_sums(f, state.fields, n̂, ws)
    return cumsum!(shells, shells)
end

"""
    kinetic_energy_flux(state, problem, workspace) -> Vector{T}

The kinetic-energy flux, SPECTRAL diagnostic `Pi_kin`, at the upper shell edges
``K_n = (n + \\tfrac12)\\Delta k`` (`grid.shell_edge`). In the incompressible model, with
``\\mathbf{N} = (\\mathbf{u}\\cdot\\nabla)\\mathbf{u}`` (specification, incompressible model, fluxes):

```math
\\Pi_E(K_n) = \\sum_{|\\mathbf{q}| < K_n} \\mathrm{Re}\\big[\\hat{\\mathbf{u}}^*(\\mathbf{q}) \\cdot
\\hat{\\mathbf{N}}(\\mathbf{q})\\big] .
```

Entry `n + 1` sums the shells `0, …, n`. It equals the specification's ``\\Pi_E(K)`` for every `K` from
the largest ``|\\mathbf{q}|`` in these shells up to ``K_n``, without ``K_n`` itself: a mode with
``|\\mathbf{q}| = K_n`` belongs to shell `n + 1` (WP-8). ``\\hat{\\mathbf{N}}`` is the advection of
[`rhs!`](@ref), also for `nonlinear = false`. The function overwrites the buffers `physical`,
`vorticity` and `diagnostics.nonlinear` of `workspace`.
"""
kinetic_energy_flux(state, problem::IncompressibleProblem, ws) =
    cumulative_flux((kx, ky, ax, ay, bx, by) -> -velocity_product(ax, ay, bx, by), state, problem, ws)

"""
    enstrophy_flux(state, problem, workspace) -> Vector{T}

The enstrophy flux, SPECTRAL diagnostic `Pi_Z`, at the upper shell edges ``K_n = (n + \\tfrac12)\\Delta k``
(`grid.shell_edge`). In the incompressible model (specification, incompressible model, fluxes):

```math
\\Pi_Z(K_n) = \\sum_{|\\mathbf{q}| < K_n} \\mathrm{Re}\\big[\\hat\\omega^*(\\mathbf{q})\\,
\\widehat{\\mathbf{u}\\cdot\\nabla\\omega}(\\mathbf{q})\\big] .
```

``\\widehat{\\mathbf{u}\\cdot\\nabla\\omega}`` is the curl of ``\\hat{\\mathbf{N}}``. The shells and the
buffers are those of [`kinetic_energy_flux`](@ref).
"""
enstrophy_flux(state, problem::IncompressibleProblem, ws) =
    cumulative_flux((kx, ky, ax, ay, bx, by) -> -vorticity_product(kx, ky, ax, ay, bx, by), state, problem, ws)

"""
    record_injection!(problem::Problem{<:Grid,<:Incompressible}, state, workspace, Δt) -> workspace

Store the injection rates of ``E_{\\rm kin}`` and ``Z`` for the step `Δt` as `inj_E_kin` and `inj_Z`.
With ``\\hat{\\mathbf{u}}`` and ``\\hat{\\mathbf{f}}`` at the start of the step and the forcing constant
over the step (WP-8, ADR-0015):

```math
\\mathrm{inj}_{E} = \\sum_{\\mathbf{k}} \\Big(\\mathrm{Re}\\big[\\hat{\\mathbf{u}}^* \\cdot \\hat{\\mathbf{f}}\\big]
+ \\tfrac12 \\Delta t\\, |\\hat{\\mathbf{f}}|^2\\Big), \\qquad
\\mathrm{inj}_{Z} = \\sum_{\\mathbf{k}} k^2 \\Big(\\mathrm{Re}\\big[\\hat{\\mathbf{u}}^* \\cdot \\hat{\\mathbf{f}}\\big]
+ \\tfrac12 \\Delta t\\, |\\hat{\\mathbf{f}}|^2\\Big) .
```

The sums run over the full plane: the stored modes have the half-plane weight 1 on the `kx = 0` column
and the Nyquist column, 2 elsewhere. The ``\\tfrac12\\Delta t|\\hat{\\mathbf{f}}|^2`` term removes the
``O(\\Delta t/\\tau_f)`` bias of the plain product: with forcing only, the rate equals
``(E_{\\rm kin}(t_{n+1}) - E_{\\rm kin}(t_n))/\\Delta t``.
"""
function record_injection!(problem::IncompressibleProblem, state::State{T}, ws, Δt) where {T}
    h = T(Δt) / 2
    rate(kx, ky, ux, uy, fx, fy) = velocity_product(ux, uy, fx, fy) + h * velocity_product(fx, fy, fx, fy)
    rates = ws.diagnostics.injection
    rates[:inj_E_kin] = mode_sum(rate, state.fields, ws.forcing, ws)
    rates[:inj_Z] = mode_sum((kx, ky, ux, uy, fx, fy) -> (kx^2 + ky^2) * rate(kx, ky, ux, uy, fx, fy),
                             state.fields, ws.forcing, ws)
    return ws
end

const INCOMPRESSIBLE_DEFAULTS =
    (GLOBAL = (:E_kin, :Z, :P, :inj_E_kin, :inj_Z, :dissS_E_kin, :dissS_Z, :dissL_E_kin, :dissL_Z),
     SPECTRAL = (:E_kin_inc_k, :Z_k, :Pi_kin, :Pi_Z))

default_diagnostics(class::Symbol, ::Incompressible) = (class_registry(class); INCOMPRESSIBLE_DEFAULTS[class])

register_diagnostic!(:GLOBAL, :E_kin, kinetic_energy)
register_diagnostic!(:GLOBAL, :Z, enstrophy)
register_diagnostic!(:GLOBAL, :P, palinstrophy)
register_diagnostic!(:GLOBAL, :inj_E_kin, stored_injection(:inj_E_kin))
register_diagnostic!(:GLOBAL, :inj_Z, stored_injection(:inj_Z))
register_diagnostic!(:GLOBAL, :dissS_E_kin, small_scale_energy_dissipation)
register_diagnostic!(:GLOBAL, :dissS_Z, small_scale_enstrophy_dissipation)
register_diagnostic!(:GLOBAL, :dissL_E_kin, large_scale_energy_dissipation)
register_diagnostic!(:GLOBAL, :dissL_Z, large_scale_enstrophy_dissipation)
register_diagnostic!(:SPECTRAL, :E_kin_inc_k, kinetic_energy_spectrum_inc)
register_diagnostic!(:SPECTRAL, :Z_k, enstrophy_spectrum)
register_diagnostic!(:SPECTRAL, :Pi_kin, kinetic_energy_flux; coordinate = :shell_edge)
register_diagnostic!(:SPECTRAL, :Pi_Z, enstrophy_flux; coordinate = :shell_edge)
