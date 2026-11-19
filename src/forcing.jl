# Forcing supports, the Kolmogorov, OU and random-phase forcings, and their per-step update (WP-9, ADR-0007).

"""
    support(pred)

The forcing support of the wavevectors `\\mathbf{k}` with `pred(kx, ky) == true`, in physical
wavenumbers `k_x = 2\\pi i/L_x`, `k_y = 2\\pi j/L_y` of the float type of the grid.

A forcing is a real field, so its support holds `-\\mathbf{k}` with every `\\mathbf{k}`: a mode is in
the support when `pred` is true for `\\mathbf{k}` or for `-\\mathbf{k}`. A forcing constructor binds
the support to its grid and throws an `ArgumentError` when the support is empty, holds the mode
`\\mathbf{k} = 0`, or holds a mode outside the Galerkin truncation. See also [`band`](@ref) and
[`modes`](@ref).
"""
support(pred) = PredicateSupport(pred)

struct PredicateSupport{P}
    predicate::P
end

"""
    band(k_start, k_end)

The forcing band `k_{\\rm start} \\le |\\mathbf{k}| \\le k_{\\rm end}`, in physical wavenumbers (not
mode indices). The bounds are inclusive, with a relative tolerance of `8\\,\\mathrm{eps}(T)` on
`|\\mathbf{k}|^2` for the round-off of the grid wavenumbers. On `L_x = 4\\pi`, `band(0.5, 1)` holds
the mode indices `(1, 0)` and `(2, 0)`. The function throws an `ArgumentError` unless
`0 ≤ k_start ≤ k_end`. See [`support`](@ref) for the rules of a forcing support.
"""
function band(k_start::Real, k_end::Real)
    0 <= k_start <= k_end || throw(ArgumentError("band needs 0 ≤ k_start ≤ k_end, got $k_start, $k_end"))
    return support((kx, ky) -> in_band(kx^2 + ky^2, k_start, k_end))
end

function in_band(k2, k_start, k_end)
    slack = 8 * eps(typeof(k2))
    return k_start^2 * (1 - slack) <= k2 <= k_end^2 * (1 + slack)
end

"""
    modes(list)

The forcing support of the listed integer mode indices `(i, j)`, the wavevectors
`(2\\pi i/L_x, 2\\pi j/L_y)`. Each entry stands for the pair `\\pm(i, j)`, and a repeated pair counts
once. See [`support`](@ref) for the rules of a forcing support.
"""
modes(list) = ModeSupport([(Int(i), Int(j)) for (i, j) in list])

struct ModeSupport
    modes::Vector{Tuple{Int,Int}}
end

# Throw an ArgumentError when the mode (i, j) cannot be forced on grid.
function check_forced_mode(grid, i, j)
    (i, j) == (0, 0) && throw(ArgumentError("the forcing support holds the mode k = 0"))
    retained(i, j, grid.Nx, grid.Ny) ||
        throw(ArgumentError("mode ($i, $j) of the forcing support lies outside the Galerkin truncation"))
    return nothing
end

# Linear positions in the real-FFT layout of the stored modes of a support, in increasing order. The
# kx = 0 column stores both (0, j) and (0, −j), so the positions hold both.
function support_positions(grid, s::PredicateSupport)
    kx, ky = Array(grid.kx), Array(grid.ky)
    positions = Int[]
    for I in CartesianIndices(kx)
        (s.predicate(kx[I], ky[I]) || s.predicate(-kx[I], -ky[I])) || continue
        check_forced_mode(grid, I[1] - 1, signed_index(I[2], grid.Ny))
        push!(positions, LinearIndices(kx)[I])
    end
    return positions
end

function support_positions(grid, s::ModeSupport)
    L = LinearIndices((grid.Nx ÷ 2 + 1, grid.Ny))
    row(j) = j < 0 ? grid.Ny + j + 1 : j + 1
    positions = Int[]
    for (i, j) in s.modes
        check_forced_mode(grid, i, j)
        a, b = i < 0 ? (-i, -j) : (i, j)
        push!(positions, L[a + 1, row(b)])
        a == 0 && push!(positions, L[1, row(-b)])
    end
    return sort!(unique!(positions))
end

# A forcing support bound to a grid: the device index list and unit wavevectors of the stored modes.
struct ForcedModes{I<:AbstractVector{Int32},V<:AbstractVector}
    index::I
    k̂x::V
    k̂y::V
    count::Int
end

# Bind the support s to grid. Returns the ForcedModes and, as host vectors, |k| and the full-plane
# weight of each stored mode (1 on the kx = 0 column, 2 elsewhere; no Nyquist mode is retained).
function bind_support(grid::Grid, s)
    positions = support_positions(grid, s)
    isempty(positions) && throw(ArgumentError("the forcing support holds no mode of the grid"))
    kx, ky = Array(grid.kx)[positions], Array(grid.ky)[positions]
    k = sqrt.(kx .^ 2 .+ ky .^ 2)
    weight = [(p - 1) % (grid.Nx ÷ 2 + 1) == 0 ? 1 : 2 for p in positions]
    on_device(a) = to_device(grid.backend, a)
    bound = ForcedModes(on_device(Int32.(positions)), on_device(kx ./ k), on_device(ky ./ k), sum(weight))
    return bound, k, weight
end

convert_float(::Type{T}, m::ForcedModes) where {T<:AbstractFloat} = ForcedModes(m.index, T.(m.k̂x), T.(m.k̂y), m.count)

# Amplitudes A = amp(|k|) √(total / Σ w amp²) of the stored modes, so that the full-plane sum Σ w A² is total.
function mode_amplitudes(amp, k::Vector{T}, weight, total) where {T}
    a = T.(amp.(k))
    norm2 = sum(weight .* a .^ 2)
    norm2 > 0 || throw(ArgumentError("the amplitude profile amp is zero on every forced mode"))
    return a .* sqrt(T(total) / norm2)
end

# Solenoidal and compressive amplitudes √(1 − χ) A and √χ A on the device, for the compressive fraction χ.
function split_amplitudes(backend, A::Vector{T}, χ) where {T}
    0 <= χ <= 1 || throw(ArgumentError("compressive_fraction must lie in [0, 1], got $χ"))
    return to_device(backend, sqrt(1 - T(χ)) .* A), to_device(backend, sqrt(T(χ)) .* A)
end

# A forcing that acts on the index list of its ForcedModes `modes`. On stored mode p it writes
#     f̂ = A_s z_s ê_s + A_c z_c ê_c,   ê_s = (−k̂y, k̂x),   ê_c = (k̂x, k̂y),
# with the amplitudes A_s = `solenoidal_amplitude[p]`, A_c = `compressive_amplitude[p]` and the complex
# coefficients z_s = `z[p]`, z_c = `z[m + p]` (m stored modes).
abstract type ModeForcing end

"""
    KolmogorovForcing(grid::Grid{T}; F, n)

Kolmogorov forcing, the single-mode shear acceleration (ADR-0007)

```math
\\mathbf{f} = F \\sin(k_f y)\\, \\hat{\\mathbf{x}}, \\qquad k_f = n\\, \\frac{2\\pi}{L_y}, \\qquad
\\hat f_x(0, \\pm k_f) = \\mp \\frac{iF}{2},
```

with support `(0, \\pm k_f)`. `F` is an acceleration (velocity per time); `n` is the mode index
along `y`, as in [`kolmogorov`](@ref). The forcing is solenoidal and has no state. The constructor
throws an `ArgumentError` for `n < 1` and for a mode outside the Galerkin truncation.
"""
struct KolmogorovForcing{T<:AbstractFloat,M<:ForcedModes,V<:AbstractVector{T},C<:AbstractVector{Complex{T}}} <: ModeForcing
    F::T
    n::Int
    modes::M
    solenoidal_amplitude::V
    compressive_amplitude::V
    z::C
end

function KolmogorovForcing(grid::Grid{T}; F::Real, n::Integer) where {T}
    n >= 1 || throw(ArgumentError("mode index n must be at least 1, not $n"))
    bound, _, _ = bind_support(grid, modes([(0, n)]))
    # Stored modes (0, n) and (0, −n), with ê_s = (−1, 0) and (1, 0): z_s = i and A_s = F/2 give the formula.
    amplitude = to_device(grid.backend, fill(T(F) / 2, 2))
    z = to_device(grid.backend, Complex{T}[im, im, 0, 0])
    return KolmogorovForcing(T(F), Int(n), bound, amplitude, zero(amplitude), z)
end

"""
    OUForcing(grid::Grid{T}, support; ε, τ_f, seed, compressive_fraction = 0, amp = Returns(1))

OU forcing on the forcing `support` ([`band`](@ref), [`modes`](@ref) or [`support`](@ref)): each
forced mode follows an independent complex Ornstein–Uhlenbeck process with correlation time `τ_f`
(ADR-0007). [`update_forcing!`](@ref) advances it once per step `Δt` with the exact update

```math
\\hat{\\mathbf{f}}(\\mathbf{k}) \\leftarrow e^{-\\Delta t/\\tau_f}\\, \\hat{\\mathbf{f}}(\\mathbf{k})
+ \\sqrt{1 - e^{-2\\Delta t/\\tau_f}}\\; \\sigma_k
\\big(\\sqrt{1-\\chi}\\, \\xi_s \\hat{\\mathbf{e}}_s + \\sqrt{\\chi}\\, \\xi_c \\hat{\\mathbf{e}}_c\\big),
\\qquad \\hat{\\mathbf{e}}_s = \\frac{(-k_y, k_x)}{|\\mathbf{k}|}, \\quad
\\hat{\\mathbf{e}}_c = \\frac{(k_x, k_y)}{|\\mathbf{k}|},
```

with `\\chi` = `compressive_fraction` and independent complex normal `\\xi_s, \\xi_c`,
`\\langle|\\xi|^2\\rangle = 1`. The stationary variance is `\\langle|\\hat{\\mathbf{f}}(\\mathbf{k})|^2\\rangle
= \\sigma_k^2`, with

```math
\\sigma_k^2 = \\frac{\\varepsilon}{\\tau_f}\\,
\\frac{\\mathrm{amp}(|\\mathbf{k}|)^2}{\\sum_{\\mathbf{q}} w_{\\mathbf{q}}\\, \\mathrm{amp}(|\\mathbf{q}|)^2},
```

where the sum runs over the stored forced modes and `w_{\\mathbf{q}}` is the half-plane weight (1 on
the `kx = 0` column, 2 elsewhere), so that `\\sum_{\\mathbf{k}} \\sigma_k^2 = \\varepsilon/\\tau_f` over
the full plane. For a velocity driven by the forcing alone, `\\partial_t\\hat{\\mathbf{u}} = \\hat{\\mathbf{f}}`,
the stationary mean injection rate is
`\\sum_{\\mathbf{k}} \\mathrm{Re}\\langle\\hat{\\mathbf{u}}^*\\cdot\\hat{\\mathbf{f}}\\rangle = \\tau_f \\sum_{\\mathbf{k}}\\sigma_k^2 = \\varepsilon`;
in a flow, this holds when `τ_f` is short against the time over which the flow changes `\\hat{\\mathbf{u}}`.
With the forcing held constant over a step, the injection `⟨u·f⟩ + ½Δt⟨|f|²⟩` of the
time-discrete process is `\\varepsilon\\, (\\Delta t/2\\tau_f) \\coth(\\Delta t/2\\tau_f) =
\\varepsilon\\,(1 + (\\Delta t/\\tau_f)^2/12 + \\dots)`.

Units: `ε` is an energy injection rate per unit mass (velocity² per time), `τ_f` a time. `amp(k)`,
a function of the physical `|\\mathbf{k}|`, is the dimensionless relative amplitude profile; only its
ratios between modes matter. The forcing starts at `\\hat{\\mathbf{f}} = 0` and reaches its
stationary variance after a few `τ_f`.

`seed` and the step number `n` seed the random numbers of each step (see [`update_forcing!`](@ref)).
A nonzero `compressive_fraction` needs the [`Compressible`](@ref) model; [`Problem`](@ref) checks it.
The constructor throws an `ArgumentError` for a bad support (see [`support`](@ref)), `ε < 0`,
`τ_f ≤ 0`, a `compressive_fraction` outside `[0, 1]`, or an `amp` that is zero on every forced mode.
"""
struct OUForcing{T<:AbstractFloat,M<:ForcedModes,V<:AbstractVector{T},C<:AbstractVector{Complex{T}}} <: ModeForcing
    ε::T
    τ_f::T
    compressive_fraction::T
    seed::Int
    modes::M
    solenoidal_amplitude::V
    compressive_amplitude::V
    z::C
    noise::C
    noise_host::Vector{Complex{T}}
end

function OUForcing(grid::Grid{T}, support; ε::Real, τ_f::Real, seed::Integer, compressive_fraction::Real = 0,
                   amp = Returns(1)) where {T}
    ε >= 0 || throw(ArgumentError("ε must be nonnegative, got $ε"))
    τ_f > 0 || throw(ArgumentError("τ_f must be positive, got $τ_f"))
    bound, k, weight = bind_support(grid, support)
    σ = mode_amplitudes(amp, k, weight, T(ε) / T(τ_f))
    a_s, a_c = split_amplitudes(grid.backend, σ, compressive_fraction)
    z = KernelAbstractions.zeros(grid.backend, Complex{T}, 2 * length(k))
    return OUForcing(T(ε), T(τ_f), T(compressive_fraction), Int(seed), bound, a_s, a_c, z, similar(z),
                     zeros(Complex{T}, length(z)))
end

"""
    RandomPhaseForcing(grid::Grid{T}, support; f_rms, τ_r, seed, compressive_fraction = 0, amp = Returns(1))

Random-phase forcing on the forcing `support` ([`band`](@ref), [`modes`](@ref) or [`support`](@ref)):
fixed amplitudes and random phases, all redrawn together every `τ_r` and constant in between
(ADR-0007). On each forced mode

```math
\\hat{\\mathbf{f}}(\\mathbf{k}) = a_k \\big(\\sqrt{1-\\chi}\\, e^{i\\varphi_s} \\hat{\\mathbf{e}}_s
+ \\sqrt{\\chi}\\, e^{i\\varphi_c} \\hat{\\mathbf{e}}_c\\big), \\qquad
a_k^2 = f_{\\rm rms}^2\\, \\frac{\\mathrm{amp}(|\\mathbf{k}|)^2}{\\sum_{\\mathbf{q}} w_{\\mathbf{q}}\\, \\mathrm{amp}(|\\mathbf{q}|)^2},
```

with `\\chi` = `compressive_fraction`, phases `\\varphi` uniform in `[0, 2\\pi)`, and the unit
vectors and half-plane weights `w_{\\mathbf{q}}` of [`OUForcing`](@ref). Then
`\\langle|\\mathbf{f}|^2\\rangle = \\sum_{\\mathbf{k}} |\\hat{\\mathbf{f}}(\\mathbf{k})|^2 = f_{\\rm rms}^2`.

Units: `f_rms` is an acceleration (velocity per time), `τ_r` a time; `amp(k)` is the dimensionless
relative amplitude profile of [`OUForcing`](@ref).

[`update_forcing!`](@ref) redraws the phases at the multiples `m τ_r` of the reset period: at the start
of the first step, and at each step that starts within `Δt/2` of the next multiple. A step length that
divides `τ_r` gives a redraw exactly every `τ_r/Δt` steps. A nonzero `compressive_fraction` needs the
[`Compressible`](@ref) model; [`Problem`](@ref) checks it. The constructor throws an `ArgumentError`
for a bad support (see [`support`](@ref)), `f_rms < 0`, `τ_r ≤ 0`, a `compressive_fraction` outside
`[0, 1]`, or an `amp` that is zero on every forced mode.
"""
mutable struct RandomPhaseForcing{T<:AbstractFloat,M<:ForcedModes,V<:AbstractVector{T},C<:AbstractVector{Complex{T}}} <: ModeForcing
    const f_rms::T
    const τ_r::T
    const compressive_fraction::T
    const seed::Int
    const modes::M
    const solenoidal_amplitude::V
    const compressive_amplitude::V
    const z::C
    const z_host::Vector{Complex{T}}
    t_reset::T
end

function RandomPhaseForcing(grid::Grid{T}, support; f_rms::Real, τ_r::Real, seed::Integer,
                            compressive_fraction::Real = 0, amp = Returns(1)) where {T}
    f_rms >= 0 || throw(ArgumentError("f_rms must be nonnegative, got $f_rms"))
    τ_r > 0 || throw(ArgumentError("τ_r must be positive, got $τ_r"))
    bound, k, weight = bind_support(grid, support)
    a = mode_amplitudes(amp, k, weight, T(f_rms)^2)
    a_s, a_c = split_amplitudes(grid.backend, a, compressive_fraction)
    z = KernelAbstractions.zeros(grid.backend, Complex{T}, 2 * length(k))
    return RandomPhaseForcing(T(f_rms), T(τ_r), T(compressive_fraction), Int(seed), bound, a_s, a_c, z,
                              zeros(Complex{T}, length(z)), T(-Inf))
end

"""
    mode_count(forcing) -> Int

The number of wavevectors in the forcing support of `forcing`, over the full wavenumber plane: a mode
and its opposite `-\\mathbf{k}` count as two.
"""
mode_count(f::ModeForcing) = f.modes.count

# SplitMix64 step (Steele, Lea and Flood, 2014): the next state and a well-mixed 64-bit output.
function splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    z = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return x, z ⊻ (z >> 31)
end

# Generator of step n for seed: a Xoshiro state from SplitMix64 over (seed, n). A restart then needs only
# (seed, n) (ADR-0007), and each step draws an independent sequence.
function step_rng(seed::Integer, n::Integer)
    _, h = splitmix64(seed % UInt64)
    x, s0 = splitmix64(h ⊻ (n % UInt64))
    x, s1 = splitmix64(x)
    x, s2 = splitmix64(x)
    _, s3 = splitmix64(x)
    return Xoshiro(s0, s1, s2, s3)
end

@inline function write_mode!(fx, fy, I, k̂x, k̂y, s, c)
    @inbounds fx[I] = c * k̂x - s * k̂y
    @inbounds fy[I] = c * k̂y + s * k̂x
    return nothing
end

@kernel function write_kernel!(fx, fy, z, a_s, a_c, index, k̂x, k̂y)
    p = @index(Global)
    m = length(index)
    @inbounds write_mode!(fx, fy, index[p], k̂x[p], k̂y[p], a_s[p] * z[p], a_c[p] * z[m + p])
end

@kernel function ou_kernel!(fx, fy, z, ξ, a_s, a_c, index, k̂x, k̂y, decay, kick)
    p = @index(Global)
    m = length(index)
    @inbounds begin
        zs = decay * z[p] + kick * ξ[p]
        zc = decay * z[m + p] + kick * ξ[m + p]
        z[p], z[m + p] = zs, zc
        write_mode!(fx, fy, index[p], k̂x[p], k̂y[p], a_s[p] * zs, a_c[p] * zc)
    end
end

@kernel function add_kernel!(dux, duy, fx, fy, index)
    p = @index(Global)
    @inbounds begin
        I = index[p]
        dux[I] += fx[I]
        duy[I] += fy[I]
    end
end

# Make the forcing buffer of the workspace Hermitian on the kx = 0 column (ADR-0015).
function enforce_hermitian_forcing!(ws)
    grid = ws.transforms.grid
    enforce_hermitian!(ws.forcing.ux, grid)
    enforce_hermitian!(ws.forcing.uy, grid)
    return nothing
end

# Write f̂ of the coefficients z of f into the forcing buffer of ws.
function write_forcing!(f::ModeForcing, ws)
    m = f.modes
    write_kernel!(ws.transforms.grid.backend)(ws.forcing.ux, ws.forcing.uy, f.z, f.solenoidal_amplitude,
                                              f.compressive_amplitude, m.index, m.k̂x, m.k̂y;
                                              ndrange = length(m.index))
    enforce_hermitian_forcing!(ws)
    return nothing
end

"""
    update_forcing!(forcing, state, workspace, Δt) -> workspace

Advance the state of `forcing` once, at the start of the step `Δt` from `state.t` (step `state.n`), and
write the forcing ``\\hat{\\mathbf{f}}`` of that step into `workspace.forcing`. The time stepping calls
it once per step, before the first stage; ``\\hat{\\mathbf{f}}`` stays constant over the step (ADR-0005).

- [`NoForcing`](@ref): no change; `workspace.forcing` stays zero.
- [`KolmogorovForcing`](@ref): writes the fixed field.
- [`OUForcing`](@ref): applies the exact OU update with step `Δt`.
- [`RandomPhaseForcing`](@ref): redraws the phases when the step starts a new period `τ_r`.

The random numbers of step `n` come from a `Random.Xoshiro` generator on the host, seeded from
`(seed, n)`; the host copies them to the device. Every backend then gets the same sequence, and a
restart needs only the seed, the step and [`forcing_state`](@ref). Every written
``\\hat{\\mathbf{f}}`` passes through [`enforce_hermitian!`](@ref). The update touches only the forced
modes of `workspace.forcing`: its cost scales with the number of forced modes.
"""
update_forcing!(::NoForcing, state, ws, Δt) = ws

function update_forcing!(f::KolmogorovForcing, state, ws, Δt)
    write_forcing!(f, ws)
    return ws
end

function update_forcing!(f::OUForcing{T}, state, ws, Δt) where {T}
    randn!(step_rng(f.seed, state.n), f.noise_host)
    copyto!(f.noise, f.noise_host)
    r = T(Δt) / f.τ_f
    m = f.modes
    ou_kernel!(ws.transforms.grid.backend)(ws.forcing.ux, ws.forcing.uy, f.z, f.noise, f.solenoidal_amplitude,
                                           f.compressive_amplitude, m.index, m.k̂x, m.k̂y, exp(-r),
                                           sqrt(-expm1(-2r)); ndrange = length(m.index))
    enforce_hermitian_forcing!(ws)
    return ws
end

function update_forcing!(f::RandomPhaseForcing{T}, state, ws, Δt) where {T}
    t = state.t + T(Δt) / 2
    if t >= f.t_reset + f.τ_r
        rng = step_rng(f.seed, state.n)
        for p in eachindex(f.z_host)
            f.z_host[p] = cis(2 * T(π) * rand(rng, T))
        end
        copyto!(f.z, f.z_host)
        f.t_reset = f.τ_r * floor(t / f.τ_r)
    end
    write_forcing!(f, ws)
    return ws
end

"""
    add_forcing!(dû, forcing::Union{KolmogorovForcing,OUForcing,RandomPhaseForcing}, state, workspace) -> dû

Add the forcing ``\\hat{\\mathbf{f}}`` of the current step, stored in `workspace.forcing` by
[`update_forcing!`](@ref), to the velocity tendency `dû.ux`, `dû.uy` on the forced modes only. A
density tendency `dû.rho` stays unchanged. The incompressible model needs a solenoidal
``\\hat{\\mathbf{f}}``: with `compressive_fraction = 0` the forcing is solenoidal by construction,
``\\mathbf{k}\\cdot\\hat{\\mathbf{f}} = 0``, so it equals its projection, and [`Problem`](@ref) rejects
a nonzero compressive fraction with the [`Incompressible`](@ref) model.
"""
function add_forcing!(dû, f::ModeForcing, state, ws)
    add_kernel!(ws.transforms.grid.backend)(dû.ux, dû.uy, ws.forcing.ux, ws.forcing.uy, f.modes.index;
                                            ndrange = length(f.modes.index))
    return dû
end

# The incompressible model rejects a compressive forcing; Problem calls this check.
check_forcing(::Incompressible, f::Union{OUForcing,RandomPhaseForcing}) =
    iszero(f.compressive_fraction) ? nothing :
    throw(ArgumentError("the incompressible model takes no compressive forcing; compressive_fraction = $(f.compressive_fraction)"))

"""
    forcing_state(forcing) -> NamedTuple

The state of `forcing` that a restart needs besides the seed and the step, as host values (ADR-0007,
ADR-0009). [`set_forcing_state!`](@ref) restores it.

- [`NoForcing`](@ref), [`KolmogorovForcing`](@ref): `(;)`, no state.
- [`OUForcing`](@ref): `(; z)`, the unit-variance OU coefficients of the stored forced modes, with
  ``\\hat{\\mathbf{f}} = \\sigma_k(\\sqrt{1-\\chi}\\, z_s \\hat{\\mathbf{e}}_s + \\sqrt{\\chi}\\, z_c \\hat{\\mathbf{e}}_c)``:
  a `Vector{Complex{T}}` with the `z_s` of the `m` stored modes followed by their `z_c`.
- [`RandomPhaseForcing`](@ref): `(; z, t_reset)`, the current phase factors `e^{i\\varphi}` in the
  same layout, and the time of the last phase reset.
"""
forcing_state(::Union{NoForcing,KolmogorovForcing}) = (;)
forcing_state(f::OUForcing) = (; z = Array(f.z))
forcing_state(f::RandomPhaseForcing) = (; z = Array(f.z), t_reset = f.t_reset)

"""
    set_forcing_state!(forcing, s::NamedTuple) -> forcing

Restore the state `s` of [`forcing_state`](@ref) into `forcing`, which has the same parameters and
support as the forcing that `s` comes from. Throws an `ArgumentError` when the length of `s.z`
differs from that of the forcing.
"""
set_forcing_state!(f::Union{NoForcing,KolmogorovForcing}, ::NamedTuple{()}) = f

function set_forcing_state!(f::OUForcing, s::NamedTuple{(:z,)})
    copy_coefficients!(f.z, s.z)
    return f
end

function set_forcing_state!(f::RandomPhaseForcing, s::NamedTuple{(:z, :t_reset)})
    copy_coefficients!(f.z, s.z)
    f.t_reset = s.t_reset
    return f
end

function copy_coefficients!(z, saved)
    length(saved) == length(z) ||
        throw(ArgumentError("the saved forcing state has $(length(saved)) coefficients, the forcing $(length(z))"))
    copyto!(z, saved)
    return z
end

convert_float(::Type{<:AbstractFloat}, f::NoForcing) = f

function convert_float(::Type{T}, f::KolmogorovForcing) where {T<:AbstractFloat}
    return KolmogorovForcing(T(f.F), f.n, convert_float(T, f.modes), T.(f.solenoidal_amplitude),
                             T.(f.compressive_amplitude), Complex{T}.(f.z))
end

function convert_float(::Type{T}, f::OUForcing) where {T<:AbstractFloat}
    z = Complex{T}.(f.z)
    return OUForcing(T(f.ε), T(f.τ_f), T(f.compressive_fraction), f.seed, convert_float(T, f.modes),
                     T.(f.solenoidal_amplitude), T.(f.compressive_amplitude), z, similar(z),
                     zeros(Complex{T}, length(z)))
end

function convert_float(::Type{T}, f::RandomPhaseForcing) where {T<:AbstractFloat}
    return RandomPhaseForcing(T(f.f_rms), T(f.τ_r), T(f.compressive_fraction), f.seed, convert_float(T, f.modes),
                              T.(f.solenoidal_amplitude), T.(f.compressive_amplitude), Complex{T}.(f.z),
                              zeros(Complex{T}, length(f.z)), T(f.t_reset))
end
