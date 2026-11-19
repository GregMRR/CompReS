# Real FFTs with the specification normalisation and the spectral operators on the real-FFT layout (WP-2, ADR-0001).

"""
    Transforms(grid::Grid)

FFT plans, derivative wavenumbers and scratch storage for one [`Grid`](@ref).
Build it once per grid and pass it to [`to_spectral!`](@ref), [`to_physical!`](@ref), [`∂x!`](@ref),
[`∂y!`](@ref), [`helmholtz_split!`](@ref) and [`project_solenoidal!`](@ref).

A physical field is a real `Nx × Ny` device array. A spectral field is a complex
`(Nx ÷ 2 + 1) × Ny` device array in the real-FFT layout of the grid.

The fields are:

- `grid`: the grid.
- `forward`, `inverse`: the real-to-complex and complex-to-real FFT plans (FFTW with
  `Threads.nthreads()` threads on the CPU, cuFFT or rocFFT on a GPU). On the CPU, `inverse` also
  holds one spectral scratch array.
- `kx_deriv` (length `Nx ÷ 2 + 1`), `ky_deriv` (length `Ny`): the derivative wavenumbers. They equal
  the grid wavenumbers, except that `kx_deriv` is zero on the Nyquist column `kx = Nx/2` and
  `ky_deriv` is zero on the Nyquist row `ky = Ny/2`. The derivative multipliers are
  `i k_x^{\\rm d}` and `i k_y^{\\rm d}`.

The Nyquist mode of a real field has no sign: it stands for both `+N/2` and `-N/2`. A nonzero
multiplier there makes the discrete derivative of a real field complex. The zero multiplier keeps
the derivative real and skew-adjoint, `\\sum f\\,\\partial_x g = -\\sum g\\,\\partial_x f`.
"""
struct Transforms{T<:AbstractFloat,G<:Grid{T},F,B,V<:AbstractVector{T}}
    grid::G
    forward::F
    inverse::B
    kx_deriv::V
    ky_deriv::V
end

# FFTW plans on the CPU use all Julia threads and the MEASURE planner; cuFFT and rocFFT take no options.
plan_forward(f::Array) = plan_rfft(f; flags = FFTW.MEASURE, num_threads = Threads.nthreads())
plan_forward(f) = plan_rfft(f)

# Inverse FFT on the CPU in two stages: a complex FFT along y from the input into `scratch`, then a
# complex-to-real FFT along x. A single 2D complex-to-real FFTW transform overwrites its input,
# and the two stages are faster than a copy followed by that transform.
struct TwoStageInverse{Y,X,S}
    along_y::Y
    along_x::X
    scratch::S
end

function plan_inverse(f̂::Array, Nx)
    scratch = similar(f̂)
    threads = Threads.nthreads()
    return TwoStageInverse(FFTW.plan_bfft(f̂, 2; flags = FFTW.MEASURE, num_threads = threads),
                           plan_brfft(scratch, Nx, 1; flags = FFTW.MEASURE, num_threads = threads), scratch)
end
# cuFFT and rocFFT copy the input of a complex-to-real transform into a buffer of the plan.
plan_inverse(f̂, Nx) = plan_brfft(f̂, Nx)

inverse!(f, p::TwoStageInverse, f̂) = mul!(f, p.along_x, mul!(p.scratch, p.along_y, f̂))
inverse!(f, p, f̂) = mul!(f, p, f̂)

function Transforms(grid::Grid{T}) where {T}
    (; Nx, Ny, backend) = grid
    f = KernelAbstractions.zeros(backend, T, Nx, Ny)
    f̂ = KernelAbstractions.zeros(backend, Complex{T}, Nx ÷ 2 + 1, Ny)
    kx = Array(grid.kx)[:, 1]
    ky = Array(grid.ky)[1, :]
    kx[end] = 0
    ky[Ny ÷ 2 + 1] = 0
    kx_deriv = KernelAbstractions.allocate(backend, T, length(kx))
    ky_deriv = KernelAbstractions.allocate(backend, T, length(ky))
    copyto!(kx_deriv, kx)
    copyto!(ky_deriv, ky)
    return Transforms(grid, plan_forward(f), plan_inverse(f̂, Nx), kx_deriv, ky_deriv)
end

@kernel function scale_kernel!(f̂, s)
    I = @index(Global)
    @inbounds f̂[I] *= s
end

"""
    to_spectral!(f̂, f, tr::Transforms) -> f̂

Write the Fourier coefficients of the real physical field `f` into the spectral field `f̂`:

```math
\\hat f(\\mathbf{k}) = \\langle f\\, e^{-i\\mathbf{k}\\cdot\\mathbf{x}} \\rangle
= \\frac{1}{N_x N_y} \\sum_{\\mathbf{x}} f(\\mathbf{x})\\, e^{-i\\mathbf{k}\\cdot\\mathbf{x}}
```

(specification, Conventions). The function leaves `f` unchanged and does not truncate `f̂`.
"""
function to_spectral!(f̂, f, tr::Transforms{T}) where {T}
    mul!(f̂, tr.forward, f)
    scale_kernel!(tr.grid.backend)(f̂, inv(T(tr.grid.Nx * tr.grid.Ny)); ndrange = length(f̂))
    return f̂
end

"""
    to_physical!(f, f̂, tr::Transforms) -> f

Write the real physical field `f(\\mathbf{x}) = \\sum_{\\mathbf{k}} \\hat f(\\mathbf{k}) e^{i\\mathbf{k}\\cdot\\mathbf{x}}`
of the spectral field `f̂` into `f` (specification, Conventions).

The function leaves `f̂` unchanged on every backend, although a complex-to-real FFT overwrites its
input: on the CPU, the inverse FFT runs in two stages through a scratch array; cuFFT and rocFFT
transform a copy in a buffer of the plan.
The transform uses only the `ky ≥ 0` half of the `kx = 0` column and of the Nyquist column;
[`enforce_hermitian!`](@ref) makes the other half consistent with it.
"""
function to_physical!(f, f̂, tr::Transforms)
    inverse!(f, tr.inverse, f̂)
    return f
end

@kernel function truncate_kernel!(f̂, keep)
    I = @index(Global)
    @inbounds f̂[I] = ifelse(keep[I], f̂[I], zero(eltype(f̂)))
end

"""
    truncate!(f̂, grid::Grid) -> f̂

Set to zero every mode of the spectral field `f̂` outside the Galerkin truncation
`|\\mathbf{k}| < k_c` (`grid.truncation`).
"""
function truncate!(f̂, grid::Grid)
    truncate_kernel!(grid.backend)(f̂, grid.truncation; ndrange = size(f̂))
    return f̂
end

# Product i·k·z for real k, written without a complex multiply.
@inline i_times(k, z) = Complex(-k * imag(z), k * real(z))

@kernel function ∂x_kernel!(df̂, f̂, kx)
    c, r = @index(Global, NTuple)
    @inbounds df̂[c, r] = i_times(kx[c], f̂[c, r])
end

@kernel function ∂y_kernel!(df̂, f̂, ky)
    c, r = @index(Global, NTuple)
    @inbounds df̂[c, r] = i_times(ky[r], f̂[c, r])
end

"""
    ∂x!(df̂, f̂, tr::Transforms) -> df̂

Write the `x` derivative `\\widehat{\\partial_x f} = i k_x^{\\rm d} \\hat f` of the spectral field `f̂`
into `df̂` (specification, Conventions). `df̂` may be `f̂`.
The multiplier is zero on the Nyquist column (see [`Transforms`](@ref)).
"""
function ∂x!(df̂, f̂, tr::Transforms)
    ∂x_kernel!(tr.grid.backend)(df̂, f̂, tr.kx_deriv; ndrange = size(f̂))
    return df̂
end

"""
    ∂y!(df̂, f̂, tr::Transforms) -> df̂

Write the `y` derivative `\\widehat{\\partial_y f} = i k_y^{\\rm d} \\hat f` of the spectral field `f̂`
into `df̂` (specification, Conventions). `df̂` may be `f̂`.
The multiplier is zero on the Nyquist row (see [`Transforms`](@ref)).
"""
function ∂y!(df̂, f̂, tr::Transforms)
    ∂y_kernel!(tr.grid.backend)(df̂, f̂, tr.ky_deriv; ndrange = size(f̂))
    return df̂
end

# Compressive part k (k·u)/k² of the mode (ux, uy) with derivative wavenumbers (kx, ky); zero at k = 0.
@inline function compressive_part(kx, ky, ux, uy)
    k2 = kx^2 + ky^2
    d = ifelse(iszero(k2), zero(ux), (kx * ux + ky * uy) / ifelse(iszero(k2), one(k2), k2))
    return kx * d, ky * d
end

@kernel function helmholtz_kernel!(usx, usy, ucx, ucy, ux, uy, kx, ky)
    c, r = @index(Global, NTuple)
    @inbounds begin
        x, y = ux[c, r], uy[c, r]
        cx, cy = compressive_part(kx[c], ky[r], x, y)
        usx[c, r], usy[c, r] = x - cx, y - cy
        ucx[c, r], ucy[c, r] = cx, cy
    end
end

@kernel function project_kernel!(ux, uy, kx, ky)
    c, r = @index(Global, NTuple)
    @inbounds begin
        x, y = ux[c, r], uy[c, r]
        cx, cy = compressive_part(kx[c], ky[r], x, y)
        ux[c, r], uy[c, r] = x - cx, y - cy
    end
end

"""
    helmholtz_split!(ûs, ûc, û, tr::Transforms) -> (ûs, ûc)

Split the spectral vector field `û = (ûx, ûy)` into its solenoidal part `ûs` and its compressive part
`ûc` (specification, compressible model, Helmholtz split):

```math
\\hat u^s_i = P_{ij}\\hat u_j = \\hat u_i - \\frac{k_i k_j \\hat u_j}{k^2}, \\qquad
\\hat u^c_i = \\frac{k_i k_j \\hat u_j}{k^2}.
```

Each argument is a pair of spectral fields. The function uses the derivative wavenumbers of
[`Transforms`](@ref), so `∂x!` and `∂y!` give an exactly zero divergence of `ûs` and an exactly zero
curl of `ûc`. Where these wavenumbers vanish (the `k = 0` mode and the corner mode
`(Nx/2, Ny/2)`), the whole mode is solenoidal: `ûs = û`, `ûc = 0`. `ûs` or `ûc` may be `û`.
"""
function helmholtz_split!(ûs, ûc, û, tr::Transforms)
    helmholtz_kernel!(tr.grid.backend)(ûs[1], ûs[2], ûc[1], ûc[2], û[1], û[2], tr.kx_deriv, tr.ky_deriv;
                                       ndrange = size(û[1]))
    return ûs, ûc
end

"""
    project_solenoidal!(û, tr::Transforms) -> û

Replace the spectral vector field `û = (ûx, ûy)` by its solenoidal part `P_{ij}\\hat u_j`
(specification, compressible model, Helmholtz split). The projection is the one of
[`helmholtz_split!`](@ref).
"""
function project_solenoidal!(û, tr::Transforms)
    project_kernel!(tr.grid.backend)(û[1], û[2], tr.kx_deriv, tr.ky_deriv; ndrange = size(û[1]))
    return û
end

@kernel function hermitian_kernel!(f̂, Nx, Ny)
    t = @index(Global)
    nc, nr = Nx ÷ 2 + 1, Ny ÷ 2 + 1
    @inbounds begin
        if t <= Ny
            f̂[nc, t] = zero(eltype(f̂))
            if t == 1
                f̂[1, 1] = real(f̂[1, 1])
            elseif t < nr
                f̂[1, Ny + 2 - t] = conj(f̂[1, t])
            end
        end
        if t <= nc
            f̂[t, nr] = zero(eltype(f̂))
        end
    end
end

"""
    enforce_hermitian!(f̂, grid::Grid) -> f̂

Make the spectral field `f̂` the transform of a real field with zero Nyquist modes.
On the `kx = 0` column the function copies the source half `ky > 0` to the other half,
`\\hat f(0, -k_y) = \\hat f(0, k_y)^*`, makes `\\hat f(0, 0)` real, and sets the Nyquist column
`kx = Nx/2` and the Nyquist row `ky = Ny/2` to zero.

The real-FFT layout stores both `(0, k_y)` and `(0, -k_y)`. Every routine that generates spectral
values directly (initial conditions, forcing) calls this function, so the spectral field and its
physical field agree (ADR-0015).
"""
function enforce_hermitian!(f̂, grid::Grid)
    hermitian_kernel!(grid.backend)(f̂, grid.Nx, grid.Ny; ndrange = max(grid.Nx ÷ 2 + 1, grid.Ny))
    return f̂
end
