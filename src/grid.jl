# Periodic rectangular grid, circular Galerkin truncation and spectral shells (WP-1, ADR-0003).

"""
    Grid{T}

Periodic rectangle `Lx × Ly` with `Nx × Ny` points, equal spacing and float type `T`.
Build it with `Grid(; Lx, Ly, Nx, Ny, T = Float64, backend = CPU())`.

The spectral arrays use the real-FFT layout `(Nx ÷ 2 + 1, Ny)`, with `kx` along the first dimension
and the second dimension in FFT order. The fields are:

- `kx`, `ky`, `k2`: wavenumbers `\\hat{k}_x, \\hat{k}_y` and `|\\mathbf{k}|^2` (specification, Conventions).
- `truncation`: `true` for the modes kept by the circular Galerkin truncation `|\\mathbf{k}| < k_c`,
  `k_c = (N_x/3)(2\\pi/L_x)` (specification, Galerkin truncation).
- `shell_index`: `n + 1` for the shell `n` of every grid mode, retained or not.
- `k_max`: the largest retained `|\\mathbf{k}|`.
- `Δk`: the shell width `2\\pi/\\max(L_x, L_y)`.
- `shell_k`, `shell_edge`, `shell_count`, `shell_count_all` (host vectors, entry `n + 1` for shell `n`,
  up to the largest shell of any grid mode): the shell centre `n Δk`, the upper edge `(n + 1/2) Δk`,
  the number of retained modes, and the number of all grid modes. Both counts use the full wavenumber
  plane: weight 1 for the columns `kx = 0` and `kx = Nx/2`, weight 2 elsewhere.

Shell `n` holds the modes with `n - 1/2 ≤ |\\mathbf{k}|/Δk < n + 1/2` (specification, Conventions).
The constructor throws an `ArgumentError` for an odd `Nx` or `Ny`, and for unequal spacing
`|Lx/Nx - Ly/Ny| > 10^{-12} Lx/Nx`.
"""
struct Grid{T<:AbstractFloat,B,A<:AbstractMatrix{T},M<:AbstractMatrix{Bool},I<:AbstractMatrix{Int32}}
    Lx::T
    Ly::T
    Nx::Int
    Ny::Int
    Δk::T
    k_max::T
    backend::B
    kx::A
    ky::A
    k2::A
    truncation::M
    shell_index::I
    shell_k::Vector{T}
    shell_edge::Vector{T}
    shell_count::Vector{Int}
    shell_count_all::Vector{Int}
end

# Signed integer index of row `r` (1-based) in FFT order for `N` points.
signed_index(r, N) = r - 1 < N ÷ 2 ? r - 1 : r - 1 - N

# Smallest shell `n` of the mode `(i, j)`, from integer indices only (ADR-0003, ADR-0016).
# An integer square root gives the start value; the exact integer test corrects it.
function shell_number(i, j, Nx, Ny)
    Nm = Int128(max(Nx, Ny))
    q = 4 * Nm^2 * (Int128(i)^2 * Ny^2 + Int128(j)^2 * Nx^2)
    d = Int128(Nx)^2 * Ny^2
    n = Int(isqrt(q ÷ d)) ÷ 2
    while n > 0 && q < (2n - 1)^2 * d
        n -= 1
    end
    while q >= (2n + 1)^2 * d
        n += 1
    end
    return n
end

# Whether the mode `(i, j)` lies strictly inside the disc `|k| < k_c`, from integer indices only.
retained(i, j, Nx, Ny) = 9 * (Int64(i)^2 * Ny^2 + Int64(j)^2 * Nx^2) < Int64(Nx)^2 * Ny^2

function Grid(; Lx, Ly, Nx::Integer, Ny::Integer, T::Type{<:AbstractFloat} = Float64, backend = CPU())
    (iseven(Nx) && iseven(Ny) && Nx > 0 && Ny > 0) || throw(ArgumentError("Nx and Ny must be positive and even"))
    hx, hy = Float64(Lx) / Nx, Float64(Ly) / Ny
    abs(hx - hy) <= 1e-12 * hx || throw(ArgumentError("grid spacing must be equal: Lx/Nx = $hx, Ly/Ny = $hy"))

    nkx = Nx ÷ 2 + 1
    kx, ky, k2 = (Matrix{T}(undef, nkx, Ny) for _ in 1:3)
    keep = Matrix{Bool}(undef, nkx, Ny)
    shell_index = Matrix{Int32}(undef, nkx, Ny)
    nshell = shell_number(Nx ÷ 2, Ny ÷ 2, Nx, Ny) + 1
    shell_count = zeros(Int, nshell)
    shell_count_all = zeros(Int, nshell)
    k2_max = 0.0
    for r in 1:Ny, c in 1:nkx
        i, j = c - 1, signed_index(r, Ny)
        n = shell_number(i, j, Nx, Ny)
        weight = (i == 0 || i == Nx ÷ 2) ? 1 : 2
        kxf, kyf = 2π * i / Float64(Lx), 2π * j / Float64(Ly)
        kx[c, r], ky[c, r], k2[c, r] = kxf, kyf, kxf^2 + kyf^2
        keep[c, r] = retained(i, j, Nx, Ny)
        shell_index[c, r] = n + 1
        shell_count_all[n + 1] += weight
        if keep[c, r]
            shell_count[n + 1] += weight
            k2_max = max(k2_max, kxf^2 + kyf^2)
        end
    end

    Δk = 2π / max(Float64(Lx), Float64(Ly))
    on_device(a) = (d = KernelAbstractions.allocate(backend, eltype(a), size(a)); copyto!(d, a); d)
    return Grid(T(Lx), T(Ly), Int(Nx), Int(Ny), T(Δk), T(sqrt(k2_max)), backend,
                on_device(kx), on_device(ky), on_device(k2), on_device(keep), on_device(shell_index),
                T.(Δk .* (0:(nshell - 1))), T.(Δk .* ((0:(nshell - 1)) .+ 1 // 2)),
                shell_count, shell_count_all)
end
