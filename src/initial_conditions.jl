# Initial conditions as spectral fields that a State takes (WP-4, ADR-0002).

# Device copy of the host array `a` on `backend`.
to_device(backend, a) = copyto!(KernelAbstractions.allocate(backend, eltype(a), size(a)), a)

# Host array of zeros for one spectral field of `grid`.
host_spectral_zeros(grid::Grid{T}) where {T} = zeros(Complex{T}, grid.Nx ÷ 2 + 1, grid.Ny)

# Array position of the mode with integer indices `(i, j)`, `i ≥ 0`, in the real-FFT layout.
# Throws an `ArgumentError` when the mode lies outside the Galerkin truncation.
function retained_position(grid::Grid, i, j)
    retained(i, j, grid.Nx, grid.Ny) ||
        throw(ArgumentError("mode ($i, $j) lies outside the Galerkin truncation of the grid"))
    return (i + 1, j < 0 ? grid.Ny + j + 1 : j + 1)
end

# Device spectral field of `host`, with the Hermitian symmetry of a real field.
device_field(grid, host) = enforce_hermitian!(to_device(grid.backend, host), grid)

"""
    kolmogorov(grid::Grid{T}; A, n) -> (; ux, uy)

The Kolmogorov flow `\\mathbf{u} = A \\sin(k_n y)\\, \\hat{\\mathbf{x}}`, `k_n = n\\, 2\\pi/L_y`, as spectral
fields of float type `T` (specification, Conventions):

```math
\\hat u_x(0, \\pm k_n) = \\mp \\frac{iA}{2}, \\qquad \\hat u_y = 0.
```

`n` is the mode index along `y`, the same index as in [`crossed_kolmogorov`](@ref).
The function throws an `ArgumentError` for `n < 1` and for a mode outside the Galerkin truncation.
"""
function kolmogorov(grid::Grid{T}; A, n::Integer) where {T}
    n >= 1 || throw(ArgumentError("mode index n must be at least 1, not $n"))
    ux = host_spectral_zeros(grid)
    ux[retained_position(grid, 0, n)...] = -im * T(A) / 2
    return (ux = device_field(grid, ux), uy = device_field(grid, host_spectral_zeros(grid)))
end

"""
    crossed_kolmogorov(grid::Grid{T}; A1, k1, A2, k2) -> (; ux, uy)

The crossed Kolmogorov flow of the user acceptance test (ADR-0011),

```math
\\mathbf{u} = A_1 \\sin(k_1 y)\\, \\hat{\\mathbf{x}} + A_2 \\sin(k_2 x)\\, \\hat{\\mathbf{y}},
```

as spectral fields of float type `T`. `k1` and `k2` are integer mode indices, as `n` in
[`kolmogorov`](@ref): the wavenumbers are `k1\\, 2\\pi/L_y` and `k2\\, 2\\pi/L_x`.
The function throws an `ArgumentError` for an index below 1 and for a mode outside the Galerkin
truncation.
"""
function crossed_kolmogorov(grid::Grid{T}; A1, k1::Integer, A2, k2::Integer) where {T}
    (k1 >= 1 && k2 >= 1) || throw(ArgumentError("mode indices k1, k2 must be at least 1, not $k1, $k2"))
    ux, uy = host_spectral_zeros(grid), host_spectral_zeros(grid)
    ux[retained_position(grid, 0, k1)...] = -im * T(A1) / 2
    uy[retained_position(grid, k2, 0)...] = -im * T(A2) / 2
    return (ux = device_field(grid, ux), uy = device_field(grid, uy))
end

"""
    uniform_density(grid::Grid{T}) -> (; rho)

The uniform density `\\rho = 1` as a spectral field of float type `T`: `\\hat\\rho(\\mathbf{0}) = 1`
and every other mode is zero. Merge it with a velocity field to get the fields of a compressible
[`State`](@ref), for example `merge(kolmogorov(grid; A, n), uniform_density(grid))`.
"""
function uniform_density(grid::Grid{T}) where {T}
    rho = host_spectral_zeros(grid)
    rho[1, 1] = one(T)
    return (rho = to_device(grid.backend, rho),)
end

"""
    random_velocity(grid::Grid{T}; spectrum, E_kin, seed) -> (; ux, uy)

A random solenoidal velocity field with shell spectrum proportional to `spectrum(k)` and kinetic
energy `E_kin`, as spectral fields of float type `T`.

Each retained mode `\\mathbf{q} \\neq 0` of shell `n` gets

```math
\\hat{\\mathbf{u}}(\\mathbf{q}) = a_n\\, e^{i\\varphi(\\mathbf{q})}\\,
\\frac{(-q_y, q_x)}{|\\mathbf{q}|}, \\qquad
a_n^2 = c\\, \\frac{2\\, \\mathrm{spectrum}(k_n)}{M_n},
```

with `k_n` the shell centre, `M_n` the number of retained modes in shell `n` over the full plane
(`grid.shell_count`) and a phase `\\varphi` uniform in `[0, 2\\pi)`. The shell spectrum
`E(k_n) = \\tfrac12 \\sum_{\\mathbf{q} \\in k_n} |\\hat{\\mathbf{u}}(\\mathbf{q})|^2` is then
`c\\,\\mathrm{spectrum}(k_n)` (specification, incompressible model, spectra), and the factor `c`
sets `E_{\\rm kin} = \\tfrac12 \\langle |\\mathbf{u}|^2 \\rangle = \\sum_n E(k_n)`.
The field is solenoidal, `\\mathbf{q} \\cdot \\hat{\\mathbf{u}} = 0`, lies inside the Galerkin
truncation and has no mean flow. [`enforce_hermitian!`](@ref) sets the half `ky < 0` of the
`kx = 0` column.

A `Random.Xoshiro(seed)` generator on the host draws the phases, and the host computes the field
before the copy to the device, so one seed gives the same field on every backend, bit for bit.
The function throws an `ArgumentError` when `spectrum` is zero on every retained mode.
"""
function random_velocity(grid::Grid{T}; spectrum, E_kin, seed) where {T}
    kx, ky, keep, shell = Array(grid.kx), Array(grid.ky), Array(grid.truncation), Array(grid.shell_index)
    phase = rand(Xoshiro(seed), T, size(kx))
    ux, uy = host_spectral_zeros(grid), host_spectral_zeros(grid)
    energy = zero(T)
    for I in CartesianIndices(ux)
        k = sqrt(kx[I]^2 + ky[I]^2)
        (keep[I] && k > 0) || continue
        n = shell[I]
        a = sqrt(2 * T(spectrum(grid.shell_k[n])) / grid.shell_count[n]) * cis(2 * T(π) * phase[I])
        ux[I], uy[I] = -a * ky[I] / k, a * kx[I] / k
        # Full-plane weight of the stored mode: 1 on the kx = 0 column, 2 elsewhere (no Nyquist mode is retained).
        energy += (I[1] == 1 ? 1 : 2) * (abs2(ux[I]) + abs2(uy[I])) / 2
    end
    energy > 0 || throw(ArgumentError("spectrum is zero on every retained mode"))
    c = sqrt(T(E_kin) / energy)
    return (ux = device_field(grid, ux .* c), uy = device_field(grid, uy .* c))
end
