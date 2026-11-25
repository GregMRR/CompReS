# r4-probe-hermitian: amended enforce_hermitian! (ROADMAP:114): copy kx = 0 column from ky > 0 to −ky
# (conj), make (0,0) real, zero the Nyquist column (kx = Nx/2) and row (ky = −Ny/2).
# Test (ROADMAP:115): output from random untruncated input survives irfft/rfft unchanged to round-off.
# Rectangular and square grids, seed 1. Run: julia --project=<env with FFTW> r4-probe-hermitian.jl
using FFTW, Random
function enforce!(f, Nx, Ny)
    ky = [0:Ny÷2-1; -Ny÷2:-1]
    for j in 2:Ny÷2                     # source ky = 1 … Ny/2−1 (index j), target −ky at index Ny+2−j
        f[1, Ny + 2 - j] = conj(f[1, j])
    end
    f[1, 1] = real(f[1, 1]); f[end, :] .= 0; f[:, Ny÷2+1] .= 0; f
end
rng = Xoshiro(1)
for (Nx, Ny) in ((16, 16), (8, 24), (32, 8))
    f = enforce!(randn(rng, ComplexF64, Nx ÷ 2 + 1, Ny), Nx, Ny)
    g = rfft(irfft(f, Nx))
    println("Nx=$Nx Ny=$Ny: max |rfft(irfft(f)) − f| = ", maximum(abs.(g .- f)))
end
