# r3-probe-hermitian: ROADMAP:112-113 defines enforce_hermitian! on the kx = 0 column only
# (f̂(0,−ky) = conj f̂(0,ky), f̂(0,0) real) and requires its output to "survive an irfft/rfft round
# trip unchanged". Input: random complex array in the real-FFT layout (N = 16, seed 1),
# (i) untruncated, (ii) truncated to the disc. Reports the max change and where it sits.
# Run: julia --project=<env with FFTW> r3-probe-hermitian.jl
using FFTW, Random
N = 16; rng = Xoshiro(1)
ky = [0:N÷2-1; -N÷2:-1]; kx = 0:N÷2
function enforce_spec!(f)                      # exactly as specified in ROADMAP:112
    for j in 2:N, jj in 2:N
        ky[jj] == -ky[j] && ky[j] > 0 && (f[1, jj] = conj(f[1, j]))
    end
    f[1, 1] = real(f[1, 1]); f
end
for (name, trunc) in (("untruncated", false), ("truncated  ", true))
    f = randn(rng, ComplexF64, N ÷ 2 + 1, N)
    trunc && (f .*= [9 * (i^2 + j^2) < N^2 for i in kx, j in ky])
    enforce_spec!(f)
    g = rfft(irfft(f, N))
    d = abs.(g .- f); i = argmax(d)
    println(name, ": max change ", round(maximum(d), sigdigits = 3), " at (kx, ky) = (", kx[i[1]], ", ", ky[i[2]], ")",
            ";  changed entries: ", count(>(1e-12), d), "  (kx=Nx/2 column: ", count(>(1e-12), d[end, :]),
            ", (0, -N/2): ", d[1, N÷2+1] > 1e-12, ")")
end
