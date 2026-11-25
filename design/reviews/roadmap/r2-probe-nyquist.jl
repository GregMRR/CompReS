# r2-probe-nyquist: do the two planned Nyquist tests of WP-2 (ROADMAP:111) fail when the
# derivative multiplier is NOT zeroed on the Nyquist column?
# Mutant: ∂x f = irfft(i*kx .* rfft(f)) with kx(N/2) = -N/2 (FFTW frequency layout), not zeroed.
# Run: julia --project=<env with FFTW> r2-probe-nyquist.jl
using FFTW, Random
N = 16; x = (0:N-1) .* (2π/N); X = [xi for xi in x, yj in x]
kx = [0:N÷2;]                       # rfft half axis, Nyquist at index N÷2+1
kx_mut = Float64.(kx); kx_mut[end] = -N/2   # mutant: Nyquist kept
kx_ok  = Float64.(kx); kx_ok[end]  = 0.0    # plan: Nyquist zeroed
dx(f, k) = irfft(im .* k .* rfft(f), N)
spec_dx(f, k) = im .* k .* rfft(f)
f = cos.(N/2 .* X)
println("test 1, d/dx cos(N x/2), physical max |.|: zeroed = ", maximum(abs, dx(f, kx_ok)),
        "  mutant = ", maximum(abs, dx(f, kx_mut)))
println("         spectral max |.|:          zeroed = ", maximum(abs, spec_dx(f, kx_ok)),
        "  mutant = ", maximum(abs, spec_dx(f, kx_mut)))
Random.seed!(3); g = randn(N, N); h = randn(N, N)
for (name, k) in (("zeroed", kx_ok), ("mutant", kx_mut))
    s = sum(g .* dx(h, k)) + sum(h .* dx(g, k))
    println("test 2, Σ g∂h + Σ h∂g (", name, ") = ", s)
end
# Consequence in spectral space: the mutant output on the Nyquist column is not Hermitian.
# A spectral diagnostic (Parseval with half-plane weights) then disagrees with the grid mean.
w = fill(2.0, N÷2+1); w[1] = 1; w[end] = 1
for (name, k) in (("zeroed", kx_ok), ("mutant", kx_mut))
    d̂ = spec_dx(h, k) ./ N^2
    println("Parseval of (∂x h)² (", name, "): spectral = ", sum(w .* abs2.(d̂)),
            "  grid mean = ", sum(abs2, dx(h, k)) / N^2)
end
