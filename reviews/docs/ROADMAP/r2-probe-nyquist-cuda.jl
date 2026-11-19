# r2-probe-nyquist-cuda: same mutant as r2-probe-nyquist.jl, with cuFFT (needs a functional CUDA).
# Run: julia --project=<env with CUDA> r2-probe-nyquist-cuda.jl
# Output on the development machine (RTX 5070): test 1 zeroed = 0.0, mutant = 0.0;
# test 2 zeroed = mutant = −2.13e-14. Both planned tests pass on the mutant with cuFFT too.
using CUDA, CUDA.CUFFT, Random
N = 16; x = (0:N-1) .* (2π/N); X = [xi for xi in x, yj in x]
kx = Float64.([0:N÷2;]); kxm = copy(kx); kxm[end] = -N/2; kx[end] = 0
dx(f, k) = Array(irfft(im .* CuArray(k) .* rfft(CuArray(f)), N))
f = cos.(N/2 .* X)
println("cuFFT test 1: zeroed = ", maximum(abs, dx(f, kx)), "  mutant = ", maximum(abs, dx(f, kxm)))
Random.seed!(3); g = randn(N, N); h = randn(N, N)
for (n, k) in (("zeroed", kx), ("mutant", kxm))
    println("cuFFT test 2 (", n, "): ", sum(g .* dx(h, k)) + sum(h .* dx(g, k)))
end
