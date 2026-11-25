# r2-E20: the planned OU normalisation test (ROADMAP:156) sums stored half-plane
# modes with weights w_k. A forcing generator that draws independent values for the stored
# pairs (0, ky) and (0, -ky) on the kx = 0 column is not Hermitian. The spectral injection
# then differs from the injection into the physical field, and the planned test cannot see it.
# Setup: N = 32, L = 2π, band(4, 6), real-FFT layout (N÷2+1, N), u' = f, f constant over Δt,
# the u and f entries drawn as in a stationary state with <u f> = per-mode target.
# Run: julia --project=<env with FFTW> r2-probe-hermitian.jl
using FFTW, Random
N = 32
kxs = 0:N÷2; kys = [0:N÷2-1; -N÷2:-1]
supp = [(i, j) for i in eachindex(kxs), j in eachindex(kys) if 4 <= hypot(kxs[i], kys[j]) <= 6]
w(i) = (kxs[i] == 0 || kxs[i] == N÷2) ? 1.0 : 2.0
nkx0 = count(p -> kxs[p[1]] == 0, supp)
println("stored support modes: ", length(supp), "  of which on kx = 0: ", nkx0)
function trial(rng; hermitian)
    û = zeros(ComplexF64, N÷2+1, N); f̂ = zeros(ComplexF64, N÷2+1, N)
    for (i, j) in supp
        z = randn(rng, ComplexF64); e = randn(rng, ComplexF64)
        f̂[i, j] = z; û[i, j] = z + e        # <Re(û* f̂)> = 1 per stored mode
    end
    if hermitian                             # enforce f̂(0,-ky) = conj f̂(0,ky)
        for (i, j) in supp
            kxs[i] == 0 && kys[j] < 0 || continue
            jp = findfirst(==(-kys[j]), kys)
            f̂[i, j] = conj(f̂[i, jp]); û[i, j] = conj(û[i, jp])
        end
    end
    spec = sum(w(i) * real(conj(û[i, j]) * f̂[i, j]) for (i, j) in supp)
    u = irfft(û .* N^2, N); f = irfft(f̂ .* N^2, N)
    phys = sum(u .* f) / N^2
    spec, phys
end
for herm in (true, false)
    rng = MersenneTwister(11); M = 4000
    r = [trial(rng; hermitian = herm) for _ in 1:M]
    s = sum(first, r) / M; p = sum(last, r) / M
    println(herm ? "Hermitian generator    " : "non-Hermitian generator",
            ": mean spectral injection = ", round(s, digits = 2),
            "   mean physical <u f> = ", round(p, digits = 2), "   ratio = ", round(p / s, digits = 3))
end
