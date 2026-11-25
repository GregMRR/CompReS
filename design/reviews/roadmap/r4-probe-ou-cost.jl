# r4-probe-ou-cost: lower bound on the CPU time of the WP-9 OU normalisation test (ROADMAP:164):
# M = 4000 realisations × 4000 steps (t ≤ 40τ_f, Δt = 0.01τ_f) = 1.6e7 calls of update_forcing!.
# Per step, the minimum work the plan prescribes: reseed a generator by (seed, step) (ADR-0007:17),
# draw the OU noise for the forced modes, apply the exact OU update. The forcing state is held per forced
# mode (index list, ADR-0007:5), so no FFT, no projection, no enforce_hermitian!, no u update are counted.
# Run: julia r4-probe-ou-cost.jl   (single thread)
using Random
modes = [(i, j) for i in 0:6, j in -6:6 if 16 <= i^2 + j^2 <= 36]     # band(4,6), half plane, N = 32
P = length(modes)
f = zeros(ComplexF64, 2, P); a = exp(-0.01); b = sqrt(1 - a^2)
function update!(f, seed, step, a, b)
    rng = Xoshiro(seed, step, 0x9e3779b97f4a7c15, 0x2545f4914f6cdd1d)   # raw state words: cheapest seeding
    @inbounds for k in axes(f, 2), c in 1:2
        f[c, k] = a * f[c, k] + b * complex(randn(rng), randn(rng))
    end
end
function update_sha!(f, seed, step, a, b)
    rng = Xoshiro(hash((seed, step)))            # integer seed: Julia hashes it with SHA-256
    @inbounds for k in axes(f, 2), c in 1:2
        f[c, k] = a * f[c, k] + b * complex(randn(rng), randn(rng))
    end
end
for (name, g) in (("Xoshiro(raw state words)", update!), ("Xoshiro(hash(seed, step))", update_sha!))
    g(f, 1, 1, a, b); n = 200_000
    t = @elapsed for s in 1:n; g(f, 1, s, a, b); end
    println(name, ": ", round(t / n * 1e9, digits = 0), " ns/step  -> 1.6e7 steps: ",
            round(t / n * 1.6e7, digits = 1), " s   (P = $P forced modes, 2 components)")
end
# Full-layout variant: scatter into f̂ (17×32×2), project, enforce_hermitian! (copy kx = 0, zero Nyquist).
N = 32; F = zeros(ComplexF64, N ÷ 2 + 1, N, 2)
kxv = [0:N÷2;]; kyv = [0:N÷2-1; -N÷2:-1]
idx = [(i + 1, j >= 0 ? j + 1 : N + j + 1) for (i, j) in modes]
function full!(F, f, idx, kxv, kyv, N)
    fill!(F, 0)
    @inbounds for (k, (a, b)) in enumerate(idx)
        kx = kxv[a]; ky = kyv[b]; k2 = kx^2 + ky^2
        fx, fy = f[1, k], f[2, k]; d = (kx * fx + ky * fy) / k2
        F[a, b, 1] = fx - kx * d; F[a, b, 2] = fy - ky * d
    end
    @inbounds for c in 1:2
        for j in 2:N÷2; F[1, N + 2 - j, c] = conj(F[1, j, c]); end
        F[1, 1, c] = real(F[1, 1, c]); F[end, :, c] .= 0; F[:, N÷2+1, c] .= 0
    end
end
full!(F, f, idx, kxv, kyv, N); n = 200_000
t = @elapsed for s in 1:n; update!(f, 1, s, a, b); full!(F, f, idx, kxv, kyv, N); end
println("update + full-layout scatter/projection/enforce_hermitian!: ", round(t / n * 1e9), " ns/step -> 1.6e7 steps: ",
        round(t / n * 1.6e7, digits = 1), " s")
