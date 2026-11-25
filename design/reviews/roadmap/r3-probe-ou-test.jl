# r3-probe-ou-test: power and sampling-point sensitivity of the planned OU normalisation test
# (ROADMAP:162): N = 32, L = 2π, band(4, 6), Δt = 0.01τ_f, M = 4000, u' = f from u = 0,
# injection = ⟨u·f⟩ + ½Δt⟨|f|²⟩ averaged over t ∈ [20τ_f, 40τ_f], pass if inside the 95% interval around ε.
# A Hermitian field: physical ⟨u·f⟩ = Σ_{all modes} Re[û*f̂] = 2 Σ_{pairs ±q} Re[û*f̂]. One complex OU per pair.
# Reports the ensemble mean / ε and its 95% half-width for:
#   (a) û sampled at the start of the step (before u += Δt f);  (b) û sampled after the step;
#   (c) mutant: kx = 0 column not Hermitian (each of its pairs injects half, as in r2-E20).
# Run: julia -t 8 r3-probe-ou-test.jl   (seed 3)
using Random, Statistics, Printf
pairs = [(i, j) for i in -6:6, j in -6:6 if 16 <= i^2 + j^2 <= 36 && (i > 0 || (i == 0 && j > 0))]
P = length(pairs); P0 = count(p -> p[1] == 0, pairs)
τ = 1.0; Δt = 0.01τ; a = exp(-Δt / τ); ε = 1.0
σ2 = ε / (2P) * (1 - a) / (Δt / 2 * (1 + a))   # per-pair E|f|², exact discrete rate = ε
nwin0 = round(Int, 20τ / Δt); nend = round(Int, 40τ / Δt); M = 4000
w = [p[1] == 0 ? 1.0 : 1.0 for p in pairs]    # weight per pair (both ±q are physical modes)
function realisation(seed)
    rng = Xoshiro(seed); f = sqrt(σ2 / 2) .* complex.(randn(rng, P), randn(rng, P)); u = zeros(ComplexF64, P)
    sa = 0.0; sb = 0.0; sc = 0.0; n = 0; b = sqrt((1 - a^2) * σ2 / 2)
    for step in 0:nend-1
        if step >= nwin0
            ff = sum(abs2, f)
            ia = 2 * sum(real.(conj.(u) .* f)); u2 = u .+ Δt .* f; ib = 2 * sum(real.(conj.(u2) .* f))
            # mutant: kx = 0 pairs inject half (independent ±ky draws, Hermitian part keeps half)
            ic = ia - sum(real(conj(u[k]) * f[k]) for k in 1:P if pairs[k][1] == 0)
            fc = ff - 0.5 * sum(abs2(f[k]) for k in 1:P if pairs[k][1] == 0)
            sa += ia + Δt * ff; sb += ib + Δt * ff; sc += ic + Δt * fc; n += 1
        end
        u .+= Δt .* f
        f .= a .* f .+ b .* complex.(randn(rng, P), randn(rng, P))
    end
    (sa / n, sb / n, sc / n)
end
res = Vector{NTuple{3,Float64}}(undef, M)
Threads.@threads for r in 1:M; res[r] = realisation(3_000_000 + r); end
@printf("P = %d pairs (%d physical modes), %d pairs on kx = 0\n", P, 2P, P0)
for (k, name) in enumerate(("(a) start of step", "(b) after step   ", "(c) mutant non-Hermitian kx=0"))
    v = [r[k] for r in res]; m = mean(v); hw = 1.96 * std(v) / sqrt(M)
    @printf("%s: mean/ε = %.4f ± %.4f (95%%)  -> %s\n", name, m / ε, hw / ε, abs(m - ε) <= hw ? "pass" : "FAIL")
end
