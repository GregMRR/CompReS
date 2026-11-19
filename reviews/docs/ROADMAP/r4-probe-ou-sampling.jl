# r4-probe-ou-sampling: the amended WP-9 OU test (ROADMAP:164): sample ⟨u·f⟩ + ½Δt⟨|f|²⟩ at the start of a
# step once per τ_f over t ∈ [20τ_f, 40τ_f] (21 samples), M = 4000, 95% interval from the spread of the
# per-realisation means. Compares with sampling at every step (r3 set-up), and checks the claim
# "samples closer than τ_f add no information" via the ratio of the two standard errors.
# Mutants: (c) kx = 0 column non-Hermitian (r2-E20), (d) OU variance ×2, (e) ×½.
# Model as r3-probe-ou-test.jl: one complex OU per ±q pair, Hermitian field, L = 2π, N = 32, band(4, 6).
# Run: julia -t 16 r4-probe-ou-sampling.jl   (seed 3)
using Random, Statistics, Printf
pairs = [(i, j) for i in -6:6, j in -6:6 if 16 <= i^2 + j^2 <= 36 && (i > 0 || (i == 0 && j > 0))]
P = length(pairs); on0 = [p[1] == 0 for p in pairs]
τ = 1.0; Δt = 0.01τ; a = exp(-Δt / τ); ε = 1.0
σ2 = ε / (2P) * (1 - a) / (Δt / 2 * (1 + a))
n0 = round(Int, 20τ / Δt); n1 = round(Int, 40τ / Δt); every = round(Int, τ / Δt)
function realisation(seed, vf)
    rng = Xoshiro(seed); f = sqrt(vf * σ2 / 2) .* complex.(randn(rng, P), randn(rng, P)); u = zeros(ComplexF64, P)
    b = sqrt((1 - a^2) * vf * σ2 / 2)
    sA = 0.0; nA = 0; sT = 0.0; nT = 0; sC = 0.0; nC = 0
    for step in 0:n1
        if step >= n0
            ia = 2 * sum(real.(conj.(u) .* f)) + Δt * sum(abs2, f)
            ic = ia - sum(real(conj(u[k]) * f[k]) for k in 1:P if on0[k]) - 0.5Δt * sum(abs2(f[k]) for k in 1:P if on0[k])
            sA += ia; nA += 1
            if (step - n0) % every == 0; sT += ia; nT += 1; sC += ic; nC += 1; end
        end
        u .+= Δt .* f
        f .= a .* f .+ b .* complex.(randn(rng, P), randn(rng, P))
    end
    (sA / nA, sT / nT, sC / nC, nT)
end
M = 4000
for (vf, lab) in ((1.0, "correct"), (2.0, "variance ×2"), (0.5, "variance ×½"))
    res = Vector{NTuple{4,Float64}}(undef, M)
    Threads.@threads for r in 1:M; res[r] = realisation(3_000_000 + r, vf); end
    for (k, name) in ((1, "every step       "), (2, "once per τ_f     "), (3, "once per τ_f, non-Hermitian kx=0"))
        (vf != 1.0 && k != 2) && continue
        v = [r[k] for r in res]; m = mean(v); hw = 1.96 * std(v) / sqrt(M)
        @printf("%-12s %s: mean/ε = %.4f ± %.4f -> %s\n", lab, name, m, hw, abs(m - ε) <= hw ? "pass" : "FAIL")
    end
    vf == 1.0 && @printf("samples per realisation once per τ_f: %d; SE ratio (τ_f / every step) = %.2f\n",
                         Int(res[1][4]), std([r[2] for r in res]) / std([r[1] for r in res]))
end
