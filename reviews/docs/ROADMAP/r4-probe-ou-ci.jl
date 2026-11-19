# r4-probe-ou-ci: in the WP-9 OU test (ROADMAP:164) the 21 samples of one realisation (once per τ_f over
# [20τ_f, 40τ_f]) are not independent: u is a random walk (u' = f, no damping). Compares the 95% half-width
# from the per-realisation means (correct) with the one that treats all M·21 samples as independent, and
# counts how often a correct generator fails each interval over R = 96 seeds (M = 1000 per seed; the fail rate does not depend on M).
# Run: julia -t 16 r4-probe-ou-ci.jl
using Random, Statistics, Printf
pairs = [(i, j) for i in -6:6, j in -6:6 if 16 <= i^2 + j^2 <= 36 && (i > 0 || (i == 0 && j > 0))]
P = length(pairs); τ = 1.0; Δt = 0.01τ; a = exp(-Δt / τ); ε = 1.0
σ2 = ε / (2P) * (1 - a) / (Δt / 2 * (1 + a)); b = sqrt((1 - a^2) * σ2 / 2)
n0, n1, ev = 2000, 4000, 100
function realisation!(out, rng)
    f = sqrt(σ2 / 2) .* complex.(randn(rng, P), randn(rng, P)); u = zeros(ComplexF64, P); m = 0
    for step in 0:n1
        if step >= n0 && (step - n0) % ev == 0
            m += 1; out[m] = 2 * sum(real.(conj.(u) .* f)) + Δt * sum(abs2, f)
        end
        u .+= Δt .* f; f .= a .* f .+ b .* complex.(randn(rng, P), randn(rng, P))
    end
end
M = 1000; R = 96; failC = Threads.Atomic{Int}(0); failN = Threads.Atomic{Int}(0); ratio = zeros(R)
Threads.@threads for r in 1:R
    rng = Xoshiro(10_000 + r); X = zeros(21, M)
    for i in 1:M; realisation!(view(X, :, i), rng); end
    mC = mean(X); hwC = 1.96 * std(vec(mean(X, dims = 1))) / sqrt(M); hwN = 1.96 * std(vec(X)) / sqrt(21M)
    abs(mC - ε) > hwC && Threads.atomic_add!(failC, 1); abs(mC - ε) > hwN && Threads.atomic_add!(failN, 1)
    ratio[r] = hwC / hwN
end
@printf("half-width ratio (per-realisation / naive) = %.2f\n", mean(ratio))
@printf("correct generator fails: per-realisation interval %d/%d, naive interval %d/%d\n", failC[], R, failN[], R)
