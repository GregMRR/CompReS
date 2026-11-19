# r1-E3: the WP-9 OU acceptance "mean injection within 5% of ε ... with white-noise velocity".
# One real OU mode f with correlation time τ = 0.05 and variance σ² = ε/τ (ε = 1), exact OU update.
# (a) velocity = independent white noise: ⟨u f⟩ → 0, not ε.
# (b) velocity driven by the forcing, u' = -γu + f: ⟨u f⟩ = σ²τ/(1+γτ); it depends on the
#     response rate γ of the forced mode, so "within 5% of ε" holds only if γτ < 0.053.
# Ensemble of M realisations, seed 7. Run: julia r1-E3.jl
using Random, Statistics
rng = MersenneTwister(7)
ε = 1.0; τ = 0.05; σ2 = ε / τ; dt = 5e-4; M = 100_000
a = exp(-dt / τ); b = sqrt(σ2 * (1 - a^2))
f = sqrt(σ2) .* randn(rng, M)
println("(a) white-noise velocity: <u f> = ", round(mean(randn(rng, M) .* f); digits=4), "  (target ε = 1)")
for γ in (1.0, 10.0)
    local f = sqrt(σ2) .* randn(rng, M); local u = zeros(M)
    e = exp(-γ * dt); nsteps = round(Int, 8 / γ / dt + 10τ / dt)
    for n in 1:nsteps
        u .= e .* u .+ f .* (1 - e) / γ          # exact for f constant in the step
        f .= a .* f .+ b .* randn(rng, M)
    end
    # sample <u f> as an ensemble mean at the final time (f decorrelated from u only through the OU memory)
    est = mean(u .* f); se = std(u .* f) / sqrt(M)
    println("(b) γ = $γ: <u f> = ", round(est; digits=3), " ± ", round(2se; digits=3),
            " (95%)   theory σ²τ/(1+γτ) = ", round(σ2 * τ / (1 + γ * τ); digits=3))
end
