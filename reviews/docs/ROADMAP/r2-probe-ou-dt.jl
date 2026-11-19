# r2-probe-ou-dt: in the planned OU test (ROADMAP:156) the forcing is constant over a step
# (ADR-0005:9) and u' = f. The measured "injection Re[û* f̂]" depends on whether û is taken
# before or after the step. The bias is ∓Δt/(2τ_f). The plan fixes neither Δt/τ_f nor the
# sampling point. Exact discrete values: before = (Δt/τ) a/(1-a), after = (Δt/τ)/(1-a),
# a = exp(-Δt/τ); the true energy input rate Δ(½|u|²)/Δt averages exactly to σ²τ.
# Run: julia r2-probe-ou-dt.jl   (seed 5; Monte Carlo std about 0.016)
using Random
function run(rng, dt_over_tau; M = 200000, nsteps_tau = 25)
    τ = 1.0; Δt = dt_over_tau * τ; a = exp(-Δt / τ); σ2 = 1.0   # complex OU, E|f|² = σ2
    n = round(Int, nsteps_tau / dt_over_tau)
    pre = 0.0; post = 0.0; dE = 0.0
    for _ in 1:M
        f = sqrt(σ2 / 2) * randn(rng, ComplexF64) * sqrt(2); u = 0.0im
        for _ in 1:n-1
            u += Δt * f
            f = a * f + sqrt((1 - a^2) * σ2 / 2) * randn(rng, ComplexF64) * sqrt(2)
        end
        pre += real(conj(u) * f); u2 = u + Δt * f
        post += real(conj(u2) * f); dE += (abs2(u2) - abs2(u)) / (2Δt)
    end
    # normalise by the continuous-time value σ²τ/2 per complex mode with ½|u|² energy convention
    pre / M / σ2, post / M / σ2, dE / M / σ2
end
rng = MersenneTwister(5)
for r in (0.01, 0.1, 0.5)
    a = exp(-r)
    pre, post, dE = run(rng, r)
    println("Δt/τ_f = ", r, ":  before-step ", round(pre, digits = 3), " (exact ", round(r * a / (1 - a), digits = 3),
            "),  after-step ", round(post, digits = 3), " (exact ", round(r / (1 - a), digits = 3),
            "),  energy rate ", round(dE, digits = 3), " (exact 1)")
end
