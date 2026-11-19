# r5-probe-injtol: the WP-7 injection test (ROADMAP:147) compares inj_E_kin with (E_kin(t_{n+1}) − E_kin(t_n))/Δt
# "to round-off". The difference of two energies cancels: the relative error grows like eps(T)·E/(Δt·inj) ≈ eps(T)·n.
# This probe measures that error in Float64 and Float32 (spectral sums over half-plane modes, weight 2, as in
# ROADMAP:154), relative to |dE/dt| (which crosses 0) and to the scale eps(T)·E/Δt, and the errors of three mutants of record_injection!. Seed 7. Run: julia r5-probe-injtol.jl
using Random, Printf
function run(T, nsteps; mut = :none, seed = 7)
    rng = MersenneTwister(seed); P = 40; τ = T(1); Δt = T(0.01); a = exp(-Δt / τ); ε = T(1)
    σ2 = ε / (2P) / (τ * (1 + T(0) )) ; b = sqrt((1 - a^2) * σ2 / 2)   # mean u·f ≈ ε for t ≫ τ
    f = sqrt(σ2 / 2) .* complex.(randn(rng, T, P), randn(rng, T, P)); fold = copy(f); u = zeros(Complex{T}, P)
    E(u) = sum(abs2, u)                                        # ½ Σ w |û|², w = 2
    worst = zero(T); worstS = zero(T)
    for n in 1:nsteps
        fold .= f; f .= a .* f .+ b .* complex.(randn(rng, T, P), randn(rng, T, P))   # update_forcing!
        ff = mut == :fold ? fold : f
        inj = 2 * sum(real.(conj.(u) .* ff)) + Δt * sum(abs2, ff)   # Σ w (Re[û*f̂] + ½Δt|f̂|²)
        mut == :nohalf && (inj = 2 * sum(real.(conj.(u) .* ff)))
        E0 = E(u); u .+= Δt .* f; E1 = E(u)
        mut == :post && (inj = 2 * sum(real.(conj.(u) .* ff)) + Δt * sum(abs2, ff))
        n > 10 && (worst = max(worst, abs(inj - (E1 - E0) / Δt) / abs((E1 - E0) / Δt)))
        n > 10 && (worstS = max(worstS, abs(inj - (E1 - E0) / Δt) / (eps(T) * E1 / Δt)))
    end
    worst, worstS
end
for T in (Float64, Float32), n in (100, 1000, 4000)
    r, sc = run(T, n); @printf("%-8s n = %4d  correct: max |err|/|dE/dt| = %.1e (= %.0f eps);  max |err|/(eps E/Δt) = %.2f\n", T, n, r, r / eps(T), sc)
end
for mut in (:post, :nohalf, :fold)
    r, sc = run(Float32, 1000; mut); @printf("mutant %-7s Float32 n = 1000: max |err|/|dE/dt| = %.1e;  max |err|/(eps E/Δt) = %.0f\n", mut, r, sc)
end
