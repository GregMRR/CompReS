# r6-probe-injtol-first-step: the WP-7 tolerance 10·eps(T)·E_kin/Δt (ROADMAP:147) does not say E_kin at t_n or t_{n+1}.
# From rest (u = 0), E_kin(t_n) = 0 at step 1, so the tolerance is 0 and a correct implementation fails ("Inf").
# With E_kin(t_{n+1}) every scheme form stays below 4.2 eps. Half-plane weight 2, OU forcing, seed 7. Run: julia r6-probe-injtol-first-step.jl
using Random, Printf
# forcing-only step u_{n+1} = scheme(u_n, f_n), OU forcing, half-plane weight 2; inj = Σw(Re[u*f] + ½Δt|f|²)
function stepu!(u, f, Δt, scheme)
    T = real(eltype(u))
    if scheme == :euler
        u .+= Δt .* f
    elseif scheme == :rk4seq        # u += Δt b_i k_i, stage by stage
        for b in (T(1)/6, T(1)/3, T(1)/3, T(1)/6); u .+= (Δt * b) .* f; end
    elseif scheme == :ssprk3        # Shu–Osher form
        u1 = u .+ Δt .* f
        u2 = T(3)/4 .* u .+ T(1)/4 .* (u1 .+ Δt .* f)
        u .= T(1)/3 .* u .+ T(2)/3 .* (u2 .+ Δt .* f)
    elseif scheme == :rk2heun
        u1 = u .+ Δt .* f
        u .= u .+ (Δt/2) .* (f .+ f)
    end
end
function run(T, nsteps, P, scheme; seed = 7)
    rng = MersenneTwister(seed); τ = T(1); Δt = T(0.01); a = exp(-Δt/τ); σ2 = T(1)/(2P); b = sqrt((1-a^2)*σ2/2)
    f = sqrt(σ2/2) .* complex.(randn(rng,T,P), randn(rng,T,P)); u = zeros(Complex{T}, P)
    E(u) = sum(abs2, u)
    wn = zero(T); wn1 = zero(T); first = NaN; wn_after10 = zero(T)
    for n in 1:nsteps
        f .= a .* f .+ b .* complex.(randn(rng,T,P), randn(rng,T,P))
        inj = 2*sum(real.(conj.(u) .* f)) + Δt*sum(abs2, f)
        E0 = E(u); stepu!(u, f, Δt, scheme); E1 = E(u)
        err = abs(inj - (E1-E0)/Δt)
        sn = E0 == 0 ? (err == 0 ? 0 : Inf) : err/(eps(T)*E0/Δt)
        sn1 = err/(eps(T)*E1/Δt)
        n == 1 && (first = sn)
        wn = max(wn, sn); wn1 = max(wn1, sn1); n > 10 && (wn_after10 = max(wn_after10, sn))
    end
    first, wn, wn_after10, wn1
end
for T in (Float64, Float32), P in (2, 31), s in (:euler, :rk2heun, :ssprk3, :rk4seq)
    r = run(T, 1000, P, s)
    @printf("%-8s P=%2d %-8s  step1 err/(eps E(t_n)/Δt)=%-6.3g  max over n of err/(eps E(t_n)/Δt)=%-8.3g (n>10: %-6.3g)  max err/(eps E(t_{n+1})/Δt)=%.3g\n", T, P, s, r...)
end
