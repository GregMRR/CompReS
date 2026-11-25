# r5-E45: the only acceptance item on the compressible injection rates (ROADMAP:177, WP-11) compares them with the
# incompressible ones at ρ = 1. At ρ = 1 every ρ factor is 1, so mutants that drop ρ survive. A forcing-only step
# with nonuniform ρ (ρ held fixed, u_{n+1} = u_n + Δt f) separates them: ΔE_K/Δt and ΔΩ/Δt are exact grid identities.
# Fields are few-mode trigonometric sums, so curl is exact on the grid (no FFT needed). Run: julia r5-E45.jl
using Random, Printf
rng = MersenneTwister(5); N = 16; x = (0:N-1) .* (2π / N); X = [xi for xi in x, yj in x]; Y = [yj for xi in x, yj in x]
modes = [(1, 2), (2, -1), (3, 1), (0, 2), (2, 0)]
function field(rng)   # returns (ax, ay, curl a) of a random few-mode vector field
    ax = zeros(N, N); ay = zeros(N, N); c = zeros(N, N)
    for (p, q) in modes
        A, B, C, D = randn(rng, 4); s = sin.(p .* X .+ q .* Y); co = cos.(p .* X .+ q .* Y)
        ax .+= A .* s .+ B .* co; ay .+= C .* s .+ D .* co
        c .+= p .* (C .* co .- D .* s) .- q .* (A .* co .- B .* s)   # ∂x ay − ∂y ax
    end
    ax, ay, c
end
ux, uy, ω = field(rng); fx, fy, cf = field(rng); Δt = 1e-3
m(a) = sum(a) / length(a)
EK(ρ, ux, uy) = m(0.5 .* ρ .* (ux .^ 2 .+ uy .^ 2)); Om(ρ, ω) = m(0.5 .* ω .^ 2 ./ ρ)
rateE(ρ)     = m(ρ .* (ux .* fx .+ uy .* fy)) + Δt / 2 * m(ρ .* (fx .^ 2 .+ fy .^ 2))      # WP-11 text
rateE_mut(ρ) = m(ux .* fx .+ uy .* fy) + Δt / 2 * m(fx .^ 2 .+ fy .^ 2)                    # mutant: ρ dropped
rateO(ρ)     = m((ω ./ ρ) .* cf) + Δt / 2 * m(cf .^ 2 ./ ρ)                                 # WP-11 text
rateO_mut(ρ) = m(ω .* cf) + Δt / 2 * m(cf .^ 2)                                            # mutant: q -> ω, /ρ dropped
for (name, ρ) in (("rho = 1", ones(N, N)), ("rho = 1 + 0.1 cos(x+y)", 1 .+ 0.1 .* cos.(X .+ Y)))
    dE = (EK(ρ, ux .+ Δt .* fx, uy .+ Δt .* fy) - EK(ρ, ux, uy)) / Δt
    dO = (Om(ρ, ω .+ Δt .* cf) - Om(ρ, ω)) / Δt
    @printf("%-24s E_K: |text - dE/dt| = %.1e  |mutant - dE/dt| = %.1e   Omega: |text - dO/dt| = %.1e  |mutant - dO/dt| = %.1e\n",
            name, abs(rateE(ρ) - dE), abs(rateE_mut(ρ) - dE), abs(rateO(ρ) - dO), abs(rateO_mut(ρ) - dO))
    @printf("%-24s rate scale: E %.2f, Omega %.2f\n", "", abs(dE), abs(dO))
end
