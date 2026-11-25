# r6-E50: power of the WP-11 forcing-only test (ROADMAP:178) against formulas that drop ρ, by forcing support.
# u starts from rest, so u stays in the span of the support. ρ = 1 + 0.1 cos(x+y), N = 32, Δt = 1e-3, seed 11.
# The E_K mutant ⟨u·f⟩ + ½Δt⟨|f|²⟩ survives unless the support holds two modes that differ by ±(1,1). Run: julia r6-E50.jl
# Power of the WP-11 forcing-only test (ROADMAP:178) against the E_K and Omega mutants that drop ρ,
# as a function of the forcing support, with u started from rest (u stays in the span of the support).
using Random, Printf
N = 32; x = (0:N-1) .* (2π/N); X = [a for a in x, b in x]; Y = [b for a in x, b in x]
m(a) = sum(a)/length(a)
function field(rng, S)   # random vector field (solenoidal + compressive parts) on support S; returns ax, ay, curl
    ax = zeros(N,N); ay = zeros(N,N); c = zeros(N,N)
    for (p,q) in S
        A,B,C,D = randn(rng,4); s = sin.(p.*X .+ q.*Y); co = cos.(p.*X .+ q.*Y)
        ax .+= A.*s .+ B.*co; ay .+= C.*s .+ D.*co; c .+= p.*(C.*co .- D.*s) .- q.*(A.*co .- B.*s)
    end; ax, ay, c
end
band(k1,k2) = [(p,q) for p in 0:8, q in -8:8 if k1^2 <= p^2+q^2 <= k2^2 && (p > 0 || q > 0)]
ρ = 1 .+ 0.1 .* cos.(X .+ Y); Δt = 1e-3
for (name, S) in (("single mode (4,0)", [(4,0)]), ("modes (4,0),(0,4)", [(4,0),(0,4)]), ("band(4,6)", band(4,6)))
    rng = MersenneTwister(11); ux,uy,ω = field(rng,S); fx,fy,cf = field(rng,S)
    dE = (m(0.5 .* ρ .* ((ux.+Δt.*fx).^2 .+ (uy.+Δt.*fy).^2)) - m(0.5 .* ρ .* (ux.^2 .+ uy.^2)))/Δt
    dO = (m(0.5 .* (ω.+Δt.*cf).^2 ./ ρ) - m(0.5 .* ω.^2 ./ ρ))/Δt
    EK = m(0.5 .* ρ .* (ux.^2 .+ uy.^2)); Om = m(0.5 .* ω.^2 ./ ρ)
    mutE = m(ux.*fx .+ uy.*fy) + Δt/2*m(fx.^2 .+ fy.^2)
    mutO = m(ω.*cf) + Δt/2*m(cf.^2)
    textE = m(ρ.*(ux.*fx .+ uy.*fy)) + Δt/2*m(ρ.*(fx.^2 .+ fy.^2))
    @printf("%-20s Omega mutant relative to rate: %.2g\n", name, abs(mutO-dO)/abs(dO))
    @printf("%-20s text: |err|/(eps E/Δt) = %.2f   E_K mutant: |err|/(eps E/Δt) = %.3g (%.2g of rate)   Omega mutant: |err|/(eps Ω/Δt) = %.3g\n",
        name, abs(textE-dE)/(eps()*EK/Δt), abs(mutE-dE)/(eps()*EK/Δt), abs(mutE-dE)/abs(dE), abs(mutO-dO)/(eps()*Om/Δt))
end
