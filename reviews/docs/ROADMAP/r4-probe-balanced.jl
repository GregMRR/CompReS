# r4-probe-balanced: WP-14 (ROADMAP:192) bound 1e-9 at N = 64, and its claim that the wrong procedure
# (shift ρ after the inversion) gives about 1e-5. Cellular flow, Ma = 0.1, u_rms = 1, L = 2π.
# Correct: ρ = ρ(h + h₀), h₀ by bisection with ⟨ρ⟩ = 1. Wrong: ρ = ρ(h + c) − ⟨ρ(h + c)⟩ + 1 for a
# constant c: c = h(1) (h zero-mean), and c = h(1) + δ for δ chosen so the wrong mean is near 1.
# Run: julia --project=<env with FFTW> r4-probe-balanced.jl
using FFTW, Printf
struct Iso; cs2::Float64; end
struct Poly; K::Float64; γ::Float64; end
hfun(e::Iso, ρ) = e.cs2 * log(ρ);  ρfun(e::Iso, h) = exp(h / e.cs2)
hfun(e::Poly, ρ) = e.K * e.γ / (e.γ - 1) * ρ^(e.γ - 1)
ρfun(e::Poly, h) = ((e.γ - 1) * h / (e.K * e.γ))^(1 / (e.γ - 1))
function setup(N)
    x = (0:N-1) .* (2π/N); X = [a for a in x, b in x]; Y = [b for a in x, b in x]
    kx = [0:N÷2;]; ky = [0:N÷2-1; -N÷2:-1]
    KX = [a for a in kx, b in ky]; KY = [b for a in kx, b in ky]; K2 = KX.^2 .+ KY.^2
    tr = 9 .* K2 .< N^2
    F(f) = rfft(f) ./ N^2; B(f̂) = irfft(f̂ .* N^2, N)
    ux = sin.(X) .* cos.(Y); uy = -cos.(X) .* sin.(Y)
    s = 1 / sqrt(sum(ux.^2 + uy.^2) / N^2); ux .*= s; uy .*= s
    ûx = F(ux) .* tr; ûy = F(uy) .* tr
    Q = (KX.^2 .* F(ux.*ux) .+ 2 .* KX .* KY .* F(ux.*uy) .+ KY.^2 .* F(uy.*uy)) .* tr
    ĥ = zeros(ComplexF64, size(Q)); m = K2 .> 0; ĥ[m] = -Q[m] ./ K2[m]
    (; N, KX, KY, tr, F, B, ux, uy, ûx, ûy, Q, h = B(ĥ))
end
function resid(S, eos, ρ)
    (; KX, KY, tr, F, B, ux, uy, ûx, ûy, Q) = S
    ρ̂ = F(ρ) .* tr; hh = F(hfun.(Ref(eos), B(ρ̂)))
    dux = B(im .* KX .* ûx); duy = B(im .* KY .* ûx); dvx = B(im .* KX .* ûy); dvy = B(im .* KY .* ûy)
    Nx̂ = F(ux .* dux .+ uy .* duy); Nŷ = F(ux .* dvx .+ uy .* dvy)
    tx = (-Nx̂ .- im .* KX .* hh) .* tr; ty = (-Nŷ .- im .* KY .* hh) .* tr
    maximum(abs, im .* KX .* tx .+ im .* KY .* ty) / maximum(abs, Q)
end
cs = 10.0
for (name, eos) in (("Isothermal", Iso(cs^2)), ("Polytropic(1.4)", Poly(cs^2 / 1.4, 1.4)))
    S = setup(64); h = S.h; N = S.N
    g(h0) = sum(ρfun.(Ref(eos), h .+ h0)) / N^2 - 1
    lo, hi = hfun(eos, 0.5), hfun(eos, 2.0)
    for _ in 1:200; mid = (lo + hi) / 2; (g(mid) > 0 ? (hi = mid) : (lo = mid)); end
    h0 = (lo + hi) / 2
    ρc = ρfun.(Ref(eos), h .+ h0)
    @printf("%-16s correct: %.1e | wrong, c = h(1): ", name, resid(S, eos, ρc))
    for c in (hfun(eos, 1.0), h0 + 1e-3 * (h0 - hfun(eos, 1.0)), h0 + 1e-6 * abs(h0 - hfun(eos, 1.0)))
        ρw = ρfun.(Ref(eos), h .+ c); ρw = ρw .- sum(ρw) / N^2 .+ 1
        @printf("%.1e  ", resid(S, eos, ρw))
    end
    @printf(" (h0 − h(1) = %.2e)\n", h0 - hfun(eos, 1.0))
end
