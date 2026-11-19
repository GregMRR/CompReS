# r3-E35: WP-14 (ROADMAP:190) requires, for the cellular flow at Ma = 0.1 with Isothermal and
# Polytropic(γ = 1.4), max|∂_tθ̂(0)| / max|k_ik_j (u_iu_j)^| < 1e-12 at N = 64.
# This probe follows the amended ADR-0013 procedure exactly: solve ∇²h = −∂i∂j(uiuj) spectrally,
# root-find h₀ so that ⟨ρ(h + h₀)⟩ = 1, ρ = density_from_enthalpy(h + h₀), truncate ρ̂ to the disc.
# Then the discrete ∂_tθ̂ = ik·trunc[−(u·∇u)^ − ik ĥ(ρ)], with h(ρ) evaluated on the grid from ρ̂.
# Ma = u_rms/c_s(1) (spec), L = 2π, c_s(1) chosen so that u_rms = 1 (scale-free otherwise).
# Run: julia --project=<env with FFTW> r3-E35.jl
using FFTW, Random, Printf
struct Iso; cs2::Float64; end
struct Poly; K::Float64; γ::Float64; end
hfun(e::Iso, ρ) = e.cs2 * log(ρ);              ρfun(e::Iso, h) = exp(h / e.cs2)
hfun(e::Poly, ρ) = e.K * e.γ / (e.γ - 1) * ρ^(e.γ - 1)
ρfun(e::Poly, h) = ((e.γ - 1) * h / (e.K * e.γ))^(1 / (e.γ - 1))
function probe(N, eos; Ma = 0.1)
    x = (0:N-1) .* (2π/N); X = [a for a in x, b in x]; Y = [b for a in x, b in x]
    kx = [0:N÷2;]; ky = [0:N÷2-1; -N÷2:-1]
    KX = [a for a in kx, b in ky]; KY = [b for a in kx, b in ky]; K2 = KX.^2 .+ KY.^2
    tr = 9 .* K2 .< N^2
    F(f) = rfft(f) ./ N^2; B(f̂) = irfft(f̂ .* N^2, N)
    ux = sin.(X) .* cos.(Y); uy = -cos.(X) .* sin.(Y)
    s = 1 / sqrt(sum(ux.^2 + uy.^2) / N^2); ux .*= s; uy .*= s          # u_rms = 1
    ûx = F(ux) .* tr; ûy = F(uy) .* tr
    Q = (KX.^2 .* F(ux.*ux) .+ 2 .* KX .* KY .* F(ux.*uy) .+ KY.^2 .* F(uy.*uy)) .* tr
    ĥ = zeros(ComplexF64, size(Q)); m = K2 .> 0; ĥ[m] = -Q[m] ./ K2[m]
    h = B(ĥ)
    # h₀ by bisection on ⟨ρ(h + h₀)⟩ = 1, around h(1)
    g(h0) = sum(ρfun.(Ref(eos), h .+ h0)) / N^2 - 1
    lo, hi = hfun(eos, 0.5), hfun(eos, 2.0)
    for _ in 1:200; mid = (lo + hi) / 2; (g(mid) > 0 ? (hi = mid) : (lo = mid)); end
    ρ = ρfun.(Ref(eos), h .+ (lo + hi) / 2)
    ρ̂ = F(ρ) .* tr
    hh = F(hfun.(Ref(eos), B(ρ̂)))
    dux = B(im .* KX .* ûx); duy = B(im .* KY .* ûx); dvx = B(im .* KX .* ûy); dvy = B(im .* KY .* ûy)
    Nx̂ = F(ux .* dux .+ uy .* duy); Nŷ = F(ux .* dvx .+ uy .* dvy)
    tx = (-Nx̂ .- im .* KX .* hh) .* tr; ty = (-Nŷ .- im .* KY .* hh) .* tr
    maximum(abs, im .* KX .* tx .+ im .* KY .* ty) / maximum(abs, Q)
end
cs = 1 / 0.1
for (name, eos) in (("Isothermal", Iso(cs^2)), ("Polytropic(γ=1.4)", Poly(cs^2 / 1.4, 1.4)))
    print(rpad(name, 18))
    for N in (16, 32, 64, 128, 256)
        @printf("  N=%-3d %.1e", N, probe(N, eos))
    end
    println()
end
