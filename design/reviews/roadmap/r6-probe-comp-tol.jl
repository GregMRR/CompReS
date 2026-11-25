# r6-probe-comp-tol: margin of the WP-11 tolerance 10·eps(T)·E_kin/Δt (ROADMAP:178) for a correct implementation:
# spectral state, physical fields by irfft, grid means, band(4,6) from rest, ρ = 1 + 0.1 cos(x+y), 1000 steps.
# Needs FFTW: julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("FFTW"); include("r6-probe-comp-tol.jl")'
# As comp32.jl, but the state is spectral (û += Δt f̂ in Complex{T}) and physical fields come from irfft,
# as a WP-11 implementation computes them. E_K and inj are grid means of the physical fields.
using FFTW, Random, Printf
function run(T, N, nsteps; seed = 3)
    x = T.((0:N-1) .* (2π/N)); X = [a for a in x, b in x]; Y = [b for a in x, b in x]
    ρ = 1 .+ T(0.1) .* cos.(X .+ Y); m(a) = sum(a)/length(a); Δt = T(0.01)
    S = [(p, q) for p in 0:8, q in -8:8 if 16 <= p^2+q^2 <= 36 && (p > 0 || q > 0)]
    idx(p, q) = (p + 1, mod(q, N) + 1)
    rng = MersenneTwister(seed); ûx = zeros(Complex{T}, N÷2+1, N); ûy = copy(ûx); f̂x = copy(ûx); f̂y = copy(ûx)
    c = randn(rng, Complex{T}, 2, length(S)); a = exp(-Δt); bb = sqrt(1 - a^2); worst = zero(T)
    phys(û) = irfft(û, N) .* T(N^2)  # f̂ = ⟨f e^{-ikx}⟩ normalisation
    for n in 1:nsteps
        c .= a .* c .+ bb .* randn(rng, Complex{T}, size(c)); f̂x .= 0; f̂y .= 0
        for (j, (p, q)) in enumerate(S)
            i, k = idx(p, q); f̂x[i, k] = c[1, j]; f̂y[i, k] = c[2, j]
            if p == 0; i2, k2 = idx(0, -q); f̂x[i2, k2] = conj(c[1, j]); f̂y[i2, k2] = conj(c[2, j]); end
        end
        ux, uy, fx, fy = phys(ûx), phys(ûy), phys(f̂x), phys(f̂y)
        inj = m(ρ .* (ux.*fx .+ uy.*fy)) + Δt/2 * m(ρ .* (fx.^2 .+ fy.^2))
        E0 = m(T(0.5) .* ρ .* (ux.^2 .+ uy.^2))
        ûx .+= Δt .* f̂x; ûy .+= Δt .* f̂y
        E1 = m(T(0.5) .* ρ .* (phys(ûx).^2 .+ phys(ûy).^2))
        worst = max(worst, abs(inj - (E1 - E0)/Δt) / (eps(T)*E1/Δt))
    end
    worst
end
for T in (Float64, Float32), N in (16, 32, 64), s in (3, 4)
    @printf("%-8s N = %3d seed %d  max over 1000 steps |inj - ΔE_K/Δt| / (eps(T) E_K(t_{n+1})/Δt) = %.2f\n", T, N, s, run(T, N, 1000; seed = s))
end
