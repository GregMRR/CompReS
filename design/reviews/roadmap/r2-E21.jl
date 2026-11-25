# r2-E21: WP-14 (ROADMAP:183) requires ∂_tθ(0) = 0 "to round-off" for μ = ζ = κ = 0.
# ADR-0013:5: solve ∇²h = −∂_i∂_j(u_iu_j) spectrally, apply density_from_enthalpy, "shift the
# constant so that ⟨ρ⟩ = 1". Isothermal EOS, c_s = 10, L = 2π, flow scaled to Ma = 0.1.
# Variants: (A) additive shift ρ += 1 − ⟨ρ⟩ (literal reading); (B) shift of h, i.e. ρ ./= ⟨ρ⟩
# (exact for isothermal). Then the state ρ̂ is (t) truncated or (u) not truncated.
# Discrete ∂_tθ̂ = ik·trunc[−(u·∇u)^ − ik ĥ(ρ)], with h(ρ) evaluated on the grid from the state ρ̂.
# Reported: max|∂_tθ̂| / max|k_ik_j (u_iu_j)^|.
# Run: julia --project=<env with FFTW> r2-probe-balance.jl
using FFTW, Random
function probe(N, uxf, uyf; Ma = 0.1, cs = 10.0)
    x = (0:N-1) .* (2π/N); X = [a for a in x, b in x]; Y = [b for a in x, b in x]
    kx = [0:N÷2;]; ky = [0:N÷2-1; -N÷2:-1]
    KX = [a for a in kx, b in ky]; KY = [b for a in kx, b in ky]; K2 = KX.^2 .+ KY.^2
    tr = K2 .< (N/3)^2
    F(f) = rfft(f) ./ N^2; B(f̂) = irfft(f̂ .* N^2, N)
    ûx = F(uxf.(X, Y)) .* tr; ûy = F(uyf.(X, Y)) .* tr
    ux = B(ûx); uy = B(ûy); s = Ma * cs / sqrt(sum(ux.^2 + uy.^2) / N^2)
    ûx .*= s; ûy .*= s; ux .*= s; uy .*= s
    Q = KX.^2 .* F(ux.*ux) .+ 2 .* KX .* KY .* F(ux.*uy) .+ KY.^2 .* F(uy.*uy)   # k_ik_j (u_iu_j)^
    Q .*= tr
    ĥ = zeros(ComplexF64, size(Q)); ĥ[K2 .> 0] = -Q[K2 .> 0] ./ K2[K2 .> 0]    # ∇²h = −∂i∂j(uu)
    h = B(ĥ); ρ0 = exp.(h ./ cs^2)
    dux = B(im .* KX .* ûx); duy = B(im .* KY .* ûx); dvx = B(im .* KX .* ûy); dvy = B(im .* KY .* ûy)
    Nx̂ = F(ux .* dux .+ uy .* duy); Nŷ = F(ux .* dvx .+ uy .* dvy)
    out = Dict{String,Float64}()
    for (vn, ρ) in (("A additive", ρ0 .+ (1 - sum(ρ0) / N^2)), ("B scaled  ", ρ0 ./ (sum(ρ0) / N^2)))
        for (tn, trunc) in (("trunc", true), ("untrunc", false))
            ρ̂ = F(ρ); trunc && (ρ̂ .*= tr)
            hh = F(cs^2 .* log.(B(ρ̂)))
            tx = (-Nx̂ .- im .* KX .* hh) .* tr; ty = (-Nŷ .- im .* KY .* hh) .* tr
            out[vn * " " * tn] = maximum(abs, im .* KX .* tx .+ im .* KY .* ty) / maximum(abs, Q)
        end
    end
    out
end
cell(x, y) = sin(x) * cos(y); cellv(x, y) = -cos(x) * sin(y)
two(x, y) = sin(x) * cos(y) + 0.6 * cos(2y + 1.0); twov(x, y) = -cos(x) * sin(y) + 0.4 * sin(3x)
for (name, f, g) in (("cellular flow", cell, cellv), ("three-mode flow", two, twov)), N in (16, 32, 64)
    r = probe(N, f, g)
    println(rpad(name, 16), " N = ", lpad(N, 2), ":  ", join(["$k: $(round(v, sigdigits = 2))" for (k, v) in sort(collect(r))], "  "))
end
