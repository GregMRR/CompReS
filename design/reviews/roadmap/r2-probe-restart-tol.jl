# r2-probe-restart-tol: WP-12 (ROADMAP:172) and ADR-0009:15 require the state rebuilt from a
# snapshot to equal the saved state "to relative 1e-14", with no norm and no float type stated.
# Rebuild: û2 = truncate(rfft(irfft(û))). Field: random phases, |û|² ∝ k^4 exp(−k²/8)/k (target API).
# Run: julia --project=<env with FFTW> r2-probe-restart-tol.jl   (seed 4)
using FFTW, Random
function rt(T, N)
    rng = MersenneTwister(4)
    kx = [0:N÷2;]; ky = [0:N÷2-1; -N÷2:-1]
    K = [hypot(a, b) for a in kx, b in ky]; tr = K .< N/3
    amp = [k == 0 ? 0.0 : sqrt(k^3 * exp(-k^2 / 8)) for k in K]
    û = Complex{T}.(amp .* cis.(2π .* rand(rng, size(K)...)) .* tr)
    û = rfft(irfft(û, N)) .* tr          # make it an exact Hermitian, truncated state first
    û2 = rfft(irfft(û, N)) .* tr
    nrm = sqrt(sum(abs2, û2 .- û) / sum(abs2, û))
    nz = abs.(û) .> 0
    permode = maximum(abs.(û2 .- û)[nz] ./ abs.(û)[nz])
    T, N, nrm, permode
end
for (T, N) in ((Float64, 64), (Float64, 512), (Float32, 64), (Float32, 512))
    t, n, a, b = rt(T, N)
    println(rpad(string(t), 8), " N = ", lpad(n, 3), ":  norm-relative ", round(a, sigdigits = 2),
            "   max per-mode relative ", round(b, sigdigits = 2))
end
