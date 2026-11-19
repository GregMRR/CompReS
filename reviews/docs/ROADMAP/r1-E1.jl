# r1-E1: a restart from a physical-space snapshot is not bit-for-bit (WP-12 acceptance, ADR-0009).
# The uninterrupted run keeps the spectral state û. The restarted run rebuilds it from the
# snapshot fields: û2 = truncate(rfft(irfft(û))). The script counts entries where û2 != û.
# Run: julia --project=<env with FFTW> r1-E1.jl
using FFTW, Random
N = 16; L = 2π
Random.seed!(1)
kx = [0:N÷2;]; ky = [0:N÷2-1; -N÷2:-1]
kmax = N / 3
mask = [sqrt(kx[i]^2 + ky[j]^2) <= kmax for i in eachindex(kx), j in eachindex(ky)]
f = randn(N, N)
uh = rfft(f) ./ N^2 .* mask                      # spectral state, normalisation f̂ = ⟨f e^{-ik·x}⟩
snapshot = irfft(uh .* N^2, N)                   # physical-space field written to NetCDF (exact in Float64)
uh2 = rfft(snapshot) ./ N^2 .* mask              # restart: to_spectral! + truncate!
ndiff = count(uh2 .!= uh)
println("retained modes: ", count(mask), "  modes that differ bitwise after restart: ", ndiff,
        "  max |diff| = ", maximum(abs.(uh2 .- uh)))
println(ndiff == 0 ? "bit-for-bit: yes" : "bit-for-bit: NO")
