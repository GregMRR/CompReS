# r1-E6: the disc |k| <= k_max = N/3 (ADR-0003, WP-1) is not alias-free when 3 divides N.
# 1D cut along x is enough: the triad (K,0)+(K,0) -> (2K,0) folds onto (-K,0), which is retained.
# Run: julia --project=<env with FFTW> r1-E6.jl
using FFTW
for N in (12, 16, 384)
    K = N / 3
    x = 2π .* (0:N-1) ./ N
    Kint = floor(Int, K)                 # largest retained integer wavenumber on the axis
    u = cos.(Kint .* x)                  # retained mode
    p = rfft(u .* u) ./ N                # pseudospectral product
    # exact product cos^2 = 1/2 + cos(2Kint x)/2 has nothing at Kint; any value there is aliasing
    println("N = $N  k_max = $(round(K; digits=3))  retained axis mode $Kint  |aliased coeff at k=$Kint| = ",
            abs(p[Kint+1]))
end

# Incompressible model: the aliased term of the self-triad is along k, so the projection removes it.
# Energy and enstrophy fluxes at k_max for a random solenoidal truncated field, N = 12 (3 divides N).
using Random
function inc_fluxes(N; seed=3)
    Random.seed!(seed)
    kx = [0:N÷2;]; ky = [0:N÷2-1; -N÷2:-1]
    KX = [kx[i] for i in eachindex(kx), j in eachindex(ky)]; KY = [ky[j] for i in eachindex(kx), j in eachindex(ky)]
    K2 = KX .^ 2 .+ KY .^ 2; mask = sqrt.(K2) .<= N / 3
    ψ = rfft(randn(N, N)) .* mask ./ N^2                    # streamfunction -> solenoidal u
    ux = im .* KY .* ψ; uy = -im .* KX .* ψ; ω = im .* (KX .* uy .- KY .* ux)
    P(a) = irfft(a .* N^2, N)
    F(a) = rfft(a) ./ N^2
    Nx = F(-P(ω) .* P(uy)); Ny = F(P(ω) .* P(ux))            # ω ẑ × u
    K2s = copy(K2); K2s[1, 1] = 1
    kd = (KX .* Nx .+ KY .* Ny) ./ K2s
    Nx = (Nx .- KX .* kd) .* mask; Ny = (Ny .- KY .* kd) .* mask # projection + truncation
    w = fill(2.0, size(mask)); w[1, :] .= 1; w[end, :] .= 1     # half-plane weights
    Pi_E = sum(w .* real.(conj.(ux) .* Nx .+ conj.(uy) .* Ny))
    Nω = im .* (KX .* Ny .- KY .* Nx)
    Pi_Z = sum(w .* real.(conj.(ω) .* Nω))
    Pi_E, Pi_Z
end
println("incompressible N = 12: Pi_E(k_max), Pi_Z(k_max) = ", inc_fluxes(12))
