# r2-probe-tie: the strict truncation |k| < k_c, k_c = (Nx/3)(2π/Lx) (ROADMAP:101) relies on the
# axis mode (0, Ny/3) giving |k| = k_c exactly. On a rectangle, ky = (Ny/3)(2π/Ly) is a different
# floating-point expression, so the tie can round either way. The grid also accepts spacings that
# differ by up to 1e-12 relative (ROADMAP:100). Search natural inputs Ly = Lx·Ny/Nx and
# Lx = 2π·a for grids the plan accepts where (0, Ny/3) is retained (aliased self-triad back).
# Also: the accept example of WP-1 in Float32.
# Run: julia r2-probe-tie.jl
accepted(Lx, Nx, Ly, Ny) = abs(Lx/Nx - Ly/Ny) <= 1e-12 * Lx/Nx
hits = Tuple[]
for a in (1.0, 2.0, 0.5, 0.3, 0.7, 1.1, 1.3, 2.5, 3.0), Nx in 6:6:600, Ny in 6:6:600
    Nx == Ny && continue
    Lx = 2π * a; Ly = Lx * Ny / Nx
    accepted(Lx, Nx, Ly, Ny) || continue
    kc = (Nx / 3) * (2π / Lx); ky = (Ny ÷ 3) * (2π / Ly)
    if ky^2 < kc^2 || ky < kc
        push!(hits, (a, Nx, Ny, ky - kc, ky^2 - kc^2))
    end
end
println("accepted rectangles with (0, Ny/3) retained: ", length(hits))
for h in hits[1:min(5, end)]
    println("  Lx = 2π·", h[1], ", Nx = ", h[2], ", Ly = Lx·Ny/Nx, Ny = ", h[3],
            ":  ky − k_c = ", h[4], ",  ky² − k_c² = ", h[5])
end
# A decimal-literal rectangle next to the WP-1 accept example:
let Lx = 0.3, Nx = 30, Ly = 0.9, Ny = 90
    kc = (Nx / 3) * (2π / Lx); ky = (Ny ÷ 3) * (2π / Ly)
    println("Lx = 0.3, Nx = 30, Ly = 0.9, Ny = 90: accepted = ", accepted(Lx, Nx, Ly, Ny),
            ", mode (0, 30) retained = ", ky < kc, " (ky − k_c = ", ky - kc, ")")
    println("  self-triad (0,30)+(0,30) = (0,60) folds to (0,", 60 - Ny, "), retained = ", ky < kc)
end
# Self-triad of q = (0, Ny/3): 2q = (0, 2Ny/3) folds to (0, −Ny/3), which is then retained.
# WP-1 accept example in Float32 (the plan does not state the precision of the spacing check):
for T in (Float64, Float32)
    Lx, Ly = T(0.3), T(0.1)
    println(T, ": Lx = 0.3, Nx = 30, Ly = 0.1, Ny = 10 accepted = ", accepted(Lx, 30, Ly, 10),
            "  (relative spacing difference ", abs(Lx/30 - Ly/10) / (Lx/30), ")")
end
