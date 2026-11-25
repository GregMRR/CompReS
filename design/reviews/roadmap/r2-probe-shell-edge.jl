# r2-probe-shell-edge: WP-8 (ROADMAP:148) reports Pi(K_n) = sum over shells 0…n at K_n = (n+½)Δk.
# The shells (spec:27) are n−½ ≤ |q|/Δk < n+½, so the shell sum is Σ_{|q| < K_n}. The spec flux
# (spec:30) is Σ_{|q| ≤ K}. They differ when a mode lies on an edge, which happens on rectangles.
# Run: julia r2-probe-shell-edge.jl
let Lx = 3π, Ly = 2π
    Δk = 2π / max(Lx, Ly); q = 2π / Ly                     # mode (0, ±1)
    n = floor(Int, q / Δk + 0.5); K1 = 1.5 * Δk
    println("Lx = 3π, Ly = 2π: Δk = ", Δk, ", mode (0,±1): |q|/Δk = ", q / Δk, ", shell n = ", n)
    println("  K_1 = 1.5Δk = ", K1, "; |q| ≤ K_1 (spec:30 includes it) = ", q <= K1,
            "; in shells 0…1 (WP-8 includes it) = ", n <= 1)
end
