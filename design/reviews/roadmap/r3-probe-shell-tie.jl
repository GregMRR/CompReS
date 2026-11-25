# r3-probe-shell-tie: ROADMAP:151 states that a mode exactly on K_n = (n+½)Δk "belongs to shell n+1"
# (GLOSSARY:62, spec:27: n − ½ ≤ |q|/Δk < n + ½). ROADMAP:103 (WP-1) builds `shell_index` with no
# evaluation rule, while ROADMAP:102 requires integer arithmetic for the truncation tie.
# Rectangle family: Lx = p·c, Ly = q·c, Nx = p·m, Ny = q·m (equal spacing), Δk = 2π/max(Lx, Ly).
# Compare, over retained modes, the exact shell (integer arithmetic) with two float rules:
#   F1: floor(|k|/Δk + 1/2)  (half-up, the glossary rule in floats)
#   F2: round(Int, |k|/Δk)   (Julia default RoundNearest, ties to even)
# Run: julia r3-probe-shell-tie.jl
function exact_shell(i, j, Nx, Ny)            # |k|/Δk = sqrt((i·Nm/Nx)² + (j·Nm/Ny)²), Nm = max(Nx, Ny)
    Nm = max(Nx, Ny); A = i^2 * Nm^2 * Ny^2 + j^2 * Nm^2 * Nx^2; D = Nx^2 * Ny^2   # (|k|/Δk)² = A/D
    n = 0; while 4A >= (2n + 1)^2 * D; n += 1; end; n            # smallest n with |k|/Δk < n + ½
end
ex = Dict{String,Int}(); tot = 0; ties = 0; bad1 = 0; bad2 = 0; first = ""
for c in (2π, 1.0, 0.3, 1.7), (p, q) in ((3, 2), (4, 3), (5, 3), (2, 1), (5, 4)), m in (6, 12, 24)
    Lx, Ly, Nx, Ny = p * c, q * c, p * m, q * m
    Δk = 2π / max(Lx, Ly)
    for i in 0:Nx÷2, j in -Ny÷2:Ny÷2-1
        9 * (i^2 * Ny^2 + j^2 * Nx^2) < Nx^2 * Ny^2 || continue
        global tot += 1
        k = sqrt((i * 2π / Lx)^2 + (j * 2π / Ly)^2); r = k / Δk
        e = exact_shell(i, j, Nx, Ny)
        Nm = max(Nx, Ny); A = i^2 * Nm^2 * Ny^2 + j^2 * Nm^2 * Nx^2; D = Nx^2 * Ny^2
        istie = 4A == (2e - 1)^2 * D; global ties += istie
        f1 = floor(Int, r + 0.5); f2 = round(Int, r)
        if f1 != e; global bad1 += 1; first == "" && (global first = "F1 c=$c Lx=$(p)c Ly=$(q)c Nx=$Nx Ny=$Ny mode ($i,$j): |k|/Δk=$r exact shell $e, float $f1"); end
        f2 != e && (global bad2 += 1)
    end
end
println("retained modes: $tot, exact ties on a shell edge: $ties")
println("F1 (float half-up) disagrees with the exact rule: $bad1;  F2 (round, ties to even): $bad2")
println(first)
println("Lx = 3π, Ly = 2π, mode (2, 1): |k|/Δk = ", sqrt((2 * 2π / (3π))^2 + (2π / (2π))^2) / (2π / (3π)),
        "; round -> ", round(Int, sqrt((2 * 2π / (3π))^2 + 1.0) / (2π / (3π))), "; exact shell 3")
