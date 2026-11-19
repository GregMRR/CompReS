# r4-probe-shell: the amended WP-1 shell rule (ROADMAP:105):
#   mode (i, j) is in the smallest n with 4·Nm²(i²Ny² + j²Nx²) < (2n+1)²Nx²Ny², Nm = max(Nx, Ny).
# (1) the acceptance example Lx = 5·2π, Ly = 4·2π, Nx = 60, Ny = 48, mode (6, −14): shell and float value;
# (2) the size of the largest product, over retained modes and over all modes of the real-FFT layout,
#     against typemax(Int64), to locate the N at which Int64 overflows.
# Run: julia r4-probe-shell.jl
function shell(i, j, Nx, Ny)
    Nm = big(max(Nx, Ny)); A = 4Nm^2 * (big(i)^2 * Ny^2 + big(j)^2 * Nx^2); D = big(Nx)^2 * Ny^2
    n = 0; while A >= (2n + 1)^2 * D; n += 1; end; (n, max(A, (2n + 1)^2 * D))
end
Lx, Ly, Nx, Ny = 5 * 2π, 4 * 2π, 60, 48
println("(6,−14): exact shell ", shell(6, -14, Nx, Ny)[1], "; float |k|/Δk = ",
        sqrt((6 * 2π / Lx)^2 + (14 * 2π / Ly)^2) / (2π / max(Lx, Ly)),
        "; retained: ", 9 * (36 * Ny^2 + 196 * Nx^2) < Nx^2 * Ny^2)
# Largest product evaluated by the 'smallest n' loop (it tests (2n+1)²D up to the final n).
for N in (1024, 1100, 1200, 1600, 1650, 1700)
    mr = big(0); ma = big(0)
    for i in (N ÷ 2 - 2):(N ÷ 2), j in (N ÷ 2 - 2):(N ÷ 2)      # corners: largest |k|
        p = shell(i, j, N, N)[2]; ma = max(ma, p)
    end
    # largest retained |k|: scan the circle edge
    for i in 0:N÷2
        j = isqrt((N^2 - 1) ÷ 9 - i^2 < 0 ? 0 : (N^2 - 1) ÷ 9 - i^2)
        9 * (i^2 + j^2) < N^2 || continue
        mr = max(mr, shell(i, j, N, N)[2])
    end
    println("N = $N: retained max product / typemax(Int64) = ", Float64(mr / typemax(Int64)),
            ";  all modes: ", Float64(ma / typemax(Int64)))
end
