# r2-probe-if-order: does the amended WP-7 IF order test (ROADMAP:142; y' = −2y + y², IF on)
# fail on integrating-factor stage-algebra mutants, for all four Butcher tables?
# Lawson IF-RK: Y_i = e^{c_i hL} y + h Σ_j a_ij e^{(c_i−c_j)hL} N(Y_j);
#               y⁺ = e^{hL} y + h Σ_j b_j e^{(1−c_j)hL} N(Y_j).
# Mutant M1: every exponential is e^{hL} (round 1). Mutant M2: stage factors e^{(c_i−c_j)hL} → 1.
# Run: julia r2-probe-if-order.jl   (y(0) = 1, t_end = 1, h = 0.1·2^−j)
tabs = Dict(
 "Euler"  => ([0.0;;], [1.0], [0.0]),
 "RK2"    => ([0 0; 1 0.0], [0.5, 0.5], [0, 1.0]),
 "SSPRK3" => ([0 0 0; 1 0 0; 0.25 0.25 0], [1/6, 1/6, 2/3], [0, 1, 0.5]),
 "RK4"    => ([0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0.0], [1/6, 1/3, 1/3, 1/6], [0, 0.5, 0.5, 1]))
L = -2.0; Nf(y) = y^2
exact(t) = 2 / (1 + exp(2t))                      # y' = −2y + y², y(0) = 1
function step(y, h, (A, b, c), mut)
    s = length(b); E(x) = mut == 1 ? exp(h * L) : exp(x * h * L)
    Y = zeros(s); F = zeros(s)
    for i in 1:s
        Y[i] = E(c[i]) * y + h * sum((mut == 2 ? 1.0 : E(c[i] - c[j])) * A[i, j] * F[j] for j in 1:i-1; init = 0.0)
        F[i] = Nf(Y[i])
    end
    E(1.0) * y + h * sum(b[j] * E(1 - c[j]) * F[j] for j in 1:s)
end
order(tab, mut) = begin
    errs = map(0:3) do j
        h = 0.1 / 2^j; y = 1.0
        for _ in 1:round(Int, 1 / h); y = step(y, h, tab, mut); end
        abs(y - exact(1.0))
    end
    round.(log2.(errs[1:end-1] ./ errs[2:end]), digits = 2)
end
for name in ("Euler", "RK2", "SSPRK3", "RK4")
    println(rpad(name, 7), " correct ", order(tabs[name], 0), "  M1 ", order(tabs[name], 1), "  M2 ", order(tabs[name], 2))
end
