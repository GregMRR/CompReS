# r1-E2: the planned blocking tests do not exercise the integrating-factor (IF) stage algebra.
# Lawson IF-RK4 for y' = L y + N(y). Mutation M: every stage uses exp(LΔt) instead of exp(c_i LΔt).
# Test A (Kolmogorov decay, N ≡ 0): correct and mutant both exact.
# Test B (conservation, ν = 0 so L = 0): correct and mutant identical.
# Test C (scalar ODE with L ≠ 0 and N ≠ 0): the mutant loses its order.
# Run: julia r1-E2.jl
function ifrk4(y, L, N, dt, nsteps; mutant=false)
    E(c) = mutant ? exp(L * dt) : exp(c * L * dt)
    for _ in 1:nsteps
        k1 = N(y)
        y2 = E(0.5) * (y + dt / 2 * k1); k2 = N(y2)
        y3 = E(0.5) * y + dt / 2 * k2;   k3 = N(y3)
        y4 = E(1.0) * y + dt * E(0.5) * k3; k4 = N(y4)
        y = E(1.0) * y + dt / 6 * (E(1.0) * k1 + 2E(0.5) * (k2 + k3) + k4)
    end
    y
end
T = 1.0
# A: Kolmogorov decay, nonlinear term zero, decay rate ν k^2 = 0.3
for m in (false, true)
    y = ifrk4(1.0, -0.3, y -> 0.0, 0.1, 10; mutant=m)
    println("A mutant=$m  error = ", abs(y - exp(-0.3T)))
end
# B: conservation test, L = 0, skew-symmetric rotation as a stand-in for the nonlinear term
for m in (false, true)
    y = ifrk4(1.0 + 0im, 0.0, y -> im * y, 0.1, 10; mutant=m)
    println("B mutant=$m  y = ", y)
end
# C: y' = -2y + y^2, y(0) = 1, exact y = 2/(1 + e^{2t})
yex = 2 / (1 + exp(2T))
for m in (false, true)
    errs = [abs(ifrk4(1.0, -2.0, y -> y^2, T / n, n; mutant=m) - yex) for n in (20, 40, 80, 160)]
    slopes = [log2(errs[i] / errs[i+1]) for i in 1:3]
    println("C mutant=$m  errors = ", errs, "  observed orders = ", round.(slopes; digits=2))
end
