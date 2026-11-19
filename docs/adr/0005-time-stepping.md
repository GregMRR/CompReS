# Time stepping

Date: 2026-11-19. Status: accepted.

One generic explicit Runge–Kutta driver reads a Butcher table. The tables are Euler, RK2 (Heun), SSP-RK3 and RK4. Each scheme runs with or without an integrating factor (IF) on the linear terms. The IF is a per-mode linear propagator `exp(LΔt)`; it is a scalar for diagonal operators and a small matrix per mode for the compressible Helmholtz split (ADR-0004).

The time step `Δt` is fixed. The driver computes the CFL number and stops the run when it exceeds 1 (ADR-0010). A fixed `Δt` keeps stochastic forcing and exact restarts simple.

Stochastic forcing updates once per step, before the RK stages, and stays constant inside the stages.

## Open question

At low Mach the acoustic CFL limit `Δt ≲ Δx/(c_s + u)` is about `1/Ma` times smaller than the advective limit. The linear acoustic operator (with `ρ₀ = 1`, `c_s(ρ₀)`) is a per-mode 2×2 map on `(ρ̂, û^c)`, so the IF propagator interface can absorb it later. Decide this after the compressible solver works and the benchmarks (WP-16) measure the cost.
