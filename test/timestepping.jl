# Tests of the time stepping (WP-7): convergence orders on a one-mode problem with and without the integrating
# factor, conservation and Kolmogorov decay (ADR-0011), linear decay of the compressible fields and the
# injection rates of a forcing-only run.
module TimesteppingTests

using Test
using KernelAbstractions
using CompReS
using CompReS: IncompressibleRates, TimeStepper, kinetic_energy, enstrophy

const SCHEMES = ((Euler(), 1), (RK2(), 2), (SSPRK3(), 3), (RK4(), 4))

# One-mode problem y' = −2y + y²: the linear rate 2 goes through the production IncompressibleRates, and rhs!
# returns the nonlinear part y² for both fields.
struct OneModeProblem{R}
    forcing::NoForcing
    rates::R
end

function CompReS.rhs!(dû, state, ::OneModeProblem, ws)
    dû.ux .= state.fields.ux .^ 2
    dû.uy .= state.fields.uy .^ 2
    return dû
end

# Exact solution of y' = −2y + y²: 1/y follows z' = 2z − 1.
exact_one_mode(y0, t) = 1 / (1 / 2 + (1 / y0 - 1 / 2) * exp(2t))

# Maximum error at t = 1 of the one-mode problem with the step Δt.
function one_mode_error(backend, scheme, integrating_factor, Δt)
    y0 = (1.0, 0.5)
    field(y) = copyto!(KernelAbstractions.allocate(backend, ComplexF64, 1, 1), fill(ComplexF64(y), 1, 1))
    state = State((ux = field(y0[1]), uy = field(y0[2])))
    problem = OneModeProblem(NoForcing(), IncompressibleRates(copyto!(KernelAbstractions.allocate(backend, Float64, 1, 1),
                                                                     fill(2.0, 1, 1))))
    ws = (stepper = TimeStepper(scheme, integrating_factor, state.fields),)
    for _ in 1:round(Int, 1 / Δt)
        step!(state, problem, ws, Δt)
    end
    @test state.t ≈ 1
    return max(abs(only(Array(state.fields.ux)) - exact_one_mode(y0[1], 1.0)),
               abs(only(Array(state.fields.uy)) - exact_one_mode(y0[2], 1.0)))
end

# Least-squares slope of log(err) against log(Δt).
function loglog_slope(Δts, errs)
    x, y = log.(Δts), log.(errs)
    return sum((x .- sum(x) / length(x)) .* y) / sum((x .- sum(x) / length(x)) .^ 2)
end

# Stability polynomial R(z) = Σ_{j ≤ r} z^j/j! of an s-stage scheme of order r = s (all four schemes).
stability_polynomial(r, z) = sum(z^j / factorial(j) for j in 0:r)

relative_error(a, b) = sqrt(sum(abs2, Array(a) .- Array(b)) / sum(abs2, Array(b)))

function run_tests(backend)
    @testset "time stepping" begin
        @testset "one-mode convergence order: $(nameof(typeof(scheme))), IF $integrating_factor" for (scheme, r) in SCHEMES,
                                                                                                    integrating_factor in (false, true)
            # The order of the two finest steps; at coarser steps RK4 with the integrating factor converges
            # faster than Δt^4 (pre-asymptotic range).
            errs = [one_mode_error(backend, scheme, integrating_factor, Δt) for Δt in (1 / 20, 1 / 40, 1 / 80, 1 / 160)]
            @test issorted(errs; rev = true)
            @test abs(log2(errs[end - 1] / errs[end]) - r) < 0.1
        end

        @testset "Kolmogorov decay: $(nameof(typeof(scheme))), IF $integrating_factor, $T" for (scheme, r) in SCHEMES,
                                                                                            integrating_factor in (false, true),
                                                                                            T in (Float64, Float32)
            # u = A sin(3y) x̂ is a steady solution of the Euler equations, so each step multiplies it by the
            # exact factor e^{−νk²Δt} with the integrating factor and by R(−νk²Δt) without it.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T, backend)
            ν, n, Δt, nsteps = 0.1, 3, 0.02, 50
            state = State(kolmogorov(g; A = 1.0, n))
            initial = map(copy, state.fields)
            ws = Workspace(g, state; scheme, integrating_factor)
            problem = Problem(g, Incompressible(; ν), NoForcing())
            E0 = kinetic_energy(state, problem, ws)
            for _ in 1:nsteps
                step!(state, problem, ws, Δt)
            end
            σΔt = ν * n^2 * Δt
            factor = integrating_factor ? exp(-σΔt * nsteps) : stability_polynomial(r, -σΔt)^nsteps
            @test state.n == nsteps && state.t ≈ nsteps * Δt
            @test relative_error(state.fields.ux, T(factor) .* Array(initial.ux)) < 10 * nsteps * eps(T)
            # The FFT round-off of the projected advection reaches uy.
            @test sqrt(sum(abs2, Array(state.fields.uy)) / sum(abs2, Array(state.fields.ux))) < 10 * nsteps * eps(T)
            @test abs(kinetic_energy(state, problem, ws) / (E0 * factor^2) - 1) < 10 * nsteps * eps(T)
        end

        @testset "compressible linear decay: $(nameof(typeof(scheme))), IF $integrating_factor" for (scheme, r) in SCHEMES,
                                                                                                 integrating_factor in (false, true)
            # u = sin(2x) x̂ is compressive, u = sin(3y) x̂ solenoidal and ρ = 1 + 0.1 cos(4y); with nonlinear =
            # false they decay at (μ + ζ)k², μk² and κk², and the mean density stays 1.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            μ, ζ, κ, Δt, nsteps = 0.05, 0.03, 0.02, 0.05, 20
            ux, uy, ρ = (zeros(ComplexF64, 9, 16) for _ in 1:3)
            ux[3, 1] = -im / 2
            ux[1, 4], ux[1, 14] = -im / 2, im / 2
            ρ[1, 1] = 1
            ρ[1, 5] = ρ[1, 13] = 0.05
            device(a) = copyto!(KernelAbstractions.allocate(backend, ComplexF64, size(a)), a)
            state = State((ux = device(ux), uy = device(uy), rho = device(ρ)))
            ws = Workspace(g, state; scheme, integrating_factor)
            problem = Problem(g, Compressible(; eos = Isothermal(cs = 1.0), μ, ζ, κ, nonlinear = false), NoForcing())
            for _ in 1:nsteps
                step!(state, problem, ws, Δt)
            end
            decay(σ) = integrating_factor ? exp(-σ * Δt * nsteps) : stability_polynomial(r, -σ * Δt)^nsteps
            expected_ux = copy(ux)
            expected_ux[3, 1] *= decay((μ + ζ) * 4)
            expected_ux[1, [4, 14]] .*= decay(μ * 9)
            expected_ρ = copy(ρ)
            expected_ρ[1, [5, 13]] .*= decay(κ * 16)
            @test relative_error(state.fields.ux, expected_ux) < 100 * eps()
            @test relative_error(state.fields.rho, expected_ρ) < 100 * eps()
            @test all(iszero, Array(state.fields.uy))
        end

        @testset "conservation order: $(nameof(typeof(scheme)))" for (scheme, r) in SCHEMES[[1, 2, 4]]
            # Unforced inviscid run: the errors of E_kin and Z at t = 1 fall at least as fast as Δt^r.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, backend)
            problem = Problem(g, Incompressible(ν = 0.0), NoForcing())
            initial = random_velocity(g; spectrum = k -> k^4 * exp(-k^2 / 8), E_kin = 0.5, seed = 3)
            Δts = r == 1 ? [1 / 400, 1 / 800, 1 / 1600, 1 / 3200] : [1 / 25, 1 / 50, 1 / 100, 1 / 200]
            errors = map(Δts) do Δt
                state = State(map(copy, initial))
                ws = Workspace(g, state; scheme, integrating_factor = false)
                E0, Z0 = kinetic_energy(state, problem, ws), enstrophy(state, problem, ws)
                for _ in 1:round(Int, 1 / Δt)
                    step!(state, problem, ws, Δt)
                end
                (abs(kinetic_energy(state, problem, ws) / E0 - 1), abs(enstrophy(state, problem, ws) / Z0 - 1))
            end
            for errs in (first.(errors), last.(errors))
                @test minimum(errs) > 1e-12
                @test loglog_slope(Δts, errs) > r - 0.2
            end
        end

        @testset "injection rate of a forcing-only run: $(nameof(typeof(scheme))), IF $integrating_factor, $T" for (scheme, r) in SCHEMES,
                                                                                                             integrating_factor in (false, true),
                                                                                                             T in (Float64, Float32)
            # With forcing only, every scheme gives u_{n+1} = u_n + Δt f_n, so the stored rates equal the
            # change of E_kin and Z over the step up to the round-off of the sums.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            Δt = 1e-3
            forcing = OUForcing(g, band(4, 6); ε = 1.0, τ_f = 0.05, seed = 1)
            problem = Problem(g, Incompressible(ν = 0.0, nonlinear = false), forcing)
            state = State(random_velocity(g; spectrum = k -> k^4 * exp(-k^2 / 8), E_kin = 0.5, seed = 2))
            ws = Workspace(g, state; scheme, integrating_factor)
            E, Z = kinetic_energy(state, problem, ws), enstrophy(state, problem, ws)
            worst_E, worst_Z = zero(T), zero(T)
            for _ in 1:1000
                step!(state, problem, ws, Δt)
                E1, Z1 = kinetic_energy(state, problem, ws), enstrophy(state, problem, ws)
                inj = ws.diagnostics.injection
                worst_E = max(worst_E, abs(inj[:inj_E_kin] - (E1 - E) / T(Δt)) / (eps(T) * E1 / T(Δt)))
                worst_Z = max(worst_Z, abs(inj[:inj_Z] - (Z1 - Z) / T(Δt)) / (eps(T) * Z1 / T(Δt)))
                E, Z = E1, Z1
            end
            @test state.n == 1000
            @test worst_E ≤ 10
            @test worst_Z ≤ 10
        end
    end
end

end
