# Tests of the compressible right-hand side (WP-10): a state at rest, an acoustic wave, mass conservation,
# the incompressible limit, and the explicit viscous remainder against a hand computation.
module Nonlinear_compressibleTests

using Test
using KernelAbstractions
using CompReS

const NAMES = (:ux, :uy, :rho)

# Array position of the mode with integer indices (i, j) in the real-FFT layout.
pos(g, i, j) = (i + 1, j < 0 ? g.Ny + j + 1 : j + 1)

# Device spectral field of the physical field f(x, y) on the grid points of g.
function spectral(g, ws, f)
    T = eltype(g.k2)
    x = [T(g.Lx) * (i - 1) / g.Nx for i in 1:g.Nx]
    y = [T(g.Ly) * (j - 1) / g.Ny for j in 1:g.Ny]
    return to_spectral!(KernelAbstractions.allocate(g.backend, Complex{T}, g.Nx ÷ 2 + 1, g.Ny),
                        copyto!(similar(ws.physical.ux), f.(x, y')), ws.transforms)
end

# Compressible fields (ux, uy, rho) of the physical functions ux(x, y), uy(x, y) and rho(x, y).
function physical_fields(g, ux, uy, rho)
    ws = Workspace(g, State(kolmogorov(g; A = 1, n = 1)))
    return (ux = spectral(g, ws, ux), uy = spectral(g, ws, uy), rho = spectral(g, ws, rho))
end

# Tendency fields filled with NaN, so a test fails when rhs! leaves a value unwritten.
nan_fields(g, names = NAMES) =
    map(_ -> fill!(KernelAbstractions.allocate(g.backend, Complex{eltype(g.k2)}, g.Nx ÷ 2 + 1, g.Ny), NaN),
        NamedTuple{names}(names))

# The tendency rhs! of the fields on g, for the model, as host arrays.
function tendency(g, model, fields)
    state = State(fields)
    dû = nan_fields(g, keys(fields))
    rhs!(dû, state, Problem(g, model, NoForcing()), Workspace(g, state))
    return map(Array, dû)
end

# A random solenoidal velocity and the density 1 + a·r, with r a random truncated field of rms about 1.
function random_fields(g, a, seed)
    u = random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed)
    r = random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed = seed + 1).ux
    rho = uniform_density(g).rho .+ eltype(g.k2)(a) .* r
    return merge(u, (; rho))
end

const VISCOUS = (μ = 0.1, ζ = 0.05)

function run_tests(backend)
    @testset "compressible nonlinear term" begin
        @testset "uniform state at rest stays at rest, $(nameof(typeof(eos)))" for eos in
                (Isothermal(cs = 2), Polytropic(K = 1, γ = 1.4))
            # u = 0 and ρ = 1 cancel every term; h(1) ≠ 0 for Polytropic, but a constant has no gradient.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            fields = merge(map(zero, kolmogorov(g; A = 1, n = 1)), uniform_density(g))
            du = tendency(g, Compressible(; eos, VISCOUS...), fields)
            @test all(f -> all(iszero, f), du)
        end

        @testset "shear flow sin(2y) x̂ with uniform ρ stays steady" begin
            # u·∇u = 0 and ∇·(ρu) = 0: the Lamb vector and ∇(½|u|²) cancel, so a wrong sign of either fails.
            # Their mode (0, 4) lies inside the Galerkin truncation of N = 16 (k_c = 16/3).
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            du = tendency(g, Compressible(; eos = Isothermal(cs = 2), VISCOUS...),
                          merge(kolmogorov(g; A = 1, n = 2), uniform_density(g)))
            @test all(f -> maximum(abs, f) < 1e-14, du)
        end

        @testset "acoustic wave has frequency c_s k, $(nameof(typeof(eos)))" for (eos, cs) in
                ((Isothermal(cs = 2), 2.0), (Polytropic(K = 1, γ = 1.4), sqrt(1.4)))
            # The right-moving wave ρ = 1 + ε cos 2x, u_x = c_s ε cos 2x solves the linearised system
            # dρ̂/dt = −ik·û, dû/dt = −ik c_s² ρ̂, so each amplitude has the tendency −i c_s k times itself.
            # The nonlinear terms add O(ε²) at other modes and O(ε³) at the wave mode: the error must
            # fall by 4 when ε halves, and a wrong sign of the continuity flux or of ∇h leaves an O(ε) error.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            model = Compressible(; eos, VISCOUS...)
            k, p = 2, CartesianIndex(pos(g, 2, 0))
            function linear_error(ε)
                du = tendency(g, model, physical_fields(g, (x, y) -> cs * ε * cos(k * x), (x, y) -> 0.0,
                                                        (x, y) -> 1 + ε * cos(k * x)))
                ρ̂, ûx = ε / 2, cs * ε / 2
                @test abs(du.rho[p] - (-im * cs * k * ρ̂)) < 10ε^3
                @test abs(du.ux[p] - (-im * cs * k * ûx)) < 10ε^3
                expected = map(zero, du)
                expected.rho[p], expected.ux[p] = -im * cs * k * ρ̂, -im * cs * k * ûx
                return maximum(n -> maximum(abs, du[n] - expected[n]), NAMES)
            end
            e1, e2 = linear_error(1e-3), linear_error(5e-4)
            @test 3.5 < e1 / e2 < 4.5
        end

        @testset "mass has zero tendency, $T" for T in (Float64, Float32)
            # The continuity tendency is −ik·(ρu)^, so the k = 0 mode, the mass ⟨ρ⟩, gets exactly zero.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            du = tendency(g, Compressible(; eos = Isothermal(cs = 2), VISCOUS...), random_fields(g, 0.05, 1))
            @test iszero(du.rho[1, 1])
            @test maximum(abs, du.rho) > 0.01
        end

        @testset "uniform ρ, solenoidal u: solenoidal part equals the incompressible rhs!, $T" for T in
                (Float64, Float32)
            # At ρ = 1 the pressure and ½|u|² terms are gradients, the remainder (1/ρ − 1)(…) is zero and
            # ∇·u = 0 holds; P(−u·∇u) is then the incompressible tendency, and the continuity flux vanishes.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            u = random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed = 2)
            state = State(merge(u, uniform_density(g)))
            ws = Workspace(g, state)
            dû = nan_fields(g)
            rhs!(dû, state, Problem(g, Compressible(; eos = Isothermal(cs = 2), VISCOUS...), NoForcing()), ws)
            project_solenoidal!((dû.ux, dû.uy), ws.transforms)
            expected = tendency(g, Incompressible(ν = VISCOUS.μ), u)
            scale = maximum(abs, expected.ux)
            @test scale > 0.1
            @test maximum(abs, Array(dû.ux) - expected.ux) < 100eps(T) * scale
            @test maximum(abs, Array(dû.uy) - expected.uy) < 100eps(T) * scale
            @test maximum(abs, Array(dû.rho)) < 100eps(T) * scale
        end

        @testset "viscous remainder is zero at ρ = 1" begin
            # (1/ρ − 1) is exactly zero for uniform ρ, so μ and ζ leave the tendency unchanged, bit for bit.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            fields = physical_fields(g, (x, y) -> sin(y) + cos(2x), (x, y) -> sin(y), (x, y) -> 1.0)
            eos = Isothermal(cs = 2)
            @test tendency(g, Compressible(; eos, VISCOUS...), fields) ==
                  tendency(g, Compressible(; eos, μ = 0, ζ = 0), fields)
        end

        @testset "viscous remainder matches the hand computation for ρ = 1 + ε cos x" begin
            # u = (sin y + cos 2x, sin y): μ∇²u + ζ∇θ = (−μ sin y − 4(μ + ζ) cos 2x, −(μ + ζ) sin y). With
            # 1/(1 + ε cos x) − 1 = Σₙ aₙ e^{inx}, aₙ = c₀(−β)^|n| − δₙ₀, c₀ = 1/√(1 − ε²), β = (1 − √(1 − ε²))/ε,
            # the remainder (1/ρ − 1)(μ∇²u + ζ∇θ) has the modes (n, ±1): ±(iμ aₙ/2, i(μ + ζ) aₙ/2) and
            # (n, 0): (−2(μ + ζ)(aₙ₋₂ + aₙ₊₂), 0). The difference of the tendencies with and without
            # viscosity isolates it; a wrong sign of the remainder or of a viscous term fails. The grid aliases
            # aₙ₋₃₂ onto aₙ; ε = 0.3 (β ≈ 0.15) keeps this error near round-off.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, backend)
            ε = 0.3
            fields = physical_fields(g, (x, y) -> sin(y) + cos(2x), (x, y) -> sin(y), (x, y) -> 1 + ε * cos(x))
            eos = Isothermal(cs = 2)
            (; μ, ζ) = VISCOUS
            du = tendency(g, Compressible(; eos, μ, ζ), fields)
            du0 = tendency(g, Compressible(; eos, μ = 0, ζ = 0), fields)
            c₀, β = 1 / sqrt(1 - ε^2), (1 - sqrt(1 - ε^2)) / ε
            a(n) = c₀ * (-β)^abs(n) - (n == 0)
            expected_x, expected_y = zero(du.ux), zero(du.uy)
            keep = Array(g.truncation)
            for n in 0:(g.Nx ÷ 2 - 1)
                for j in (1, -1)
                    p = pos(g, n, j)
                    keep[p...] || continue
                    expected_x[p...] = j * im * μ * a(n) / 2
                    expected_y[p...] = j * im * (μ + ζ) * a(n) / 2
                end
                keep[pos(g, n, 0)...] && (expected_x[pos(g, n, 0)...] = -2 * (μ + ζ) * (a(n - 2) + a(n + 2)))
            end
            @test abs(expected_x[pos(g, 1, 1)...]) > 0.005 && abs(expected_x[pos(g, 1, 0)...]) > 0.005
            @test maximum(abs, du.ux - du0.ux - expected_x) < 1e-13
            @test maximum(abs, du.uy - du0.uy - expected_y) < 1e-13
            @test du.rho == du0.rho
        end

        @testset "tendency lies inside the Galerkin truncation" begin
            # Every product is truncated after the forward FFT, so modes outside the disc stay zero.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, backend)
            du = tendency(g, Compressible(; eos = Polytropic(K = 1, γ = 1.4), VISCOUS...), random_fields(g, 0.05, 3))
            keep = Array(g.truncation)
            @test all(f -> all(iszero, f[.!keep]), du)
            @test all(f -> maximum(abs, f) > 0.01, du)
        end

        @testset "nonlinear = false gives zero without forcing" begin
            # The test switch drops advection, pressure, remainder and continuity flux, so ρ stays fixed.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            du = tendency(g, Compressible(; eos = Isothermal(cs = 2), VISCOUS..., nonlinear = false),
                          random_fields(g, 0.05, 4))
            @test all(f -> all(iszero, f), du)
        end

        @testset "rhs! is type-stable and leaves the state unchanged" begin
            # ADR-0006: a Float32 run must infer and keep Float32; the state is an input only.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T = Float32, backend)
            fields = random_fields(g, 0.05f0, 5)
            before = map(Array, fields)
            state = State(fields)
            ws = Workspace(g, state)
            dû = nan_fields(g)
            problem = Problem(g, Compressible(; eos = Polytropic(K = 1, γ = 1.4), VISCOUS...), NoForcing())
            @test (@inferred rhs!(dû, state, problem, ws)) === dû
            @test eltype(dû.rho) == ComplexF32 && all(f -> all(isfinite, Array(f)), dû)
            @test map(Array, state.fields) == before
        end
    end
end

end
