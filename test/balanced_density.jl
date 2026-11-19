# Tests of the balanced-density routine (WP-14, ADR-0013) on the cellular flow u = (sin x cos y, −cos x sin y) at
# Ma = 0.1, with μ = ζ = κ = 0. Its incompressible pressure is h = (cos 2x + cos 2y)/4. The tests measure θ(0) and
# ∂_tθ(0) with the production rhs!, compare with the wrong procedure (shift ρ after the inversion), and run the
# production step! to compare the maximum of E_kin_comp with the uniform-density run.
module Balanced_densityTests

using Test
using KernelAbstractions
using CompReS

const Ma = 0.1
# u_rms of the cellular flow is 1/√2, so Ma = u_rms/c_s fixes c_s.
const cs = 1 / (sqrt(2) * Ma)
const EOS = (Isothermal(cs = cs), Polytropic(K = cs^2 / 1.4, γ = 1.4))

# Host grid coordinates of g, as Float64 matrices X[i, j] = x_i, Y[i, j] = y_j.
function coordinates(g)
    x = [g.Lx * (i - 1) / g.Nx for i in 1:g.Nx]
    y = [g.Ly * (j - 1) / g.Ny for j in 1:g.Ny]
    return x .+ 0 .* y', 0 .* x .+ y'
end

# Device spectral field of the host physical values f, truncated to the disc.
function spectral(g, tr, f)
    T = eltype(g.k2)
    f̂ = KernelAbstractions.allocate(g.backend, Complex{T}, g.Nx ÷ 2 + 1, g.Ny)
    p = copyto!(KernelAbstractions.allocate(g.backend, T, g.Nx, g.Ny), T.(f))
    return truncate!(to_spectral!(f̂, p, tr), g)
end

# Host physical field of the spectral field f̂.
physical(g, tr, f̂) = Array(to_physical!(KernelAbstractions.allocate(g.backend, eltype(g.k2), g.Nx, g.Ny), f̂, tr))

# Grid of [0, 2π)² with N × N points, its transforms, the cellular flow (; ux, uy) and the host coordinates.
function cellular(N, backend; T = Float64)
    g = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N, T, backend)
    tr = Transforms(g)
    X, Y = coordinates(g)
    u = (ux = spectral(g, tr, sin.(X) .* cos.(Y)), uy = spectral(g, tr, -cos.(X) .* sin.(Y)))
    return g, tr, u, X, Y
end

# Inviscid compressible problem with the equation of state `eos`.
inviscid(g, eos) = Problem(g, Compressible(; eos, μ = 0, ζ = 0, κ = 0), NoForcing())

# max |ik·∂_t û| over all modes, with ∂_t û from the production rhs! (the linear terms are zero for μ = ζ = κ = α = 0).
function divergence_tendency(u, rho, problem)
    state = State((ux = copy(u.ux), uy = copy(u.uy), rho = copy(rho)))
    ws = Workspace(problem.grid, state)
    dû = map(similar, state.fields)
    rhs!(dû, state, problem, ws)
    kx, ky = Array(ws.transforms.kx_deriv), Array(ws.transforms.ky_deriv)'
    return maximum(abs, im .* kx .* Array(dû.ux) .+ im .* ky .* Array(dû.uy))
end

# ∂_tθ(0) of the density rho, relative to its value for uniform ρ, max|k_ik_j(u_iu_j)^|.
relative_divergence_tendency(u, rho, problem) =
    divergence_tendency(u, rho, problem) / divergence_tendency(u, uniform_density(problem.grid).rho, problem)

# The wrong procedure: invert the enthalpy without the constant h₀, then shift ρ to unit mean.
function shifted_density(g, tr, eos, X, Y)
    h = (cos.(2X) .+ cos.(2Y)) ./ 4
    ρ = density_from_enthalpy.(Ref(eos), h .+ enthalpy(eos, 1.0))
    return spectral(g, tr, ρ .- sum(ρ) / length(ρ) .+ 1)
end

# Maximum of E_kin_comp over t ∈ [0, 10 max(Lx, Ly)/c_s], through the production step!, for the initial density rho.
function max_compressive_energy(u, rho, problem, Δt)
    g = problem.grid
    state = State((ux = copy(u.ux), uy = copy(u.uy), rho = copy(rho)))
    ws = Workspace(g, state; scheme = RK4(), integrating_factor = true)
    E_comp = only(select_diagnostics(:GLOBAL, (:E_kin_comp,))).second
    t_end = 10 * max(g.Lx, g.Ly) / cs
    E_max = E_comp(state, problem, ws)
    for _ in 1:ceil(Int, t_end / Δt)
        step!(state, problem, ws, Δt)
        E_max = max(E_max, E_comp(state, problem, ws))
    end
    return E_max
end

function run_tests(backend)
    @testset "balanced_density" begin
        @testset "$(nameof(typeof(eos))): mean density and enthalpy shape" for eos in EOS
            g, tr, u, X, Y = cellular(64, backend)
            (; rho) = balanced_density(u, eos, g)
            @test keys(balanced_density(u, eos, g)) == (:rho,)
            # The grid sum of the root find and the forward FFT each add round-off; 100 eps bounds both.
            @test abs(Array(rho)[1, 1] - 1) < 100eps()
            @test all(iszero, Array(rho)[.!Array(g.truncation)])
            # h(ρ) equals the pressure (cos 2x + cos 2y)/4 of the cellular flow, up to the constant h₀.
            dh = enthalpy.(Ref(eos), physical(g, tr, rho)) .- (cos.(2X) .+ cos.(2Y)) ./ 4
            @test maximum(abs, dh .- sum(dh) / length(dh)) < 1e-10
        end

        @testset "$(nameof(typeof(eos))): θ(0) = 0 and ∂_tθ(0) at N = 64" for eos in EOS
            g, tr, u, X, Y = cellular(64, backend)
            problem = inviscid(g, eos)
            ws = Workspace(g, State(merge(u, uniform_density(g))))
            @test CompReS.mach_number(State(merge(u, uniform_density(g))), problem, ws) ≈ Ma
            kx, ky = Array(ws.transforms.kx_deriv), Array(ws.transforms.ky_deriv)'
            @test maximum(abs, kx .* Array(u.ux) .+ ky .* Array(u.uy)) < 10eps()
            rho = balanced_density(u, eos, g).rho
            balanced = relative_divergence_tendency(u, rho, problem)
            wrong = relative_divergence_tendency(u, shifted_density(g, tr, eos, X, Y), problem)
            println("WP-14 $(nameof(typeof(eos))), N = 64, $(nameof(typeof(backend))): ∂_tθ(0) relative = ",
                    balanced, ", wrong procedure (shift ρ) = ", wrong)
            @test balanced < 1e-9
            @test wrong > 1e-6
        end

        @testset "Float32 grid with a Float64 equation of state" begin
            g, tr, u, X, Y = cellular(16, backend; T = Float32)
            rho = balanced_density(u, EOS[2], g).rho
            @test eltype(rho) == ComplexF32
            @test abs(Array(rho)[1, 1] - 1) < 10eps(Float32)
        end

        @testset "$(nameof(typeof(eos))): max E_kin_comp at least 10 times below uniform ρ" for eos in EOS
            N = 16
            g, tr, u, X, Y = cellular(N, backend)
            problem = inviscid(g, eos)
            Δt = 0.5 / (cs * g.k_max)
            t0 = time()
            balanced = max_compressive_energy(u, balanced_density(u, eos, g).rho, problem, Δt)
            uniform = max_compressive_energy(u, uniform_density(g).rho, problem, Δt)
            println("WP-14 $(nameof(typeof(eos))), N = $N, $(nameof(typeof(backend))): max E_kin_comp balanced = ",
                    balanced, ", uniform = ", uniform, ", ratio = ", uniform / balanced,
                    ", runtime = ", round(time() - t0; digits = 2), " s")
            @test uniform >= 10balanced
        end
    end
end

end
