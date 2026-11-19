# Non-blocking tests of the compressible fluxes (WP-11, ADR-0011): Pi_kin over all grid modes is zero to
# round-off (Lamb form, specification note 2), and the aliasing floor of Pi_Omega at k_max drops as N grows.
# The file runs on the CPU and prints the measured values.
module CompressibleFluxesNonblocking

using Test
using CompReS

spectral_flux(name, state, problem, ws) = only(select_diagnostics(:SPECTRAL, (name,))).second(state, problem, ws)

# Compressible state of the physical functions on an N × N grid of [0, 2π)², with its workspace and problem.
function setup(N, ux, uy, rho)
    g = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N)
    x = [2π * (i - 1) / N for i in 1:N]
    ws0 = Workspace(g, State(kolmogorov(g; A = 1, n = 1)))
    spectral(f) = to_spectral!(similar(ws0.forcing.ux), f.(x, x'), ws0.transforms)
    state = State((ux = spectral(ux), uy = spectral(uy), rho = spectral(rho)))
    return g, state, Workspace(g, state), Problem(g, Compressible(; eos = Isothermal(cs = 3), μ = 0.1, ζ = 0), NoForcing())
end

@testset "non-blocking: compressible fluxes" begin
    @testset "Pi_kin over all grid modes is zero to round-off" begin
        # A random solenoidal field plus a compressive part, and a nonuniform density, all truncated.
        g = Grid(; Lx = 2π, Ly = 2π, Nx = 64, Ny = 64)
        u = random_velocity(g; spectrum = k -> k^2 * exp(-k / 3), E_kin = 1, seed = 31)
        r = random_velocity(g; spectrum = k -> k^2 * exp(-k / 3), E_kin = 1, seed = 32)
        fields = (ux = u.ux .+ 0.3 .* r.uy, uy = u.uy, rho = uniform_density(g).rho .+ 0.1 .* r.ux)
        state = State(fields)
        ws = Workspace(g, state)
        problem = Problem(g, Compressible(; eos = Isothermal(cs = 3), μ = 0.1, ζ = 0), NoForcing())
        Pi = spectral_flux(:Pi_kin, state, problem, ws)
        relative = abs(Pi[end]) / maximum(abs, Pi)
        println("Pi_kin over all grid modes, relative to max |Pi_kin|, N = 64: ", relative)
        @test relative < 1e-13
    end

    @testset "Pi_Omega(k_max) drops as N grows" begin
        # Fixed smooth fields: the density factors 1/√ρ and ρ^{-3/2} are not band-limited, so Pi_Omega does not
        # vanish at k_max; the floor falls with the resolution.
        ux(x, y) = 0.7sin(3y) + 0.5sin(2x) + 0.3cos(x + y)
        uy(x, y) = 0.4cos(2x - y) + 0.6sin(x)
        rho(x, y) = 1 + 0.3cos(x + 2y) + 0.2sin(2x)
        floors = map((16, 32, 64, 128)) do N
            g, state, ws, problem = setup(N, ux, uy, rho)
            Pi = spectral_flux(:Pi_Omega, state, problem, ws)
            n_max = findfirst(>(g.k_max), g.shell_edge)
            floor = abs(Pi[n_max]) / maximum(abs, Pi)
            println("Pi_Omega(k_max) relative to max |Pi_Omega|, N = $N: ", floor)
            floor
        end
        @test issorted(floors; rev = true)
        @test floors[end] < 1e-6
    end
end

end
