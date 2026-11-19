# Tests of the incompressible nonlinear term, Problem and rhs! (WP-6): exact steady flows, a hand-computed
# two-mode interaction, and the vanishing energy and enstrophy fluxes at the Galerkin truncation.
module Nonlinear_incompressibleTests

using Test
using KernelAbstractions
using CompReS

const NAMES = (:ux, :uy)

# Spectral fields of the physical velocity (fx(x, y), fy(x, y)) on the grid points of g.
function spectral_velocity(g, ws, fx, fy)
    T = eltype(g.k2)
    x = [T(g.Lx) * (i - 1) / g.Nx for i in 1:g.Nx]
    y = [T(g.Ly) * (j - 1) / g.Ny for j in 1:g.Ny]
    spectral(f) = to_spectral!(KernelAbstractions.allocate(g.backend, Complex{T}, g.Nx ÷ 2 + 1, g.Ny),
                               copyto!(similar(ws.physical.ux), f.(x, y')), ws.transforms)
    return (ux = spectral(fx), uy = spectral(fy))
end

# Tendency fields filled with NaN, so a test fails when rhs! leaves a value unwritten.
nan_fields(g) = map(_ -> fill!(KernelAbstractions.allocate(g.backend, Complex{eltype(g.k2)}, g.Nx ÷ 2 + 1, g.Ny),
                               NaN), NamedTuple{NAMES}(NAMES))

# The tendency rhs! of the velocity u on g, for the model, as host arrays.
function tendency(g, model, u)
    state = State(u)
    ws = Workspace(g, state)
    dû = nan_fields(g)
    rhs!(dû, state, Problem(g, model, NoForcing()), ws)
    return map(Array, dû)
end

# Half-plane weight of the stored column c: no Nyquist mode is retained, so only the kx = 0 column counts once.
weight(c) = c == 1 ? 1 : 2

# Weighted sums Σ w Re[a*·b] and Σ w |a||b| over the stored modes; a and b are pairs of host spectral arrays.
function weighted_products(a, b)
    s, scale = 0.0, 0.0
    for I in CartesianIndices(a[1])
        w = weight(I[1])
        s += w * sum(real(conj(a[n][I]) * b[n][I]) for n in 1:length(a))
        scale += w * sum(abs(a[n][I]) * abs(b[n][I]) for n in 1:length(a))
    end
    return s, scale
end

# Spectral vorticity i(kx ûy − ky ûx) of the host velocity (ux, uy) on g.
curl(g, ux, uy) = im .* (Array(g.kx) .* uy .- Array(g.ky) .* ux)

function run_tests(backend)
    @testset "incompressible nonlinear term" begin
        model = Incompressible(ν = 0.1)

        @testset "Problem bundles grid, model in the grid type, forcing and rates" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float32, backend)
            p = @inferred Problem(g, model, NoForcing())
            # The model takes the float type of the grid, so a Float32 run never promotes (ADR-0006).
            @test p.grid === g && p.model === convert_float(Float32, model) && p.forcing === NoForcing()
            @test Array(p.rates.velocity) == Array(linear_rates(g, model).velocity)
        end

        @testset "shear flow sin(ky) x̂ gives zero" begin
            # u·∇u = 0 for a parallel shear flow, so the advection vanishes at every mode.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            du = tendency(g, model, kolmogorov(g; A = 1.0, n = 3))
            @test maximum(abs, du.ux) < 1e-14 && maximum(abs, du.uy) < 1e-14
        end

        @testset "cellular flow ψ = sin x sin y gives zero" begin
            # ω = 2ψ is a function of ψ, so u·∇ω = 0, and u·∇u is a gradient that the projection removes.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            ws = Workspace(g, State(kolmogorov(g; A = 1.0, n = 1)))
            u = spectral_velocity(g, ws, (x, y) -> sin(x) * cos(y), (x, y) -> -cos(x) * sin(y))
            du = tendency(g, model, u)
            @test maximum(abs, du.ux) < 1e-14 && maximum(abs, du.uy) < 1e-14
        end

        @testset "two crossed modes give the hand-computed tendency" for T in (Float64, Float32)
            # u = (sin y, sin 2x): u·∇ω = −3 sin 2x sin y, so ∂ψ/∂t = (3/5) sin 2x sin y and
            # ∂u/∂t = (3/5 sin 2x cos y, −6/5 cos 2x sin y). A sign or projection error changes this value.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T, backend)
            ws = Workspace(g, State(kolmogorov(g; A = 1.0, n = 1)))
            u = spectral_velocity(g, ws, (x, y) -> sin(y), (x, y) -> sin(2x))
            expected = map(Array, spectral_velocity(g, ws, (x, y) -> 3sin(2x) * cos(y) / 5,
                                                    (x, y) -> -6cos(2x) * sin(y) / 5))
            du = tendency(g, model, u)
            @test eltype(du.ux) == Complex{T}
            @test maximum(abs, du.ux - expected.ux) < 10eps(T) && maximum(abs, du.uy - expected.uy) < 10eps(T)
            # The two coefficients of the (2, 1) mode, for a reader who checks by hand.
            @test du.ux[3, 2] ≈ -3im / 20 && du.uy[3, 2] ≈ 3im / 10
        end

        @testset "energy and enstrophy fluxes vanish on a random field" for T in (Float64, Float32)
            # Specification, Galerkin truncation: Σ Re[û*·N̂] = Σ Re[ω̂*(u·∇ω)^] = 0 over the retained modes.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 64, Ny = 64, T, backend)
            u = random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed = 3)
            du = tendency(g, model, u)
            û = map(Array, u)
            energy, energy_scale = weighted_products((û.ux, û.uy), (du.ux, du.uy))
            ω̂, dω̂ = curl(g, û.ux, û.uy), curl(g, du.ux, du.uy)
            enstrophy, enstrophy_scale = weighted_products((ω̂,), (dω̂,))
            @test energy_scale > 1 && enstrophy_scale > 1
            @test abs(energy) < 100eps(T) * energy_scale
            @test abs(enstrophy) < 100eps(T) * enstrophy_scale
            # The tendency is solenoidal and lies inside the Galerkin truncation.
            kx, ky, keep = Array(g.kx), Array(g.ky), Array(g.truncation)
            @test maximum(abs, kx .* du.ux .+ ky .* du.uy) < 100eps(T) * maximum(abs, du.ux)
            @test all(iszero, du.ux[.!keep]) && all(iszero, du.uy[.!keep])
        end

        @testset "nonlinear = false gives zero without forcing" begin
            # The test switch drops the advection; NoForcing adds nothing, so every mode of dû is zero.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            du = tendency(g, Incompressible(ν = 0.1, nonlinear = false),
                          random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed = 4))
            @test all(iszero, du.ux) && all(iszero, du.uy)
        end

        @testset "rhs! is type-stable and leaves the state unchanged" begin
            # ADR-0006: a Float32 run must infer and keep Float32; the state is an input only.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T = Float32, backend)
            u = random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed = 5)
            before = map(Array, u)
            state = State(u)
            ws = Workspace(g, state)
            dû = nan_fields(g)
            @test (@inferred rhs!(dû, state, Problem(g, model, NoForcing()), ws)) === dû
            @test map(Array, state.fields) == before
        end
    end
end

end
