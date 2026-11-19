# Tests of the initial conditions (WP-4): single Fourier modes on N = 8 or 16, and the random velocity field.
module Initial_conditionsTests

using Test
using KernelAbstractions
using CompReS

# Array position of the mode with integer indices (i, j) in the real-FFT layout.
pos(g, i, j) = (i + 1, j < 0 ? g.Ny + j + 1 : j + 1)

# Physical grid points x and y of the grid, as host arrays of size (Nx, 1) and (1, Ny).
coordinates(g) = ((0:(g.Nx - 1)) .* (g.Lx / g.Nx), ((0:(g.Ny - 1)) .* (g.Ly / g.Ny))')

# Half-plane weights of the real FFT: 1 on the kx = 0 and Nyquist columns, 2 elsewhere.
half_plane_weights(g) = [c == 1 || c == g.Nx ÷ 2 + 1 ? 1 : 2 for c in 1:(g.Nx ÷ 2 + 1), _ in 1:g.Ny]

# Physical host arrays of the state fields, through the physical buffers of the Workspace.
function physical(g, fields)
    ws = Workspace(g, State(fields))
    return map((f, f̂) -> Array(to_physical!(f, f̂, ws.transforms)), ws.physical, fields)
end

# Spectral kinetic energy ½ Σ w |û|² over the stored half-plane.
spectral_E_kin(g, u) = sum(half_plane_weights(g) .* (abs2.(Array(u.ux)) .+ abs2.(Array(u.uy)))) / 2

# Shell spectrum ½ Σ_{q ∈ shell n} |û(q)|² over the full plane, entry n + 1 for shell n.
function shell_spectrum(g, u)
    e = half_plane_weights(g) .* (abs2.(Array(u.ux)) .+ abs2.(Array(u.uy))) ./ 2
    s = zeros(length(g.shell_k))
    for (n, v) in zip(Array(g.shell_index), e)
        s[n] += v
    end
    return s
end

function run_tests(backend)
    @testset "initial conditions" begin
        @testset "kolmogorov" begin
            # A sin(y) = A (e^{iy} - e^{-iy})/(2i): with A = 2, ûx(0, 1) = -i and ûx(0, -1) = i; every other mode is zero.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, backend)
            u = kolmogorov(g; A = 2, n = 1)
            @test keys(u) == (:ux, :uy)
            expected = zeros(ComplexF64, 5, 8)
            expected[pos(g, 0, 1)...] = -im
            expected[pos(g, 0, -1)...] = im
            @test Array(u.ux) == expected
            @test all(iszero, Array(u.uy))
            # On Ly = 4π the mode index n = 2 is the wavenumber k_n = 2 · 2π/4π = 1, so ux = A sin(y).
            g = Grid(; Lx = 4π, Ly = 4π, Nx = 16, Ny = 16, backend)
            _, y = coordinates(g)
            p = physical(g, kolmogorov(g; A = 3, n = 2))
            @test maximum(abs, p.ux .- 3 .* sin.(y)) <= 10 * eps(3.0)
            @test all(iszero, p.uy)
            # On N = 8 the truncation keeps |k| < 8/3: n = 2 is kept, n = 3 is not; n = 0 gives no flow.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, backend)
            @test kolmogorov(g; A = 1, n = 2) isa NamedTuple
            @test_throws ArgumentError kolmogorov(g; A = 1, n = 3)
            @test_throws ArgumentError kolmogorov(g; A = 1, n = 0)
        end

        @testset "crossed_kolmogorov" begin
            # On Lx = 4π, Ly = 2π, mode indices k1 = 1, k2 = 2 are wavenumbers 1 in y and 2 · 2π/4π = 1 in x:
            # u = (sin y, 3 sin x), so ûx(0, ±1) = ∓i/2 and ûy(2, 0) = -3i/2 on the stored half-plane.
            g = Grid(; Lx = 4π, Ly = 2π, Nx = 16, Ny = 8, backend)
            u = crossed_kolmogorov(g; A1 = 1, k1 = 1, A2 = 3, k2 = 2)
            ux, uy = zeros(ComplexF64, 9, 8), zeros(ComplexF64, 9, 8)
            ux[pos(g, 0, 1)...] = -im / 2
            ux[pos(g, 0, -1)...] = im / 2
            uy[pos(g, 2, 0)...] = -3im / 2
            @test Array(u.ux) == ux && Array(u.uy) == uy
            x, y = coordinates(g)
            p = physical(g, u)
            @test maximum(abs, p.ux .- sin.(y)) <= 10 * eps()
            @test maximum(abs, p.uy .- 3 .* sin.(x)) <= 10 * eps(3.0)
            # Each mode must lie inside the truncation: k2 = 6 is |k| = 6 · 1/2 = 3, above k_c = (16/3)(1/2) = 8/3.
            @test_throws ArgumentError crossed_kolmogorov(g; A1 = 1, k1 = 1, A2 = 1, k2 = 6)
        end

        @testset "uniform_density" begin
            # ρ = 1 has the single coefficient ρ̂(0, 0) = 1.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, backend)
            ρ = uniform_density(g)
            @test keys(ρ) == (:rho,)
            expected = zeros(ComplexF64, 5, 8)
            expected[1, 1] = 1
            @test Array(ρ.rho) == expected
            # The velocity and the density merge into the fields of a compressible State.
            state = State(merge(kolmogorov(g; A = 1, n = 1), ρ))
            @test keys(state.fields) == (:ux, :uy, :rho)
            @test all(isone, physical(g, state.fields).rho)
        end

        @testset "random_velocity" begin
            spectrum(k) = k^4 * exp(-k^2 / 8)
            for g in (Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend), Grid(; Lx = 4π, Ly = 2π, Nx = 16, Ny = 8, backend))
                u = random_velocity(g; spectrum, E_kin = 0.5, seed = 2)
                ux, uy = Array(u.ux), Array(u.uy)
                kx, ky, keep = Array(g.kx), Array(g.ky), Array(g.truncation)
                # E_kin = ½ Σ w |û|² over the half-plane equals the target.
                @test spectral_E_kin(g, u) ≈ 0.5 rtol = 10 * eps()
                # By Parseval the grid mean ½⟨|u|²⟩ of the physical field equals the same target.
                p = physical(g, u)
                @test sum(abs2, p.ux) / 2length(p.ux) + sum(abs2, p.uy) / 2length(p.uy) ≈ 0.5 rtol = 100 * eps()
                # Solenoidal: k · û = 0 on every mode, to round-off.
                @test maximum(abs, kx .* ux .+ ky .* uy) <= 10 * eps() * maximum(sqrt.(kx .^ 2 .+ ky .^ 2) .* abs.(ux))
                # Inside the truncation: every mode outside the disc is exactly zero.
                @test all(iszero, ux[.!keep]) && all(iszero, uy[.!keep])
                # Hermitian on the kx = 0 column: û(0, -ky) == conj(û(0, ky)) exactly, and the mean flow is zero.
                for j in 1:(g.Ny ÷ 2 - 1)
                    @test ux[pos(g, 0, -j)...] == conj(ux[pos(g, 0, j)...])
                    @test uy[pos(g, 0, -j)...] == conj(uy[pos(g, 0, j)...])
                end
                @test ux[1, 1] == 0 && uy[1, 1] == 0
                # The shell spectrum is proportional to spectrum(k_n): the ratio is the same in every shell with a mode.
                s = shell_spectrum(g, u)
                shells = findall(>(0), g.shell_count)[2:end]
                ratio = s[shells] ./ spectrum.(g.shell_k[shells])
                @test all(≈(ratio[1]; rtol = 100 * eps()), ratio)
            end
        end

        @testset "random_velocity in one shell" begin
            # On a square 2π box shell 1 holds the 8 modes (±1, 0), (0, ±1), (±1, ±1); with all energy in shell 1,
            # each has |û|² = 2 E_kin / 8 = 1/4 for E_kin = 1, and every other mode is zero.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, backend)
            u = random_velocity(g; spectrum = k -> k == 1 ? 1.0 : 0.0, E_kin = 1, seed = 3)
            a2 = abs2.(Array(u.ux)) .+ abs2.(Array(u.uy))
            in_shell = Array(g.shell_index) .== 2
            @test g.shell_count[2] == 8
            @test all(≈(1 / 4; rtol = 10 * eps()), a2[in_shell])
            @test all(iszero, a2[.!in_shell])
        end

        @testset "random_velocity seeds" begin
            # The host generator makes the field a function of the seed only: the backend gives the CPU field bit for bit.
            on(b) = Grid(; Lx = 4π, Ly = 2π, Nx = 16, Ny = 8, backend = b)
            spectrum(k) = k^2 * exp(-k)
            u = random_velocity(on(backend); spectrum, E_kin = 0.5, seed = 4)
            u_cpu = random_velocity(on(CPU()); spectrum, E_kin = 0.5, seed = 4)
            @test Array(u.ux) == Array(u_cpu.ux) && Array(u.uy) == Array(u_cpu.uy)
            # Another seed gives other phases.
            @test Array(random_velocity(on(backend); spectrum, E_kin = 0.5, seed = 5).ux) != Array(u.ux)
            # A spectrum that is zero on every retained mode cannot reach a nonzero E_kin.
            @test_throws ArgumentError random_velocity(on(backend); spectrum = k -> 0.0, E_kin = 1, seed = 1)
        end

        @testset "Float32 stays Float32" begin
            # Every initial condition takes the float type of the grid, even from Float64 arguments and spectrum values.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float32, backend)
            u = @inferred random_velocity(g; spectrum = k -> 0.5 * k^2 * exp(-k^2), E_kin = 0.5, seed = 1)
            @test eltype(u.ux) == ComplexF32 && eltype(u.uy) == ComplexF32
            @test spectral_E_kin(g, u) ≈ 0.5f0 rtol = 10 * eps(Float32)
            @test eltype((@inferred kolmogorov(g; A = 1.0, n = 1)).ux) == ComplexF32
            @test eltype((@inferred crossed_kolmogorov(g; A1 = 1.0, k1 = 1, A2 = 1.0, k2 = 2)).uy) == ComplexF32
            @test eltype((@inferred uniform_density(g)).rho) == ComplexF32
        end
    end
end

end
