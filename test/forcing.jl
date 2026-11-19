# Tests of the forcing (WP-9): supports counted by hand on small grids, the Kolmogorov formula, Hermitian
# symmetry, the random-phase reset times, the solenoidal output, restart, and the OU normalisation.
module ForcingTests

using Test
using KernelAbstractions
using CompReS
using CompReS: update_forcing!, add_forcing!, forcing_state, set_forcing_state!, mode_count

spectral_zeros(g) = KernelAbstractions.zeros(g.backend, Complex{eltype(g.k2)}, g.Nx ÷ 2 + 1, g.Ny)
zero_velocity(g) = (ux = spectral_zeros(g), uy = spectral_zeros(g))

# Array position of the mode with integer indices (i, j), i ≥ 0, in the real-FFT layout.
pos(g, i, j) = (i + 1, j < 0 ? g.Ny + j + 1 : j + 1)

# A State at rest and its Workspace on g.
function rest(g)
    state = State(zero_velocity(g))
    return state, Workspace(g, state)
end

# Advance the forcing by one step of length Δt, as the time stepping does: update, then advance t and n.
function advance!(forcing, state, ws, Δt)
    update_forcing!(forcing, state, ws, Δt)
    state.t += Δt
    state.n += 1
    return state
end

# The forcing f̂ of the current step as host arrays.
forcing_field(ws) = (Array(ws.forcing.ux), Array(ws.forcing.uy))

# Advance the forcing by one step and return the new f̂ as host arrays.
function next_field!(forcing, state, ws, Δt)
    advance!(forcing, state, ws, Δt)
    return forcing_field(ws)
end

# Positions where f̂ is nonzero, after one update of the forcing on g.
function forced_positions(g, forcing)
    state, ws = rest(g)
    update_forcing!(forcing, state, ws, 0.01)
    fx, fy = forcing_field(ws)
    return Set(Tuple(I) for I in CartesianIndices(fx) if fx[I] != 0 || fy[I] != 0)
end

ou(g, s; kw...) = OUForcing(g, s; ε = 1, τ_f = 1, seed = 1, kw...)

# Whether f̂ is Hermitian on the kx = 0 column: f̂(0, −ky) == conj(f̂(0, ky)), f̂(0, 0) real.
hermitian_column(g, f) = isreal(f[1, 1]) && all(f[1, g.Ny + 2 - r] == conj(f[1, r]) for r in 2:(g.Ny ÷ 2))

# Largest |k·f̂| relative to the largest |k||f̂|.
function divergence(g, fx, fy)
    kx, ky = Array(g.kx), Array(g.ky)
    return maximum(abs.(kx .* fx .+ ky .* fy)) / maximum(sqrt.(kx .^ 2 .+ ky .^ 2) .* hypot.(abs.(fx), abs.(fy)))
end

# Full-plane sum Σ |f̂(k)|² = ⟨|f|²⟩: weight 1 on the kx = 0 column, 2 elsewhere (no Nyquist mode is forced).
mean_square(fx, fy) = sum((c == 1 ? 1 : 2) * (abs2(fx[c, r]) + abs2(fy[c, r])) for c in axes(fx, 1), r in axes(fx, 2))

# OU normalisation (CPU only). Each of M realisations drives the forced modes with u' = f from u = 0
# and samples the injection ⟨u·f⟩ + ½Δt⟨|f|²⟩ at the start of a step once per τ_f over [20τ_f, 40τ_f].
# The test keeps v = u/Δt, so that add_forcing! gives v_{n+1} = v_n + f̂_n on the forced modes only.
# Returns the M realisation means.
function ou_realisation_means(; M, ε = 1.0, τ_f = 1.0)
    g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32)
    Δt = τ_f / 100
    state, ws = rest(g)
    v = state.fields
    fx, fy = similar(ws.physical.ux), similar(ws.physical.uy)
    vx, vy = ws.physical.ux, ws.physical.uy
    means = zeros(M)
    for r in 1:M
        forcing = OUForcing(g, band(4, 6); ε, τ_f, seed = r)
        foreach(a -> fill!(a, 0), v)
        state.t, state.n = 0, 0
        total, samples = 0.0, 0
        for n in 0:4000
            update_forcing!(forcing, state, ws, Δt)
            if n >= 2000 && n % 100 == 0
                to_physical!(vx, v.ux, ws.transforms)
                to_physical!(vy, v.uy, ws.transforms)
                to_physical!(fx, ws.forcing.ux, ws.transforms)
                to_physical!(fy, ws.forcing.uy, ws.transforms)
                total += Δt * sum((vx .+ fx ./ 2) .* fx .+ (vy .+ fy ./ 2) .* fy) / length(fx)
                samples += 1
            end
            add_forcing!(v, forcing, state, ws)
            state.t += Δt
            state.n += 1
        end
        means[r] = total / samples
    end
    return means
end

function run_tests(backend)
    @testset "forcing" begin
        @testset "supports counted by hand on N = 16" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            # |k|² ∈ [1, 2.25]: (±1, 0), (0, ±1) and (±1, ±1).
            @test mode_count(ou(g, band(1, 1.5))) == 8
            # |k| = 3: (±3, 0), (0, ±3).
            @test mode_count(ou(g, band(3, 3))) == 4
            # |k| = 5: (±5, 0), (0, ±5), (±3, ±4), (±4, ±3); 9·25 < 16², so all are retained.
            @test mode_count(ou(g, band(5, 5))) == 12
            # A listed mode stands for the pair ±k; a repeated pair counts once.
            @test mode_count(ou(g, modes([(1, 2)]))) == 2
            @test mode_count(ou(g, modes([(1, 2), (-1, -2)]))) == 2
            @test mode_count(ou(g, modes([(0, 3), (2, -1)]))) == 4
            @test mode_count(ou(g, support((kx, ky) -> kx == 0 && 0 < abs(ky) <= 2))) == 4
            # A support holds −k with every k: (1, 2) alone selects the pair ±(1, 2).
            @test mode_count(ou(g, support((kx, ky) -> (kx, ky) == (1, 2)))) == 2
            # The forcing acts on the stored half of the pairs: (−1, −2) is stored as (1, 2); on the
            # kx = 0 column both (0, 3) and (0, −3) are stored.
            @test forced_positions(g, ou(g, modes([(-1, -2)]))) == Set([pos(g, 1, 2)])
            @test forced_positions(g, ou(g, modes([(0, 3)]))) == Set([pos(g, 0, 3), pos(g, 0, -3)])
            @test forced_positions(g, ou(g, band(1, 1.5))) ==
                  Set([pos(g, 1, 0), pos(g, 1, 1), pos(g, 1, -1), pos(g, 0, 1), pos(g, 0, -1)])
        end

        @testset "band selects by physical |k| on Lx = 4π, Ly = 2π" begin
            g = Grid(; Lx = 4π, Ly = 2π, Nx = 32, Ny = 16, backend)
            # kx = i/2, ky = j: |k| ∈ [0.5, 1] holds the indices (1, 0), (2, 0) and (0, ±1), but not
            # (1, ±1), whose |k| = √1.25.
            f = ou(g, band(0.5, 1))
            @test mode_count(f) == 6
            @test forced_positions(g, f) == Set([pos(g, 1, 0), pos(g, 2, 0), pos(g, 0, 1), pos(g, 0, -1)])
        end

        @testset "empty or out-of-truncation support rejected" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            # No |k|² lies in [1.21, 1.69].
            @test_throws ArgumentError ou(g, band(1.1, 1.3))
            # (2, 5) has |k|² = 29 and 9·29 > 16²: outside the Galerkin truncation.
            @test_throws ArgumentError ou(g, band(5, 6))
            @test_throws ArgumentError ou(g, modes([(6, 0)]))
            @test_throws ArgumentError ou(g, modes(Tuple{Int,Int}[]))
            # The mode k = 0 is not forced.
            @test_throws ArgumentError ou(g, band(0, 1))
            @test_throws ArgumentError ou(g, modes([(0, 0)]))
            @test_throws ArgumentError RandomPhaseForcing(g, band(5, 6); f_rms = 1, τ_r = 1, seed = 1)
            @test_throws ArgumentError KolmogorovForcing(g; F = 1, n = 6)
            @test_throws ArgumentError KolmogorovForcing(g; F = 1, n = 0)
            @test_throws ArgumentError band(2, 1)
        end

        @testset "Kolmogorov field matches the formula" begin
            # F sin(3y) = F (e^{3iy} − e^{−3iy})/(2i): with F = 2, f̂x(0, 3) = −i and f̂x(0, −3) = i.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            f = KolmogorovForcing(g; F = 2, n = 3)
            @test mode_count(f) == 2
            state, ws = rest(g)
            update_forcing!(f, state, ws, 0.01)
            fx, fy = forcing_field(ws)
            expected = zeros(ComplexF64, size(fx))
            expected[pos(g, 0, 3)...] = -im
            expected[pos(g, 0, -3)...] = im
            @test fx == expected
            @test all(iszero, fy)
            # The physical field is F sin(3y) at the grid points.
            y = (0:15) .* (2π / 16)
            physical = Array(to_physical!(ws.physical.ux, ws.forcing.ux, ws.transforms))
            @test physical ≈ repeat(2 .* sin.(3 .* y'), 16) atol = 1e-13
        end

        @testset "every generated f̂ is Hermitian on the kx = 0 column" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            for f in (ou(g, band(1, 3)), ou(g, band(1, 3); compressive_fraction = 0.5),
                      RandomPhaseForcing(g, band(1, 3); f_rms = 1, τ_r = 0.02, seed = 3),
                      RandomPhaseForcing(g, band(1, 3); f_rms = 1, τ_r = 0.02, seed = 3, compressive_fraction = 1),
                      KolmogorovForcing(g; F = 1, n = 2))
                state, ws = rest(g)
                fields = [next_field!(f, state, ws, 0.01) for _ in 1:5]
                @test all(((fx, fy),) -> hermitian_column(g, fx) && hermitian_column(g, fy), fields)
                # The kx = 0 column is forced, so the check is not empty.
                @test all(((fx, fy),) -> fx[pos(g, 0, 2)...] != 0 || fy[pos(g, 0, 2)...] != 0, fields)
            end
        end

        @testset "phases change only at multiples of τ_r" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            f = RandomPhaseForcing(g, band(1, 3); f_rms = 1, τ_r = 0.05, seed = 4)
            state, ws = rest(g)
            previous = nothing
            changed = Int[]
            for n in 0:30
                update_forcing!(f, state, ws, 0.01)
                current = forcing_field(ws)
                current == previous || push!(changed, n)
                previous = current
                state.t += 0.01
                state.n += 1
            end
            # τ_r = 5Δt: the phases are drawn at t = 0 and redrawn every 5 steps.
            @test changed == 0:5:30
            # Between redraws the amplitudes stay fixed, with ⟨|f|²⟩ = f_rms².
            @test mean_square(previous...) ≈ 1
        end

        @testset "solenoidal output with compressive fraction 0" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            for f in (ou(g, band(1, 4)), RandomPhaseForcing(g, band(1, 4); f_rms = 1, τ_r = 1, seed = 5))
                state, ws = rest(g)
                advance!(f, state, ws, 0.01)
                @test divergence(g, forcing_field(ws)...) < 1e-15
            end
            # With compressive fraction 1 the forcing is curl-free: k × f̂ = 0.
            f = RandomPhaseForcing(g, band(1, 4); f_rms = 1, τ_r = 1, seed = 5, compressive_fraction = 1)
            state, ws = rest(g)
            advance!(f, state, ws, 0.01)
            fx, fy = forcing_field(ws)
            @test maximum(abs.(Array(g.kx) .* fy .- Array(g.ky) .* fx)) < 1e-15
        end

        @testset "compressive fraction only with the compressible model" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            compressive = ou(g, band(1, 2); compressive_fraction = 0.5)
            @test_throws ArgumentError Problem(g, Incompressible(ν = 0.1), compressive)
            @test_throws ArgumentError Problem(g, Incompressible(ν = 0.1),
                                               RandomPhaseForcing(g, band(1, 2); f_rms = 1, τ_r = 1, seed = 1,
                                                                  compressive_fraction = 0.5))
            @test Problem(g, Compressible(; eos = Isothermal(cs = 1.0), μ = 0.1, ζ = 0.0), compressive).forcing ===
                  compressive
            @test Problem(g, Incompressible(ν = 0.1), ou(g, band(1, 2))) isa Problem
            @test_throws ArgumentError ou(g, band(1, 2); compressive_fraction = 1.5)
        end

        @testset "add_forcing! adds f̂ on the forced modes only" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            f = ou(g, band(2, 3); compressive_fraction = 0.3)
            state = State(merge(zero_velocity(g), uniform_density(g)))
            ws = Workspace(g, state)
            advance!(f, state, ws, 0.1)
            host = (ux = rand(ComplexF64, 9, 16), uy = rand(ComplexF64, 9, 16), rho = rand(ComplexF64, 9, 16))
            dû = map(a -> copyto!(spectral_zeros(g), a), host)
            @test add_forcing!(dû, f, state, ws) === dû
            fx, fy = forcing_field(ws)
            @test Array(dû.ux) == host.ux .+ fx && Array(dû.uy) == host.uy .+ fy
            @test Array(dû.rho) == host.rho
        end

        @testset "forcing state restores the sequence" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            for build in (() -> ou(g, band(1, 3); compressive_fraction = 0.2),
                          () -> RandomPhaseForcing(g, band(1, 3); f_rms = 1, τ_r = 0.03, seed = 2))
                f = build()
                state, ws = rest(g)
                foreach(_ -> advance!(f, state, ws, 0.01), 1:7)
                saved, t, n = forcing_state(f), state.t, state.n
                reference = [next_field!(f, state, ws, 0.01) for _ in 1:7]
                # A new forcing of the same parameters takes the saved state and repeats the sequence.
                restarted = build()
                state2, ws2 = rest(g)
                state2.t, state2.n = t, n
                @test set_forcing_state!(restarted, saved) === restarted
                @test [next_field!(restarted, state2, ws2, 0.01) for _ in 1:7] == reference
            end
            @test forcing_state(KolmogorovForcing(g; F = 1, n = 1)) == (;)
            @test forcing_state(NoForcing()) == (;)
        end

        @testset "convert_float and type stability" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float32, backend)
            state, ws = rest(g)
            for f in (ou(Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend), band(1, 3)),
                      RandomPhaseForcing(Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend), band(1, 3);
                                         f_rms = 1, τ_r = 1, seed = 1),
                      KolmogorovForcing(Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend); F = 1, n = 2))
                f32 = convert_float(Float32, f)
                # Every float field and array of the converted forcing is Float32.
                floats = [getfield(f32, k) for k in fieldnames(typeof(f32)) if getfield(f32, k) isa AbstractFloat]
                arrays = [getfield(f32, k) for k in fieldnames(typeof(f32)) if getfield(f32, k) isa AbstractArray]
                @test all(x -> x isa Float32, floats)
                @test all(a -> real(eltype(a)) == Float32, arrays)
                @inferred update_forcing!(f32, state, ws, 0.01f0)
                @test eltype(ws.forcing.ux) == ComplexF32
            end
        end

        @testset "OU normalisation" begin
            if Main.is_cpu(backend)
                # 99.9% interval over M independent realisation means (ADR-0011, ADR-0017).
                M = 4000
                means = ou_realisation_means(; M)
                mean = sum(means) / M
                std = sqrt(sum(abs2, means .- mean) / (M - 1))
                @test abs(mean - 1) <= 3.3 * std / sqrt(M)
            else
                # The device reproduces the CPU forcing sequence for the same seed.
                sequences = map((CPU(), backend)) do b
                    g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, backend = b)
                    f = OUForcing(g, band(4, 6); ε = 1, τ_f = 1, seed = 7, compressive_fraction = 0.4)
                    state = State(merge(zero_velocity(g), uniform_density(g)))
                    ws = Workspace(g, state)
                    [next_field!(f, state, ws, 0.01) for _ in 1:100]
                end
                @test all(((cpu, device),) -> all(map((a, b) -> isapprox(a, b; rtol = 1e-12), cpu, device)),
                          zip(sequences...))
            end
        end
    end
end

end
