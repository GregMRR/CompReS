# Tests of the State and the Workspace (WP-4): field names, time and step, named buffers per state field.
module StateTests

using Test
using KernelAbstractions
using CompReS

spectral_zeros(g, T) = KernelAbstractions.zeros(g.backend, Complex{T}, g.Nx ÷ 2 + 1, g.Ny)

function run_tests(backend)
    @testset "state" begin
        @testset "State holds the fields, the time and the step" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, backend)
            fields = (ux = spectral_zeros(g, Float64), uy = spectral_zeros(g, Float64))
            state = State(fields)
            # The State keeps the given arrays, not copies, and starts at t = 0, n = 0 by default.
            @test state.fields.ux === fields.ux && state.fields.uy === fields.uy
            @test state.t === 0.0 && state.n === 0
            # The time loop advances t and n in place.
            state.t += 0.5
            state.n += 1
            @test state.t === 0.5 && state.n === 1
            # The keywords set the start time and step; t takes the float type of the fields.
            restarted = State(fields; t = 2, n = 7)
            @test restarted.t === 2.0 && restarted.n === 7
        end

        @testset "State accepts the two models and rejects other field sets" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, backend)
            f̂() = spectral_zeros(g, Float64)
            # The compressible model adds rho after the velocity.
            @test keys(State((ux = f̂(), uy = f̂(), rho = f̂())).fields) == (:ux, :uy, :rho)
            # Other names or another order do not describe a model.
            @test_throws ArgumentError State((ux = f̂(),))
            @test_throws ArgumentError State((uy = f̂(), ux = f̂()))
            @test_throws ArgumentError State((ux = f̂(), uy = f̂(), T = f̂()))
            # All fields must share one size and one element type.
            small = KernelAbstractions.zeros(backend, ComplexF64, 3, 8)
            @test_throws ArgumentError State((ux = f̂(), uy = small))
            @test_throws ArgumentError State((ux = f̂(), uy = spectral_zeros(g, Float32)))
        end

        @testset "Float32 State stays Float32" begin
            # The time takes the float type of the fields, so a Float32 run never stores a Float64 time.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, T = Float32, backend)
            state = @inferred State((ux = spectral_zeros(g, Float32), uy = spectral_zeros(g, Float32)); t = 0.1)
            @test state.t === 0.1f0
        end

        @testset "Workspace has one physical buffer per state field" begin
            g = Grid(; Lx = 4π, Ly = 2π, Nx = 16, Ny = 8, T = Float32, backend)
            f̂() = spectral_zeros(g, Float32)
            for state in (State((ux = f̂(), uy = f̂())), State((ux = f̂(), uy = f̂(), rho = f̂())))
                ws = @inferred Workspace(g, state)
                # The buffers carry the names of the state fields, are real Nx × Ny arrays of type T on the
                # backend of the grid, and share no storage.
                @test keys(ws.physical) == keys(state.fields)
                @test all(b -> size(b) == (16, 8) && eltype(b) == Float32, ws.physical)
                @test all(b -> KernelAbstractions.get_backend(b) == backend, ws.physical)
                @test allunique(map(objectid, values(ws.physical)))
                # The Transforms belong to the same grid.
                @test ws.transforms.grid === g
            end
        end
    end
end

end
