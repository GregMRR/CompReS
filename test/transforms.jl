# Tests of the transforms and spectral operators (WP-2): single Fourier modes on N = 8 or 16, plus identities on random fields.
module TransformsTests

using Test
using Random
using KernelAbstractions
using CompReS

# Array position of the mode with integer indices (i, j) in the real-FFT layout.
pos(g, i, j) = (i + 1, j < 0 ? g.Ny + j + 1 : j + 1)

# Device copy of a host array.
device(backend, a) = copyto!(KernelAbstractions.allocate(backend, eltype(a), size(a)), a)

# Physical grid points x and y of the grid, as host arrays of size (Nx, 1) and (1, Ny).
coordinates(g) = ((0:(g.Nx - 1)) .* (g.Lx / g.Nx), ((0:(g.Ny - 1)) .* (g.Ly / g.Ny))')

spectral_zeros(g, T) = KernelAbstractions.zeros(g.backend, Complex{T}, g.Nx ÷ 2 + 1, g.Ny)

# Spectral field of a physical host array.
function spectral(tr, f)
    f̂ = spectral_zeros(tr.grid, eltype(f))
    return to_spectral!(f̂, device(tr.grid.backend, f), tr)
end

# Physical host array of a spectral field.
function physical(tr, f̂)
    g = tr.grid
    f = KernelAbstractions.zeros(g.backend, real(eltype(f̂)), g.Nx, g.Ny)
    return Array(to_physical!(f, f̂, tr))
end

# Half-plane weights of the real FFT: 1 on the kx = 0 and Nyquist columns, 2 elsewhere.
half_plane_weights(g) = [c == 1 || c == g.Nx ÷ 2 + 1 ? 1 : 2 for c in 1:(g.Nx ÷ 2 + 1), _ in 1:g.Ny]

function run_tests(backend)
    @testset "transforms" begin
        @testset "normalisation on single modes" begin
            # cos(3x) = (e^{3ix} + e^{-3ix})/2 and 2 sin(2y) = -i e^{2iy} + i e^{-2iy}; the stored half-plane
            # holds (3, 0) = 1/2, (0, 2) = -i and (0, -2) = i, and every other mode is zero.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, backend)
            tr = Transforms(g)
            x, y = coordinates(g)
            f̂ = Array(spectral(tr, cos.(3 .* x) .+ 2 .* sin.(2 .* y)))
            expected = zeros(ComplexF64, 5, 8)
            expected[pos(g, 3, 0)...] = 1 / 2
            expected[pos(g, 0, 2)...] = -im
            expected[pos(g, 0, -2)...] = im
            @test f̂ ≈ expected atol = 10 * eps()
        end

        @testset "round trip" begin
            # The inverse of the normalised forward transform is the plain Fourier sum, so the round trip is the identity.
            for g in (Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend), Grid(; Lx = 4π, Ly = 2π, Nx = 16, Ny = 8, backend))
                tr = Transforms(g)
                f = randn(MersenneTwister(1), g.Nx, g.Ny)
                f̂ = spectral(tr, f)
                @test physical(tr, f̂) ≈ f rtol = 10 * eps()
                # to_physical! must leave its input unchanged, although the complex-to-real FFT overwrites its input.
                saved = Array(f̂)
                physical(tr, f̂)
                @test Array(f̂) == saved
            end
        end

        @testset "derivatives of sin(3x)cos(2y)" begin
            # ∂x = 3 cos(3x) cos(2y) and ∂y = -2 sin(3x) sin(2y); both are band-limited, so the spectral result is exact.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            tr = Transforms(g)
            x, y = coordinates(g)
            f̂ = spectral(tr, sin.(3 .* x) .* cos.(2 .* y))
            df̂ = spectral_zeros(g, Float64)
            @test physical(tr, ∂x!(df̂, f̂, tr)) ≈ 3 .* cos.(3 .* x) .* cos.(2 .* y) rtol = 10 * eps()
            @test physical(tr, ∂y!(df̂, f̂, tr)) ≈ -2 .* sin.(3 .* x) .* sin.(2 .* y) rtol = 10 * eps()
            # The output may be the input: ∂x then ∂x in place gives -9 f.
            ∂x!(f̂, ∂x!(f̂, f̂, tr), tr)
            @test physical(tr, f̂) ≈ -9 .* sin.(3 .* x) .* cos.(2 .* y) rtol = 10 * eps()
        end

        @testset "derivative multipliers are zero on the Nyquist column and row" begin
            # The derivative wavenumbers equal the grid wavenumbers except at kx = Nx/2 and ky = Ny/2, where they are zero.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            tr = Transforms(g)
            @test Array(tr.kx_deriv) == [0:7; 0]
            @test Array(tr.ky_deriv) == [0:7; 0; -7:-1]
            # cos(8x) is the Nyquist column mode (8, 0) with coefficient 1; its spectral x derivative must be zero
            # before any inverse FFT, which would hide a nonzero value.
            x, y = coordinates(g)
            f̂ = spectral(tr, cos.(8 .* x) .+ 0 .* y)
            @test Array(f̂)[pos(g, 8, 0)...] ≈ 1
            @test all(iszero, Array(∂x!(spectral_zeros(g, Float64), f̂, tr)))
            # cos(8y) is the Nyquist row mode (0, -8) with coefficient 1; the same holds for ∂y!.
            f̂ = spectral(tr, 0 .* x .+ cos.(8 .* y))
            @test Array(f̂)[pos(g, 0, -8)...] ≈ 1
            @test all(iszero, Array(∂y!(spectral_zeros(g, Float64), f̂, tr)))
        end

        @testset "enforce_hermitian!" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            tr = Transforms(g)
            rng = MersenneTwister(2)
            source = randn(rng, ComplexF64, 9, 16)
            f̂ = enforce_hermitian!(device(backend, source), g)
            h = Array(f̂)
            # The kx = 0 column copies the ky > 0 half: (0, -j) = conj((0, j)) for j = 1…7, and the source half is unchanged.
            @test all(h[1, 16 + 2 - r] == conj(source[1, r]) && h[1, r] == source[1, r] for r in 2:8)
            # The k = 0 mode keeps its real part only.
            @test h[1, 1] == real(source[1, 1])
            # The Nyquist column kx = 8 and the Nyquist row ky = -8 are zero; the remaining modes are unchanged.
            @test all(iszero, h[9, :]) && all(iszero, h[:, 9])
            @test h[2:8, [1:8; 10:16]] == source[2:8, [1:8; 10:16]]
            # A spectral field of a real field with zero Nyquist modes survives the inverse and forward FFT unchanged.
            @test Array(spectral(tr, physical(tr, f̂))) ≈ h rtol = 10 * eps()
        end

        @testset "skew-adjoint derivatives on random fields" begin
            # Summation by parts holds for the discrete derivative only when the Nyquist multipliers are zero;
            # the fields are not truncated, so the Nyquist modes are nonzero.
            g = Grid(; Lx = 4π, Ly = 2π, Nx = 16, Ny = 8, backend)
            tr = Transforms(g)
            rng = MersenneTwister(3)
            f, h = randn(rng, g.Nx, g.Ny), randn(rng, g.Nx, g.Ny)
            for d! in (∂x!, ∂y!)
                df = physical(tr, d!(spectral_zeros(g, Float64), spectral(tr, f), tr))
                dh = physical(tr, d!(spectral_zeros(g, Float64), spectral(tr, h), tr))
                @test abs(sum(f .* dh) + sum(h .* df)) <= 100 * eps() * sum(abs.(f .* dh))
            end
        end

        @testset "truncate!" begin
            # Modes outside the disc become zero; modes inside keep their value.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            source = randn(MersenneTwister(4), ComplexF64, 9, 16)
            h = Array(truncate!(device(backend, source), g))
            keep = Array(g.truncation)
            @test h[keep] == source[keep]
            @test all(iszero, h[.!keep])
        end

        @testset "Helmholtz split" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            tr = Transforms(g)
            x, y = coordinates(g)
            full(a) = a .+ zeros(g.Nx, g.Ny)
            field_pair(a, b) = (spectral(tr, full(a)), spectral(tr, full(b)))
            pair() = (spectral_zeros(g, Float64), spectral_zeros(g, Float64))
            # u = (cos x, 0) is parallel to k = (±1, 0), so it is purely compressive.
            ûs, ûc = helmholtz_split!(pair(), pair(), field_pair(cos.(x), 0), tr)
            @test maximum(abs, Array(ûs[1])) <= eps() && maximum(abs, Array(ûs[2])) <= eps()
            @test physical(tr, ûc[1]) ≈ full(cos.(x)) atol = 10 * eps()
            # u = (cos y, 0) is normal to k = (0, ±1), so it is purely solenoidal.
            ûs, ûc = helmholtz_split!(pair(), pair(), field_pair(cos.(y), 0), tr)
            @test maximum(abs, Array(ûc[1])) <= eps() && maximum(abs, Array(ûc[2])) <= eps()
            @test physical(tr, ûs[1]) ≈ full(cos.(y)) atol = 10 * eps()

            # Random field, not truncated; the k = 0 mode (1, 2) is a uniform flow, divergence-free and curl-free.
            rng = MersenneTwister(5)
            û = field_pair(randn(rng, g.Nx, g.Ny), randn(rng, g.Nx, g.Ny))
            for (a, v) in zip(û, (1, 2))
                a[1:1, 1:1] .= v
            end
            ûs, ûc = helmholtz_split!(pair(), pair(), û, tr)
            # The parts add up to the field.
            @test Array(ûs[1]) .+ Array(ûc[1]) ≈ Array(û[1]) rtol = 10 * eps()
            @test Array(ûs[2]) .+ Array(ûc[2]) ≈ Array(û[2]) rtol = 10 * eps()
            # The k = 0 mode is finite and goes to the solenoidal part.
            @test Array(ûs[1])[1, 1] == 1 && Array(ûs[2])[1, 1] == 2
            @test Array(ûc[1])[1, 1] == 0 && Array(ûc[2])[1, 1] == 0
            @test all(isfinite, Array(ûs[1])) && all(isfinite, Array(ûc[2]))
            # k·ûs = 0: the divergence of the solenoidal part, built with ∂x! and ∂y!, is zero to round-off.
            div = Array(∂x!(spectral_zeros(g, Float64), ûs[1], tr)) .+ Array(∂y!(spectral_zeros(g, Float64), ûs[2], tr))
            @test maximum(abs, div) <= 100 * eps() * maximum(abs, Array(û[1]))
            # k × ûc = 0: the curl of the compressive part is zero to round-off.
            curl = Array(∂x!(spectral_zeros(g, Float64), ûc[2], tr)) .- Array(∂y!(spectral_zeros(g, Float64), ûc[1], tr))
            @test maximum(abs, curl) <= 100 * eps() * maximum(abs, Array(û[1]))
            # project_solenoidal! gives the solenoidal part of the split, and a second projection changes nothing.
            project_solenoidal!(û, tr)
            @test Array(û[1]) ≈ Array(ûs[1]) rtol = 10 * eps()
            @test Array(û[2]) ≈ Array(ûs[2]) rtol = 10 * eps()
            once = Array.(û)
            project_solenoidal!(û, tr)
            @test Array(û[1]) ≈ once[1] rtol = 10 * eps()
            @test Array(û[2]) ≈ once[2] rtol = 10 * eps()
        end

        @testset "Parseval with half-plane weights" begin
            # ⟨fg⟩ = Σ over the full plane of conj(f̂) ĝ; the stored half-plane counts each mode with kx ≠ 0, Nx/2 twice.
            g = Grid(; Lx = 4π, Ly = 2π, Nx = 16, Ny = 8, backend)
            tr = Transforms(g)
            rng = MersenneTwister(6)
            f, h = randn(rng, g.Nx, g.Ny), randn(rng, g.Nx, g.Ny)
            f̂, ĥ = Array(spectral(tr, f)), Array(spectral(tr, h))
            @test sum(f .* h) / length(f) ≈ sum(half_plane_weights(g) .* real.(conj.(f̂) .* ĥ)) rtol = 100 * eps()
        end

        @testset "Float32 stays Float32" begin
            # Every operator keeps the element type, and the round trip holds to Float32 round-off.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float32, backend)
            tr = Transforms(g)
            @test eltype(tr.kx_deriv) == Float32
            f = randn(MersenneTwister(7), Float32, g.Nx, g.Ny)
            f̂ = spectral(tr, f)
            @test eltype(f̂) == ComplexF32
            @test physical(tr, f̂) ≈ f rtol = 10 * eps(Float32)
            û = (copy(f̂), copy(f̂))
            @test eltype(∂x!(spectral_zeros(g, Float32), f̂, tr)) == ComplexF32
            @test eltype(project_solenoidal!(û, tr)[1]) == ComplexF32
            @test eltype(enforce_hermitian!(f̂, g)) == ComplexF32
        end
    end
end

end
