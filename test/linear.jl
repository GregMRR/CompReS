# Tests of the linear operators and the integrating factor (WP-5): single Fourier modes on N = 16.
module LinearTests

using Test
using Random
using KernelAbstractions
using CompReS

# Array position of the mode with integer indices (i, j) in the real-FFT layout.
pos(g, i, j) = (i + 1, j < 0 ? g.Ny + j + 1 : j + 1)

# Device copy of a host array.
device(backend, a) = copyto!(KernelAbstractions.allocate(backend, eltype(a), size(a)), a)

host_zeros(g) = zeros(Complex{eltype(g.k2)}, g.Nx ÷ 2 + 1, g.Ny)

# Device fields with the given names, each from a host spectral array.
fields(g, names, arrays) = NamedTuple{names}(Tuple(device(g.backend, a) for a in arrays))

# Device fields with the given names and zero values.
zero_fields(g, names) = fields(g, names, Tuple(host_zeros(g) for _ in names))

# Host copies of every field.
host(f) = map(Array, f)

# Fields that hold the single mode (i, j): velocity amplitude (ax, ay) and density amplitude aρ.
function single_mode(g, names, i, j, ax, ay, aρ = 0)
    arrays = Tuple(host_zeros(g) for _ in names)
    arrays[1][pos(g, i, j)...] = ax
    arrays[2][pos(g, i, j)...] = ay
    length(names) == 3 && (arrays[3][pos(g, i, j)...] = aρ)
    return fields(g, names, arrays)
end

# Physical wavenumber vector of the mode (i, j).
wavevector(g, i, j) = (2π * i / g.Lx, 2π * j / g.Ly)

# Random fields, nonzero on every mode.
random_fields(g, names, seed) =
    fields(g, names, Tuple(randn(MersenneTwister(seed + n), Complex{eltype(g.k2)}, g.Nx ÷ 2 + 1, g.Ny)
                           for n in eachindex(names)))

# Every float stored in `x`, also inside nested structs, has the type T.
holds_only(T, x::AbstractFloat) = x isa T
holds_only(T, x::Union{Integer,Bool}) = true
holds_only(T, x) = all(holds_only(T, getfield(x, f)) for f in fieldnames(typeof(x)))

const INC = (:ux, :uy)
const COMP = (:ux, :uy, :rho)

# Apply the propagator for the time τ to `u` and return host copies; also return L·u from the tendency.
function propagated_and_tendency(g, rates, u, τ)
    names = keys(u)
    out = zero_fields(g, names)
    propagate!(out, u, rates, τ)
    du = zero_fields(g, names)
    linear_tendency!(du, u, rates)
    return host(out), host(du)
end

function run_tests(backend)
    @testset "linear operators" begin
        @testset "model constructors" begin
            # The keyword defaults give no large-scale dissipation and the full equations.
            m = Incompressible(ν = 1e-3)
            @test m isa Incompressible{Float64}
            @test m.large_scale.α == 0 && m.nonlinear
            @test Incompressible(ν = 1e-3, nonlinear = false).nonlinear == false
            c = Compressible(eos = Isothermal(cs = 10.0), μ = 1e-3, ζ = 0.5e-3)
            @test c.κ == 0 && c.large_scale.α == 0 && c.nonlinear
            # ADR-0004 excludes hyperviscosity, so a negative exponent m is rejected.
            @test_throws ArgumentError LargeScaleDissipation(α = 1, m = -1)
            # A negative coefficient amplifies modes instead of dissipating them.
            @test_throws ArgumentError Incompressible(ν = -1e-3)
            @test_throws ArgumentError Compressible(eos = Isothermal(cs = 1.0), μ = 0, ζ = 0, κ = -1)
            @test_throws ArgumentError LargeScaleDissipation(α = -1, m = 0)
        end

        @testset "convert_float holds only the target type" begin
            # ADR-0006: a Float32 run must not promote, so every float inside the model must be Float32.
            ls = LargeScaleDissipation(α = 0.1, m = 1)
            for model in (Incompressible(ν = 1e-3, large_scale = ls),
                          Compressible(eos = Isothermal(cs = 10.0), μ = 1e-3, ζ = 2e-3, κ = 1e-4, large_scale = ls),
                          Compressible(eos = Polytropic(K = 1.0, γ = 1.4), μ = 1e-3, ζ = 0, large_scale = ls))
                m32 = convert_float(Float32, model)
                @test holds_only(Float32, m32)
                @test !holds_only(Float32, model)
                @test isbits(m32)
                # The conversion keeps the exponent m and the switch nonlinear unchanged.
                @test m32.large_scale.m == 1 && m32.nonlinear == model.nonlinear
            end
            @test convert_float(Float32, ls) === LargeScaleDissipation(α = 0.1f0, m = 1)
        end

        g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
        τ = 0.3

        @testset "incompressible single mode decays at νk²" begin
            # Mode (2, 1) with the solenoidal amplitude (ky, -kx) = (1, -2): k² = 5, rate νk² = 0.5.
            ν = 0.1
            rates = linear_rates(g, Incompressible(ν = ν))
            u = single_mode(g, INC, 2, 1, 1.0, -2.0)
            out, du = propagated_and_tendency(g, rates, u, τ)
            p = pos(g, 2, 1)
            @test out.ux[p...] ≈ exp(-0.5τ) && out.uy[p...] ≈ -2exp(-0.5τ)
            @test du.ux[p...] ≈ -0.5 && du.uy[p...] ≈ 1.0
            # The operator acts mode by mode, so every other mode stays zero.
            @test count(!iszero, out.ux) == 1 && count(!iszero, du.uy) == 1
        end

        @testset "large-scale dissipation decays at α k^(-2m)" begin
            # Rectangle Lx = 4π: mode (2, 1) has k = (1, 1), k² = 2, so α k^(-2m) = 0.4/2^m.
            gr = Grid(; Lx = 4π, Ly = 2π, Nx = 32, Ny = 16, backend)
            for m in (0, 1, 2)
                rates = linear_rates(gr, Incompressible(ν = 0, large_scale = LargeScaleDissipation(α = 0.4, m = m)))
                out, du = propagated_and_tendency(gr, rates, single_mode(gr, INC, 2, 1, 1.0, -1.0), τ)
                σ = 0.4 / 2^m
                @test out.ux[pos(gr, 2, 1)...] ≈ exp(-σ * τ)
                @test du.uy[pos(gr, 2, 1)...] ≈ σ
            end
            # Small-scale and large-scale dissipation add: k² = 2 gives νk² + α/k² = 0.2 + 0.2.
            rates = linear_rates(gr, Incompressible(ν = 0.1, large_scale = LargeScaleDissipation(α = 0.4, m = 1)))
            out, _ = propagated_and_tendency(gr, rates, single_mode(gr, INC, 2, 1, 1.0, -1.0), τ)
            @test out.ux[pos(gr, 2, 1)...] ≈ exp(-0.4τ)
        end

        ls = LargeScaleDissipation(α = 0.2, m = 1)
        model = Compressible(eos = Isothermal(cs = 10.0), μ = 0.1, ζ = 0.05, κ = 0.03, large_scale = ls)
        rates = linear_rates(g, model)

        @testset "compressible solenoidal mode decays at μk², compressive mode at (μ+ζ)k²" begin
            # Mode (2, 1): k = (2, 1), k² = 5, α/k² = 0.04. Solenoidal rate μk² + α/k² = 0.54.
            p = pos(g, 2, 1)
            out, du = propagated_and_tendency(g, rates, single_mode(g, COMP, 2, 1, 1.0, -2.0), τ)
            @test out.ux[p...] ≈ exp(-0.54τ) && out.uy[p...] ≈ -2exp(-0.54τ)
            @test du.ux[p...] ≈ -0.54
            # A compressive mode, parallel to k = (2, 1), has the rate (μ + ζ)k² + α/k² = 0.79.
            out, du = propagated_and_tendency(g, rates, single_mode(g, COMP, 2, 1, 2.0, 1.0), τ)
            @test out.ux[p...] ≈ 2exp(-0.79τ) && out.uy[p...] ≈ exp(-0.79τ)
            @test du.uy[p...] ≈ -0.79
            # The operator is diagonal in the Helmholtz basis: the sum of both modes decays part by part.
            out, du = propagated_and_tendency(g, rates, single_mode(g, COMP, 2, 1, 3.0, -1.0), τ)
            @test out.ux[p...] ≈ exp(-0.54τ) + 2exp(-0.79τ)
            @test out.uy[p...] ≈ -2exp(-0.54τ) + exp(-0.79τ)
            @test du.ux[p...] ≈ -0.54 - 2 * 0.79
        end

        @testset "density decays at κk²" begin
            # Mode (0, 3): k² = 9, κk² = 0.27. Large-scale dissipation acts on u only, not on ρ.
            out, du = propagated_and_tendency(g, rates, single_mode(g, COMP, 0, 3, 0, 0, 1.0), τ)
            @test out.rho[pos(g, 0, 3)...] ≈ exp(-0.27τ)
            @test du.rho[pos(g, 0, 3)...] ≈ -0.27
        end

        @testset "k = 0 untouched" begin
            # The mean flow and the mean density have k = 0: small-scale terms vanish and ADR-0004 excludes
            # the large-scale term, also for linear friction m = 0, so the propagator is exactly 1.
            for m in (0, 1)
                lsm = LargeScaleDissipation(α = 0.5, m = m)
                for (names, mdl) in ((INC, Incompressible(ν = 0.1, large_scale = lsm)),
                                     (COMP, Compressible(eos = Isothermal(cs = 1.0), μ = 0.1, ζ = 0.1, κ = 0.1,
                                                         large_scale = lsm)))
                    u = single_mode(g, names, 0, 0, 0.7, -0.3, 1.0)
                    out, du = propagated_and_tendency(g, linear_rates(g, mdl), u, τ)
                    @test out == host(u)
                    @test all(iszero, du.ux) && all(iszero, du.uy)
                end
            end
        end

        @testset "propagator equals exp of the rates" begin
            # On random fields, exp(Lτ₁)exp(Lτ₂) = exp(L(τ₁ + τ₂)) and exp(L·0) = I hold for an exponential.
            for (names, mdl) in ((INC, Incompressible(ν = 0.1, large_scale = ls)), (COMP, model))
                r = linear_rates(g, mdl)
                u = random_fields(g, names, 1)
                a, b = zero_fields(g, names), zero_fields(g, names)
                propagate!(a, u, r, 0.1)
                propagate!(a, a, r, 0.2)
                propagate!(b, u, r, 0.3)
                @test all(map((x, y) -> isapprox(x, y; rtol = 1e-14), host(a), host(b)))
                propagate!(a, u, r, 0.0)
                @test host(a) == host(u)
                # The derivative of exp(Lτ)u at τ = 0 is Lu, the explicit tendency.
                δ = 1e-6
                propagate!(a, u, r, δ)
                propagate!(b, u, r, -δ)
                du = zero_fields(g, names)
                linear_tendency!(du, u, r)
                @test all(map((x, y, z) -> isapprox((x - y) / 2δ, z; rtol = 1e-8), host(a), host(b), host(du)))
            end
            # On the incompressible model the propagator is exp(-σ τ) with the stored rate array σ.
            r = linear_rates(g, Incompressible(ν = 0.1, large_scale = ls))
            u = random_fields(g, INC, 2)
            out = zero_fields(g, INC)
            propagate!(out, u, r, τ)
            @test Array(out.ux) ≈ exp.(-Array(r.velocity) .* τ) .* Array(u.ux) rtol = 1e-14
        end

        @testset "linear_tendency! adds to the tendency" begin
            # Explicit schemes sum the linear and nonlinear tendencies, so the function adds and does not overwrite.
            u = single_mode(g, COMP, 0, 3, 0, 0, 1.0)
            du = single_mode(g, COMP, 0, 3, 0, 0, 5.0)
            linear_tendency!(du, u, rates)
            @test Array(du.rho)[pos(g, 0, 3)...] ≈ 5 - 0.27
        end

        @testset "Float32 rates and propagator" begin
            # A Float32 grid gives Float32 rates, and a Float64 time does not promote the result.
            g32 = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float32, backend)
            r32 = linear_rates(g32, convert_float(Float32, model))
            @test eltype(r32.solenoidal) == eltype(r32.compressive) == eltype(r32.density) == Float32
            out, du = propagated_and_tendency(g32, r32, single_mode(g32, COMP, 2, 1, 1.0f0, -2.0f0), τ)
            @test eltype(out.ux) == ComplexF32
            @test out.ux[pos(g32, 2, 1)...] ≈ exp(-0.54f0 * Float32(τ)) rtol = 10 * eps(Float32)
            @test du.ux[pos(g32, 2, 1)...] ≈ -0.54f0 rtol = 10 * eps(Float32)
        end
    end
end

end
