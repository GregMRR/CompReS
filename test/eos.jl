# Tests of the equation of state (WP-3): analytic values, derivative identities, round trip, device kernel.
module EosTests

using Test
using KernelAbstractions
using CompReS

@kernel function eos_kernel!(out, eos, ρ)
    i = @index(Global)
    out[i, 1] = pressure(eos, ρ[i])
    out[i, 2] = enthalpy(eos, ρ[i])
    out[i, 3] = internal_energy(eos, ρ[i])
    out[i, 4] = sound_speed²(eos, ρ[i])
    out[i, 5] = density_from_enthalpy(eos, enthalpy(eos, ρ[i]))
end

function run_tests(backend)
    @testset "Isothermal analytic values" begin
        # cs = 3, ρ = 2: p = cs²ρ = 18, h = e = cs² ln ρ = 9 ln 2, c_s² = cs² = 9.
        eos = Isothermal(cs = 3.0)
        @test pressure(eos, 2.0) == 18.0
        @test enthalpy(eos, 2.0) ≈ 9 * log(2.0)
        @test internal_energy(eos, 2.0) ≈ 9 * log(2.0)
        @test sound_speed²(eos, 2.0) == 9.0
        # The specification fixes ρ₀ = 1, so h(1) = e(1) = 0.
        @test enthalpy(eos, 1.0) == 0.0
        @test internal_energy(eos, 1.0) == 0.0
        # Inverse of h = 9 ln 2 at cs² = 9 is ρ = exp(ln 2) = 2.
        @test density_from_enthalpy(eos, 9 * log(2.0)) ≈ 2.0
    end

    @testset "Polytropic analytic values" begin
        # K = 2, γ = 2, ρ = 3: p = 18, h = Kγ/(γ-1) ρ = 12, e = K/(γ-1) ρ = 6, c_s² = Kγρ = 12.
        eos = Polytropic(K = 2.0, γ = 2.0)
        @test pressure(eos, 3.0) == 18.0
        @test enthalpy(eos, 3.0) == 12.0
        @test internal_energy(eos, 3.0) == 6.0
        @test sound_speed²(eos, 3.0) == 12.0
        @test density_from_enthalpy(eos, 12.0) ≈ 3.0
        # γ = 1 makes h and e singular; the isothermal law covers that case.
        @test_throws ArgumentError Polytropic(K = 1.0, γ = 1)
        # K = 4, γ = 3/2, ρ = 4 gives ρ^(γ-1) = 2: p = 32, h = 24, e = 16, c_s² = 12.
        eos = Polytropic(K = 4.0, γ = 1.5)
        @test pressure(eos, 4.0) ≈ 32.0
        @test enthalpy(eos, 4.0) ≈ 24.0
        @test internal_energy(eos, 4.0) ≈ 16.0
        @test sound_speed²(eos, 4.0) ≈ 12.0
        @test density_from_enthalpy(eos, 24.0) ≈ 4.0
    end

    @testset "derivative identities" begin
        # The specification gives h' = p'/ρ, e' = p/ρ² and c_s² = p'; a central difference checks each.
        δ = 1e-5
        d(f, eos, ρ) = (f(eos, ρ + δ) - f(eos, ρ - δ)) / 2δ
        for eos in (Isothermal(cs = 3.0), Polytropic(K = 2.0, γ = 1.4)), ρ in (0.5, 1.0, 2.5)
            @test d(enthalpy, eos, ρ) ≈ d(pressure, eos, ρ) / ρ rtol = 1e-7
            @test d(internal_energy, eos, ρ) ≈ pressure(eos, ρ) / ρ^2 rtol = 1e-7
            @test sound_speed²(eos, ρ) ≈ d(pressure, eos, ρ) rtol = 1e-7
        end
        # The specification gives h = e + p/ρ; the polytropic constants satisfy it exactly.
        eos = Polytropic(K = 2.0, γ = 1.4)
        @test enthalpy(eos, 2.5) ≈ internal_energy(eos, 2.5) + pressure(eos, 2.5) / 2.5
        # The isothermal constants differ by c_s² = p/ρ, which the specification allows.
        eos = Isothermal(cs = 3.0)
        @test enthalpy(eos, 2.5) - internal_energy(eos, 2.5) == 0
    end

    @testset "density_from_enthalpy inverts enthalpy" begin
        # h is strictly increasing in ρ, so the inverse recovers ρ to round-off.
        for eos in (Isothermal(cs = 10.0), Polytropic(K = 0.7, γ = 1.4)), ρ in (0.1, 0.9, 1.0, 1.1, 7.0)
            @test density_from_enthalpy(eos, enthalpy(eos, ρ)) ≈ ρ rtol = 1e-13
        end
    end

    @testset "parameters and float type" begin
        # Integer parameters become Float64; the keyword constructor takes the float type from its arguments.
        @test Isothermal(cs = 10) isa Isothermal{Float64}
        @test Polytropic(K = 1, γ = 2) isa Polytropic{Float64}
        @test Isothermal(cs = 2.0f0) isa Isothermal{Float32}
        @test Polytropic(K = 1.0f0, γ = 2.0f0) isa Polytropic{Float32}
        # convert_float changes the parameter type and keeps the values that Float32 can hold exactly.
        @test convert_float(Float32, Isothermal(cs = 2.0)) === Isothermal(cs = 2.0f0)
        @test convert_float(Float32, Polytropic(K = 0.5, γ = 1.5)) === Polytropic(K = 0.5f0, γ = 1.5f0)
        @test convert_float(Float64, Isothermal(cs = 2.0f0)) === Isothermal(cs = 2.0)
        # Float32 parameters and a Float32 density give a Float32 result for every method (no promotion).
        for eos in (Isothermal(cs = 3.0f0), Polytropic(K = 2.0f0, γ = 1.4f0)), f in
            (pressure, enthalpy, internal_energy, sound_speed², density_from_enthalpy)
            @test @inferred(f(eos, 2.0f0)) isa Float32
        end
        # Both structs are isbits, so a kernel can take them as arguments.
        @test isbitstype(Isothermal{Float32})
        @test isbitstype(Polytropic{Float64})
    end

    @testset "kernel on the device" begin
        # The kernel output equals the host value of each method; column 5 is the round trip, which returns ρ.
        for T in (Float64, Float32), eos in (Isothermal(cs = T(3)), Polytropic(K = T(2), γ = T(1.4)))
            ρ_host = T[0.5, 1, 2, 3]
            ρ = KernelAbstractions.allocate(backend, T, length(ρ_host))
            copyto!(ρ, ρ_host)
            out = KernelAbstractions.zeros(backend, T, length(ρ_host), 5)
            eos_kernel!(backend)(out, eos, ρ; ndrange = length(ρ_host))
            KernelAbstractions.synchronize(backend)
            expected = [f(eos, r) for r in ρ_host, f in (pressure, enthalpy, internal_energy, sound_speed²)]
            expected = hcat(expected, ρ_host)
            @test eltype(out) == T
            @test Array(out) ≈ expected rtol = 10 * eps(T)
        end
    end
end

end
