# Tests of the diagnostics registry and the incompressible diagnostics (WP-8): registry selection, shell spectra
# of single modes, Parseval sums, dissipation rates, fluxes and the start-of-step injection rates.
module DiagnosticsTests

using Test
using KernelAbstractions
using CompReS
using CompReS: update_forcing!, kinetic_energy, kinetic_energy_spectrum_inc, kinetic_energy_flux

# Spectral fields of the physical velocity (fx(x, y), fy(x, y)) on the grid points of g.
function spectral_velocity(g, ws, fx, fy)
    T = eltype(g.k2)
    x = [T(g.Lx) * (i - 1) / g.Nx for i in 1:g.Nx]
    y = [T(g.Ly) * (j - 1) / g.Ny for j in 1:g.Ny]
    spectral(f) = to_spectral!(KernelAbstractions.allocate(g.backend, Complex{T}, g.Nx ÷ 2 + 1, g.Ny),
                               copyto!(similar(ws.physical.ux), f.(x, y')), ws.transforms)
    return (ux = spectral(fx), uy = spectral(fy))
end

# Physical host field of the spectral field f̂.
physical(f̂, ws) = Array(to_physical!(similar(ws.physical.ux), f̂, ws.transforms))

# Spectral curl i(kx âᵧ − ky âₓ) of the spectral vector field a, on the device.
function spectral_curl(g, a)
    kx, ky = Array(g.kx), Array(g.ky)
    host = im .* (kx .* Array(a.uy) .- ky .* Array(a.ux))
    return copyto!(similar(a.ux), host)
end

# Grid means ⟨a·b⟩ of two physical vector fields and ⟨a b⟩ of two physical scalar fields.
mean_dot(a, b) = sum(a[1] .* b[1] .+ a[2] .* b[2]) / length(a[1])
mean_product(a, b) = sum(a .* b) / length(a)

# Value of the GLOBAL or SPECTRAL diagnostic `name` for the state and problem.
evaluate(class, name, state, problem, ws) = only(select_diagnostics(class, (name,))).second(state, problem, ws)

function run_tests(backend)
    @testset "diagnostics" begin
        @testset "registry: defaults, selection by name and custom diagnostics" begin
            # ADR-0009 defaults of the incompressible model; the _comp names are absent.
            model = Incompressible(ν = 0.1)
            @test default_diagnostics(:GLOBAL, model) ==
                  (:E_kin, :Z, :P, :inj_E_kin, :inj_Z, :dissS_E_kin, :dissS_Z, :dissL_E_kin, :dissL_Z)
            @test default_diagnostics(:SPECTRAL, model) == (:E_kin_inc_k, :Z_k, :Pi_kin, :Pi_Z)
            for class in (:GLOBAL, :SPECTRAL)
                names = default_diagnostics(class, model)
                @test first.(select_diagnostics(class, names)) == collect(names)
            end
            # Spectra are reported at the shell centres, fluxes at the upper shell edges.
            spectral = Dict(select_diagnostics(:SPECTRAL, default_diagnostics(:SPECTRAL, model)))
            @test spectral[:E_kin_inc_k].coordinate == spectral[:Z_k].coordinate == :shell_k
            @test spectral[:Pi_kin].coordinate == spectral[:Pi_Z].coordinate == :shell_edge
            @test_throws ArgumentError select_diagnostics(:GLOBAL, (:E_kin, :no_such_name))
            @test_throws ArgumentError select_diagnostics(:FIELDS, (:E_kin,))
            @test_throws ArgumentError register_diagnostic!(:FIELDS, :x, (s, p, w) -> 0)
            @test_throws ArgumentError register_diagnostic!(:SPECTRAL, :x, (s, p, w) -> 0; coordinate = :k)

            # A new diagnostic is one function and one registration line (ADR-0009).
            register_diagnostic!(:GLOBAL, :test_twice_E_kin, (s, p, w) -> 2 * evaluate(:GLOBAL, :E_kin, s, p, w))
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, backend)
            state = State(kolmogorov(g; A = 2.0, n = 3))
            ws = Workspace(g, state)
            problem = Problem(g, model, NoForcing())
            @test evaluate(:GLOBAL, :test_twice_E_kin, state, problem, ws) ≈ 2.0
        end

        @testset "single mode lands in its shell with its value" for T in (Float64, Float32)
            # ψ = sin(3x + 4y), u = (−∂yψ, ∂xψ): |k| = 5, so E = |k|²/4 = 25/4 in shell 5, Z = |k|²E,
            # P = |k|²Z, and the rest of the spectrum is zero.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T, backend)
            ws = Workspace(g, State(kolmogorov(g; A = 1.0, n = 1)))
            u = spectral_velocity(g, ws, (x, y) -> -4cos(3x + 4y), (x, y) -> 3cos(3x + 4y))
            state = State(u)
            ws = Workspace(g, state)
            problem = Problem(g, Incompressible(ν = 0.1), NoForcing())
            E = evaluate(:SPECTRAL, :E_kin_inc_k, state, problem, ws)
            @test E == @inferred kinetic_energy_spectrum_inc(state, problem, ws)
            Zk = evaluate(:SPECTRAL, :Z_k, state, problem, ws)
            @test E isa Vector{T} && length(E) == length(g.shell_k)
            @test E[6] ≈ T(25 / 4) && Zk[6] ≈ T(625 / 4)
            @test maximum(abs, E[[1:5; 7:end]]) < 10eps(T) && maximum(abs, Zk[[1:5; 7:end]]) < 1000eps(T)
            @test evaluate(:GLOBAL, :E_kin, state, problem, ws) ≈ T(25 / 4)
            @test evaluate(:GLOBAL, :Z, state, problem, ws) ≈ T(625 / 4)
            @test evaluate(:GLOBAL, :P, state, problem, ws) ≈ T(15625 / 4)
        end

        @testset "shells use physical wavenumbers on a rectangle" begin
            # Lx = 4π, Ly = 2π: Δk = 1/2, and uy = cos(3x/2) has |k| = 3/2, shell 3, weight 2 (kx ≠ 0).
            g = Grid(; Lx = 4π, Ly = 2π, Nx = 32, Ny = 16, backend)
            ws = Workspace(g, State(kolmogorov(g; A = 1.0, n = 1)))
            state = State(spectral_velocity(g, ws, (x, y) -> 0.0, (x, y) -> cos(3x / 2)))
            ws = Workspace(g, state)
            E = evaluate(:SPECTRAL, :E_kin_inc_k, state, Problem(g, Incompressible(ν = 0.1), NoForcing()), ws)
            @test E[4] ≈ 1 / 4 && sum(E) ≈ 1 / 4
        end

        @testset "shell sums equal E_kin and Z of the physical field" for T in (Float64, Float32)
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            state = State(random_velocity(g; spectrum = k -> k^2 * exp(-k / 4), E_kin = 0.7, seed = 11))
            ws = Workspace(g, state)
            problem = Problem(g, Incompressible(ν = 0.1), NoForcing())
            E_kin = evaluate(:GLOBAL, :E_kin, state, problem, ws)
            @test E_kin == @inferred kinetic_energy(state, problem, ws)
            Z = evaluate(:GLOBAL, :Z, state, problem, ws)
            @test E_kin isa T && Z isa T
            u = (physical(state.fields.ux, ws), physical(state.fields.uy, ws))
            ω = physical(spectral_curl(g, state.fields), ws)
            @test E_kin ≈ mean_dot(u, u) / 2 rtol = 10eps(T)
            @test Z ≈ mean_product(ω, ω) / 2 rtol = 10eps(T)
            @test sum(evaluate(:SPECTRAL, :E_kin_inc_k, state, problem, ws)) ≈ E_kin rtol = 10eps(T)
            @test sum(evaluate(:SPECTRAL, :Z_k, state, problem, ws)) ≈ Z rtol = 10eps(T)
        end

        @testset "dissipation rates equal the loss by the linear terms" for T in (Float64, Float32)
            # dissS_E_kin = 2νZ, dissS_Z = 2νP; with ν only or α only, the linear terms remove E_kin and Z at
            # the rate Σ w Re[û*·σû] (WP-5), which the test computes on the host.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            state = State(random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed = 12))
            ws = Workspace(g, state)
            û = map(Array, state.fields)
            k2 = Array(g.k2)
            w = [c == 1 ? 1 : 2 for c in 1:size(k2, 1)]
            # The host reference sums in Float64, so its own round-off stays below the Float32 tolerance.
            linear_loss(σ, weight) = T(sum(w .* Float64.(weight) .* Float64.(σ) .* (abs2.(ComplexF64.(û.ux)) .+ abs2.(ComplexF64.(û.uy)))))
            for (model, small, large) in ((Incompressible(ν = 0.3), true, false),
                                          (Incompressible(ν = 0, large_scale = LargeScaleDissipation(α = 0.2, m = 1)),
                                           false, true))
                problem = Problem(g, model, NoForcing())
                σ = Array(problem.rates.velocity)
                d = Dict(name => evaluate(:GLOBAL, name, state, problem, ws)
                         for name in (:dissS_E_kin, :dissS_Z, :dissL_E_kin, :dissL_Z, :Z, :P))
                S, L = small ? (:dissS_E_kin, :dissS_Z) : (:dissL_E_kin, :dissL_Z)
                @test d[S] ≈ linear_loss(σ, 1) rtol = 100eps(T)
                @test d[L] ≈ linear_loss(σ, k2) rtol = 100eps(T)
                # The other pair is zero.
                O = small ? (:dissL_E_kin, :dissL_Z) : (:dissS_E_kin, :dissS_Z)
                @test d[O[1]] == 0 && d[O[2]] == 0
            end
            problem = Problem(g, Incompressible(ν = 0.3), NoForcing())
            @test evaluate(:GLOBAL, :dissS_E_kin, state, problem, ws) ≈ 2 * T(0.3) * evaluate(:GLOBAL, :Z, state, problem, ws)
            @test evaluate(:GLOBAL, :dissS_Z, state, problem, ws) ≈ 2 * T(0.3) * evaluate(:GLOBAL, :P, state, problem, ws)
        end

        @testset "fluxes follow the definition and vanish at k_max" for T in (Float64, Float32)
            # Π(K_n) = Σ_{shells ≤ n} w Re[â*·N̂] with ∂t û = −N̂: rhs! without forcing gives −N̂ (WP-6), and
            # the host sums it per shell by grid.shell_index, independent of the device reduction.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 64, Ny = 64, T, backend)
            state = State(random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed = 13))
            ws = Workspace(g, state)
            problem = Problem(g, Incompressible(ν = 0.1), NoForcing())
            Pi_kin = evaluate(:SPECTRAL, :Pi_kin, state, problem, ws)
            @test Pi_kin == @inferred kinetic_energy_flux(state, problem, ws)
            Pi_Z = evaluate(:SPECTRAL, :Pi_Z, state, problem, ws)
            @test Pi_kin isa Vector{T} && Pi_Z isa Vector{T}

            dû = map(similar, state.fields)
            rhs!(dû, state, problem, ws)
            û, n̂ = map(Array, state.fields), map(Array, dû)
            ω̂, ω̂_n = Array(spectral_curl(g, state.fields)), Array(spectral_curl(g, dû))
            shell = Array(g.shell_index)
            energy, enstrophy, scale = zeros(length(g.shell_k)), zeros(length(g.shell_k)), zeros(2)
            for I in CartesianIndices(û.ux)
                w = I[1] == 1 ? 1 : 2
                e = -w * real(conj(û.ux[I]) * n̂.ux[I] + conj(û.uy[I]) * n̂.uy[I])
                z = -w * real(conj(ω̂[I]) * ω̂_n[I])
                energy[shell[I]] += e
                enstrophy[shell[I]] += z
                scale .+= abs.((e, z))
            end
            @test scale[1] > 1 && scale[2] > 1
            @test maximum(abs, Pi_kin - cumsum(energy)) < 100eps(T) * scale[1]
            @test maximum(abs, Pi_Z - cumsum(enstrophy)) < 100eps(T) * scale[2]
            # Blocking (specification, Galerkin truncation): the fluxes vanish at k_max and beyond.
            n_max = findfirst(>(g.k_max), g.shell_edge)
            @test maximum(abs, Pi_kin[n_max:end]) < 100eps(T) * scale[1]
            @test maximum(abs, Pi_Z[n_max:end]) < 100eps(T) * scale[2]
            # The flux is not trivially zero below k_max.
            @test maximum(abs, Pi_kin) > 1000eps(T) * scale[1]
        end

        @testset "injection rates on a hand-built û and f̂" for T in (Float64, Float32)
            # u = A sin(3y) x̂ and f = F sin(3y) x̂: ⟨u·f⟩ + ½Δt⟨|f|²⟩ = AF/2 + ΔtF²/4, and inj_Z = 9 inj_E_kin.
            A, F, Δt = T(0.8), T(1.5), T(0.1)
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T, backend)
            state = State(kolmogorov(g; A, n = 3))
            ws = Workspace(g, state)
            forcing = KolmogorovForcing(g; F, n = 3)
            problem = Problem(g, Incompressible(ν = 0.1), forcing)
            # Stored rates start as NaN; the default method for any other problem changes nothing.
            @test isnan(evaluate(:GLOBAL, :inj_E_kin, state, problem, ws))
            @test isnan(evaluate(:GLOBAL, :inj_Z, state, problem, ws))
            @test record_injection!(:other_problem, state, ws, Δt) === ws
            @test isnan(evaluate(:GLOBAL, :inj_E_kin, state, problem, ws))

            update_forcing!(forcing, state, ws, Δt)
            @test (@inferred record_injection!(problem, state, ws, Δt)) === ws
            inj_E = evaluate(:GLOBAL, :inj_E_kin, state, problem, ws)
            @test inj_E isa T
            @test inj_E ≈ A * F / 2 + Δt * F^2 / 4
            @test evaluate(:GLOBAL, :inj_Z, state, problem, ws) ≈ 9 * (A * F / 2 + Δt * F^2 / 4)
        end

        @testset "injection rates equal the physical-space means for a band forcing" for T in (Float64, Float32)
            # The band holds modes with kx = 0 (weight 1) and kx ≠ 0 (weight 2); the grid means
            # ⟨u·f⟩ + ½Δt⟨|f|²⟩ and ⟨ω curl f⟩ + ½Δt⟨(curl f)²⟩ need no weights.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            state = State(random_velocity(g; spectrum = k -> one(k), E_kin = 1, seed = 14))
            ws = Workspace(g, state)
            forcing = OUForcing(g, band(3, 5); ε = 1, τ_f = 0.05, seed = 2)
            problem = Problem(g, Incompressible(ν = 0.1), forcing)
            Δt = T(0.01)
            update_forcing!(forcing, state, ws, Δt)
            record_injection!(problem, state, ws, Δt)
            u = (physical(state.fields.ux, ws), physical(state.fields.uy, ws))
            f = (physical(ws.forcing.ux, ws), physical(ws.forcing.uy, ws))
            ω, curl_f = physical(spectral_curl(g, state.fields), ws), physical(spectral_curl(g, ws.forcing), ws)
            expected_E = mean_dot(u, f) + Δt / 2 * mean_dot(f, f)
            expected_Z = mean_product(ω, curl_f) + Δt / 2 * mean_product(curl_f, curl_f)
            scale_E = sqrt(mean_dot(u, u) * mean_dot(f, f))
            scale_Z = sqrt(mean_product(ω, ω) * mean_product(curl_f, curl_f))
            @test abs(evaluate(:GLOBAL, :inj_E_kin, state, problem, ws) - expected_E) < 100eps(T) * scale_E
            @test abs(evaluate(:GLOBAL, :inj_Z, state, problem, ws) - expected_Z) < 100eps(T) * scale_Z
        end
    end
end

end
