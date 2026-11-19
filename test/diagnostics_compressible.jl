# Tests of the compressible diagnostics (WP-11): registry defaults, hand values on single modes, GLOBAL values
# against grid means of analytic fields, spectra and fluxes against their definitions, equality with the
# incompressible diagnostics at ρ = 1, the injection rates of a forcing-only run and the budgets of a nonlinear
# run through the production step!.
module Diagnostics_compressibleTests

using Test
using KernelAbstractions
using CompReS
using CompReS: update_forcing!

# Host grid coordinates of g, as Float64 matrices X[i, j] = x_i, Y[i, j] = y_j.
function coordinates(g)
    x = [g.Lx * (i - 1) / g.Nx for i in 1:g.Nx]
    y = [g.Ly * (j - 1) / g.Ny for j in 1:g.Ny]
    return Float64.(x .+ 0 .* y'), Float64.(0 .* x .+ y')
end

# Device spectral field of the host physical values f, in the float type of g.
function spectral(g, ws, f)
    T = eltype(g.k2)
    return to_spectral!(KernelAbstractions.allocate(g.backend, Complex{T}, g.Nx ÷ 2 + 1, g.Ny),
                        copyto!(similar(ws.physical.ux), T.(f)), ws.transforms)
end

# Physical host field of the spectral field f̂.
physical(f̂, ws) = Array(to_physical!(similar(ws.physical.ux), f̂, ws.transforms))

# Compressible state of the host physical values (ux, uy, rho), with its workspace.
function state_of(g, ux, uy, rho; kwargs...)
    T = eltype(g.k2)
    ws0 = Workspace(g, State(merge(kolmogorov(g; A = 1, n = 1), uniform_density(g))))
    state = State((ux = spectral(g, ws0, ux), uy = spectral(g, ws0, uy), rho = spectral(g, ws0, rho)))
    return state, Workspace(g, state; kwargs...)
end

# Value of the GLOBAL or SPECTRAL diagnostic `name` for the state and problem.
evaluate(class, name, state, problem, ws) = only(select_diagnostics(class, (name,))).second(state, problem, ws)

# Grid mean of a host array.
gmean(a) = sum(a) / length(a)

# Host sums per shell of the per-mode values d over all grid modes, with the half-plane weight
# (1 on the kx = 0 column and the Nyquist column, 2 elsewhere).
function shell_reference(g, d)
    nkx = size(d, 1)
    shell = Array(g.shell_index)
    out = zeros(length(g.shell_k))
    for I in CartesianIndices(d)
        out[shell[I]] += (I[1] == 1 || I[1] == nkx ? 1 : 2) * d[I]
    end
    return out
end

# Host Helmholtz parts (solenoidal, compressive) of the spectral vector field (ax, ay), with the derivative
# wavenumbers of the transforms (helmholtz_split! convention).
function host_helmholtz(ws, ax, ay)
    kx, ky = Array(ws.transforms.kx_deriv), Array(ws.transforms.ky_deriv)'
    k2 = kx .^ 2 .+ ky .^ 2
    d = ifelse.(k2 .== 0, 0, (kx .* ax .+ ky .* ay) ./ ifelse.(k2 .== 0, 1, k2))
    return (ax .- kx .* d, ay .- ky .* d), (kx .* d, ky .* d)
end

dot_re(a, b) = real.(conj.(a[1]) .* b[1] .+ conj.(a[2]) .* b[2])

# Analytic test fields with low modes: u has a solenoidal and a compressive part, ρ = 1 + δ with two modes.
# Every product of two of these fields stays below the Nyquist wavenumber of N = 32, so spectral derivatives
# of the products are exact on the grid.
const a1, a2, a3, b1, b2 = 0.7, 0.5, 0.3, 0.4, 0.6
const r1, r2 = 0.1, 0.05
ux_f(x, y) = a1 * sin(3y) + a2 * sin(2x) + a3 * cos(x + y)
uy_f(x, y) = b1 * cos(2x - y) + b2 * sin(x)
ρ_f(x, y) = 1 + r1 * cos(x + 2y) + r2 * sin(2x)
dxux(x, y) = 2a2 * cos(2x) - a3 * sin(x + y)
dyux(x, y) = 3a1 * cos(3y) - a3 * sin(x + y)
dxuy(x, y) = -2b1 * sin(2x - y) + b2 * cos(x)
dyuy(x, y) = b1 * sin(2x - y)
dxρ(x, y) = -r1 * sin(x + 2y) + 2r2 * cos(2x)
dyρ(x, y) = -2r1 * sin(x + 2y)
lapρ(x, y) = -5r1 * cos(x + 2y) - 4r2 * sin(2x)
ω_f(x, y) = dxuy(x, y) - dyux(x, y)
θ_f(x, y) = dxux(x, y) + dyuy(x, y)
# (−∇²)^{-1} of each velocity component and of ω, mode by mode.
Lux(x, y) = a1 * sin(3y) / 9 + a2 * sin(2x) / 4 + a3 * cos(x + y) / 2
Luy(x, y) = b1 * cos(2x - y) / 5 + b2 * sin(x)
Lω(x, y) = -2b1 * sin(2x - y) / 5 + b2 * cos(x) - 3a1 * cos(3y) / 9 + a3 * sin(x + y) / 2
# Gradients of ω and θ.
dxω(x, y) = -4b1 * cos(2x - y) - b2 * sin(x) + a3 * cos(x + y)
dyω(x, y) = 2b1 * cos(2x - y) + 9a1 * sin(3y) + a3 * cos(x + y)
dxθ(x, y) = -4a2 * sin(2x) - a3 * cos(x + y) + 2b1 * cos(2x - y)
dyθ(x, y) = -a3 * cos(x + y) - b1 * cos(2x - y)

const GLOBAL_DEFAULTS = (:mass, :E_kin, :E_kin_inc, :E_kin_comp, :E_int, :E_tot, :Z, :Omega, :u_rms, :Ma, :p_theta,
                         :inj_E_kin, :inj_Omega, :dissS_E_kin, :dissS_Omega, :dissL_E_kin, :dissL_Omega,
                         :kappa_E_kin, :kappa_E_int, :kappa_Omega)
const SPECTRAL_DEFAULTS = (:E_kin_inc_k, :E_kin_comp_k, :Omega_k, :rho_k, :Pi_kin, :Pi_kin_inc, :Pi_kin_comp,
                           :Pi_Omega)
const SCHEMES = (Euler(), RK2(), SSPRK3(), RK4())

function run_tests(backend)
    @testset "compressible diagnostics" begin
        @testset "registry: compressible defaults" begin
            model = Compressible(; eos = Isothermal(cs = 1), μ = 0.1, ζ = 0)
            @test default_diagnostics(:GLOBAL, model) == GLOBAL_DEFAULTS
            @test default_diagnostics(:SPECTRAL, model) == SPECTRAL_DEFAULTS
            @test_throws ArgumentError default_diagnostics(:FIELDS, model)
            spectral_d = Dict(select_diagnostics(:SPECTRAL, SPECTRAL_DEFAULTS))
            @test all(spectral_d[n].coordinate == :shell_k for n in SPECTRAL_DEFAULTS[1:4])
            @test all(spectral_d[n].coordinate == :shell_edge for n in SPECTRAL_DEFAULTS[5:8])
            @test first.(select_diagnostics(:GLOBAL, GLOBAL_DEFAULTS)) == collect(GLOBAL_DEFAULTS)
        end

        @testset "single modes with uniform density: hand values, $T" for T in (Float64, Float32)
            # ρ = ρ₀ = 2, u = (A sin 3y + B sin 2x) x̂: A is solenoidal in shell 3, B compressive in shell 2.
            # ω = −3A cos 3y, θ = 2B cos 2x; Isothermal c_s = 3; m = 1, so (−∇²)^{-1} divides by k².
            A, B, ρ₀, cs, μ, ζ, α = 0.8, 0.6, 2.0, 3.0, 0.1, 0.05, 0.2
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T, backend)
            X, Y = coordinates(g)
            state, ws = state_of(g, A .* sin.(3Y) .+ B .* sin.(2X), zero(X), fill(ρ₀, size(X)))
            model = Compressible(; eos = Isothermal(; cs), μ, ζ, κ = 0.3, large_scale = LargeScaleDissipation(; α, m = 1))
            problem = Problem(g, model, NoForcing())
            v(name) = evaluate(:GLOBAL, name, state, problem, ws)
            s(name) = evaluate(:SPECTRAL, name, state, problem, ws)
            tol = 100eps(T)
            @test v(:mass) isa T
            @test v(:mass) ≈ ρ₀ rtol = tol
            @test v(:E_kin) ≈ ρ₀ * (A^2 + B^2) / 4 rtol = tol
            @test v(:E_kin_inc) ≈ ρ₀ * A^2 / 4 rtol = tol
            @test v(:E_kin_comp) ≈ ρ₀ * B^2 / 4 rtol = tol
            @test v(:E_int) ≈ ρ₀ * cs^2 * log(ρ₀) rtol = tol
            @test v(:E_tot) ≈ ρ₀ * (A^2 + B^2) / 4 + ρ₀ * cs^2 * log(ρ₀) rtol = tol
            @test v(:Z) ≈ 9A^2 / 4 rtol = tol
            @test v(:Omega) ≈ 9A^2 / (4ρ₀) rtol = tol
            @test v(:u_rms) ≈ sqrt((A^2 + B^2) / 2) rtol = tol
            @test v(:Ma) ≈ sqrt((A^2 + B^2) / 2) / cs rtol = tol
            @test abs(v(:p_theta)) < tol * cs^2 * ρ₀ * B
            @test v(:dissS_E_kin) ≈ μ * (9A^2 / 2 + 2B^2) + 2ζ * B^2 rtol = tol
            @test v(:dissS_Omega) ≈ μ * 81A^2 / (2ρ₀^2) rtol = tol
            @test v(:dissL_E_kin) ≈ ρ₀ * α * (A^2 / 18 + B^2 / 8) rtol = tol
            @test v(:dissL_Omega) ≈ α * A^2 / (2ρ₀) rtol = tol
            for name in (:kappa_E_kin, :kappa_E_int, :kappa_Omega)
                @test abs(v(name)) < tol
            end
            # Spectra: each mode in its shell (entry n + 1 is shell n), nothing elsewhere.
            E_inc, E_comp, Ωk, ρk = s(:E_kin_inc_k), s(:E_kin_comp_k), s(:Omega_k), s(:rho_k)
            @test E_inc isa Vector{T} && length(E_inc) == length(g.shell_k)
            @test E_inc[4] ≈ ρ₀ * A^2 / 4 rtol = tol
            @test E_comp[3] ≈ ρ₀ * B^2 / 4 rtol = tol
            @test Ωk[4] ≈ 9A^2 / (4ρ₀) rtol = tol
            @test ρk[1] ≈ ρ₀^2 / 2 rtol = tol
            @test maximum(abs, E_inc[[1:3; 5:end]]) < tol && maximum(abs, E_comp[[1:2; 4:end]]) < tol
            @test maximum(abs, Ωk[[1:3; 5:end]]) < 10tol && maximum(abs, ρk[2:end]) < tol
        end

        @testset "GLOBAL values equal grid means of the physical fields, $(nameof(typeof(eos))), $T" for eos in
                (Isothermal(cs = 3), Polytropic(K = 2, γ = 1.4)), T in (Float64, Float32)
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            X, Y = coordinates(g)
            f(fn) = fn.(X, Y)
            ux, uy, ρ, ω, θ, ∇²ρ = f(ux_f), f(uy_f), f(ρ_f), f(ω_f), f(θ_f), f(lapρ)
            state, ws = state_of(g, ux, uy, ρ)
            μ, ζ, κ, α = 0.1, 0.05, 0.03, 0.2
            model = Compressible(; eos, μ, ζ, κ, large_scale = LargeScaleDissipation(; α, m = 1))
            problem = Problem(g, model, NoForcing())
            v(name) = evaluate(:GLOBAL, name, state, problem, ws)
            e64 = CompReS.convert_float(Float64, eos)
            u2 = ux .^ 2 .+ uy .^ 2
            # Each reference is the grid mean of an integrand; the tolerance scales with the mean of its modulus.
            integrands = Dict(:mass => ρ, :E_kin => ρ .* u2 ./ 2, :E_int => ρ .* internal_energy.(Ref(e64), ρ),
                              :Z => ω .^ 2 ./ 2, :Omega => ω .^ 2 ./ (2ρ), :p_theta => pressure.(Ref(e64), ρ) .* θ,
                              :dissS_E_kin => μ .* (f(dxux) .^ 2 .+ f(dyux) .^ 2 .+ f(dxuy) .^ 2 .+ f(dyuy) .^ 2) .+
                                              ζ .* θ .^ 2,
                              :dissL_E_kin => α .* ρ .* (ux .* f(Lux) .+ uy .* f(Luy)),
                              :dissL_Omega => α .* ω .* f(Lω) ./ ρ,
                              :kappa_E_kin => κ .* u2 .* ∇²ρ ./ 2,
                              :kappa_E_int => κ .* enthalpy.(Ref(e64), ρ) .* ∇²ρ,
                              :kappa_Omega => -κ .* ω .^ 2 .* ∇²ρ ./ (2 .* ρ .^ 2))
            for (name, a) in integrands
                @test v(name) isa T
                @test abs(v(name) - gmean(a)) < 100eps(T) * gmean(abs.(a))
            end
            E_kin, E_int = gmean(integrands[:E_kin]), gmean(integrands[:E_int])
            @test v(:E_tot) ≈ E_kin + E_int rtol = 100eps(T)
            @test v(:u_rms) ≈ sqrt(gmean(u2)) rtol = 100eps(T)
            @test v(:Ma) ≈ sqrt(gmean(u2) / sound_speed²(e64, 1.0)) rtol = 100eps(T)
            @test v(:E_kin_inc) + v(:E_kin_comp) ≈ E_kin rtol = 100eps(T)
            @test v(:E_kin_comp) > 0.01 * E_kin
            # dissS_Omega = −⟨q curl g⟩ with g = (μ∇²u + ζ∇θ)/ρ, from the analytic gradients. The code takes the
            # curl of the grid values of g spectrally: 1/ρ is not band-limited, so its aliased harmonics leave a
            # relative difference of order max(r1, r2)^8 ≈ 1e-8 at N = 32.
            Vx = μ .* (-9a1 .* sin.(3Y) .- 4a2 .* sin.(2X) .- 2a3 .* cos.(X .+ Y)) .+ ζ .* f(dxθ)
            Vy = μ .* (-5b1 .* cos.(2X .- Y) .- b2 .* sin.(X)) .+ ζ .* f(dyθ)
            curlV = μ .* (10b1 .* sin.(2X .- Y) .- b2 .* cos.(X) .+ 27a1 .* cos.(3Y) .- 2a3 .* sin.(X .+ Y))  # μ∇²ω
            curl_g = curlV ./ ρ .- (f(dxρ) .* Vy .- f(dyρ) .* Vx) ./ ρ .^ 2
            dissS_Ω = -gmean(ω ./ ρ .* curl_g)
            @test abs(v(:dissS_Omega) - dissS_Ω) < max(100eps(T), 1e-7) * abs(dissS_Ω)
        end

        @testset "spectra and fluxes follow their definitions, $T" for T in (Float64, Float32)
            # Reference: ŵ, ŝ and the nonlinear terms N^w, N^s of the specification in advective form, sampled on
            # the grid and transformed, then summed per shell on the host. The code uses the Lamb form.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            X, Y = coordinates(g)
            f(fn) = fn.(X, Y)
            ux, uy, ρ, ω = f(ux_f), f(uy_f), f(ρ_f), f(ω_f)
            state, ws = state_of(g, ux, uy, ρ)
            problem = Problem(g, Compressible(; eos = Isothermal(cs = 3), μ = 0.1, ζ = 0.05), NoForcing())
            s(name) = evaluate(:SPECTRAL, name, state, problem, ws)
            divρu = f(dxρ) .* ux .+ f(dyρ) .* uy .+ ρ .* f(θ_f)
            advx = ux .* f(dxux) .+ uy .* f(dyux)
            advy = ux .* f(dxuy) .+ uy .* f(dyuy)
            Nwx = sqrt.(ρ) .* advx .+ ux .* divρu ./ (2 .* sqrt.(ρ))
            Nwy = sqrt.(ρ) .* advy .+ uy .* divρu ./ (2 .* sqrt.(ρ))
            div_uω = ux .* f(dxω) .+ uy .* f(dyω) .+ ω .* f(θ_f)
            Ns = div_uω ./ sqrt.(ρ) .- ω .* divρu ./ (2 .* ρ .^ 1.5)
            host(a) = ComplexF64.(Array(spectral(g, ws, a)))
            ŵ = (host(sqrt.(ρ) .* ux), host(sqrt.(ρ) .* uy))
            N̂w = (host(Nwx), host(Nwy))
            ŝ, N̂s = host(ω ./ sqrt.(ρ)), host(Ns)
            ŵs, ŵc = host_helmholtz(ws, ŵ...)
            # Per-mode values d and their round-off scale m = |a||b| for the product of the modes a and b.
            norm2(a) = sqrt.(abs2.(a[1]) .+ abs2.(a[2]))
            modes = Dict(:E_kin_inc_k => (dot_re(ŵs, ŵs) ./ 2, norm2(ŵ) .^ 2),
                         :E_kin_comp_k => (dot_re(ŵc, ŵc) ./ 2, norm2(ŵ) .^ 2),
                         :Omega_k => (abs2.(ŝ) ./ 2, abs2.(ŝ)), :rho_k => (abs2.(host(ρ)) ./ 2, abs2.(host(ρ))),
                         :Pi_kin => (dot_re(ŵ, N̂w), norm2(ŵ) .* norm2(N̂w)),
                         :Pi_kin_inc => (dot_re(ŵs, N̂w), norm2(ŵ) .* norm2(N̂w)),
                         :Pi_kin_comp => (dot_re(ŵc, N̂w), norm2(ŵ) .* norm2(N̂w)),
                         :Pi_Omega => (real.(conj.(ŝ) .* N̂s), abs.(ŝ) .* abs.(N̂s)))
            # The tolerance scales with Σ m over all modes.
            for (name, (d, m)) in modes
                expected = name in (:Pi_kin, :Pi_kin_inc, :Pi_kin_comp, :Pi_Omega) ? cumsum(shell_reference(g, d)) :
                           shell_reference(g, d)
                scale = sum(shell_reference(g, m))
                value = s(name)
                @test value isa Vector{T} && length(value) == length(g.shell_k)
                @test maximum(abs, expected) > 1000eps(T) * scale
                @test maximum(abs, value .- expected) < 100eps(T) * scale
            end
            # The Helmholtz parts add up, and the spectra sum to the GLOBAL values.
            @test maximum(abs, s(:Pi_kin_inc) .+ s(:Pi_kin_comp) .- s(:Pi_kin)) < 100eps(T) * maximum(abs, s(:Pi_kin))
            v(name) = evaluate(:GLOBAL, name, state, problem, ws)
            @test sum(s(:E_kin_inc_k)) ≈ v(:E_kin_inc) rtol = 100eps(T)
            @test sum(s(:E_kin_comp_k)) ≈ v(:E_kin_comp) rtol = 100eps(T)
            @test sum(s(:Omega_k)) ≈ v(:Omega) rtol = 100eps(T)
            @test sum(s(:rho_k)) ≈ gmean(ρ .^ 2) / 2 rtol = 100eps(T)
        end

        @testset "ρ = 1: compressible diagnostics equal the incompressible ones, $T" for T in (Float64, Float32)
            # Solenoidal truncated u, ρ = 1, μ = ν, the same large-scale dissipation and the same solenoidal OU
            # forcing: E_kin, Z, Omega = Z, the spectra, the fluxes and every rate agree (specification: for
            # ρ = ρ₀, Π_K = ρ₀Π_E and Π_Ω = Π_Z/ρ₀). The two states share the velocity arrays.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            u = random_velocity(g; spectrum = k -> k^2 * exp(-k / 4), E_kin = 0.7, seed = 21)
            si = State(u)
            sc = State(merge(u, uniform_density(g)))
            wi, wc = Workspace(g, si), Workspace(g, sc)
            ls = LargeScaleDissipation(α = 0.2, m = 1)
            forcing = OUForcing(g, band(3, 5); ε = 1, τ_f = 0.05, seed = 4)
            pi_ = Problem(g, Incompressible(; ν = 0.1, large_scale = ls), forcing)
            pc = Problem(g, Compressible(; eos = Isothermal(cs = 5), μ = 0.1, ζ = 0.07, κ = 0.02, large_scale = ls), forcing)
            Δt = T(0.01)
            update_forcing!(forcing, sc, wc, Δt)
            foreach(copyto!, values(wi.forcing), values(wc.forcing))
            record_injection!(pi_, si, wi, Δt)
            record_injection!(pc, sc, wc, Δt)
            gi(name) = evaluate(:GLOBAL, name, si, pi_, wi)
            gc(name) = evaluate(:GLOBAL, name, sc, pc, wc)
            for (c, i) in ((:E_kin, :E_kin), (:E_kin_inc, :E_kin), (:Z, :Z), (:Omega, :Z), (:inj_E_kin, :inj_E_kin),
                           (:inj_Omega, :inj_Z), (:dissS_E_kin, :dissS_E_kin), (:dissS_Omega, :dissS_Z),
                           (:dissL_E_kin, :dissL_E_kin), (:dissL_Omega, :dissL_Z))
                @test gc(c) isa T
                @test gc(c) ≈ gi(i) rtol = 100eps(T)
            end
            @test gc(:mass) == 1
            for name in (:E_kin_comp, :E_int, :p_theta, :kappa_E_kin, :kappa_E_int, :kappa_Omega)
                @test abs(gc(name)) < 100eps(T) * gc(:E_kin) * 25
            end
            si_(name) = evaluate(:SPECTRAL, name, si, pi_, wi)
            sc_(name) = evaluate(:SPECTRAL, name, sc, pc, wc)
            for (c, i) in ((:E_kin_inc_k, :E_kin_inc_k), (:Omega_k, :Z_k), (:Pi_kin, :Pi_kin), (:Pi_kin_inc, :Pi_kin),
                           (:Pi_Omega, :Pi_Z))
                expected = si_(i)
                @test maximum(abs, sc_(c) .- expected) < 100eps(T) * maximum(abs, expected)
            end
            @test maximum(abs, sc_(:E_kin_comp_k)) < 100eps(T) * gc(:E_kin)
            @test maximum(abs, sc_(:Pi_kin_comp)) < 100eps(T) * maximum(abs, si_(:Pi_kin))
        end

        @testset "forcing-only injection through step!: $(nameof(typeof(scheme))), IF $integrating_factor, $T" for scheme in SCHEMES,
                                                                                                             integrating_factor in (false, true),
                                                                                                             T in (Float64, Float32)
            # WP-11 set-up: with forcing only and ρ fixed, every scheme gives u_{n+1} = u_n + Δt f_n, so
            # inj_E_kin = ΔE_kin/Δt and inj_Omega = ΔOmega/Δt up to round-off. ρ couples the modes of the band
            # that differ by ±(1, 1), so a formula without ρ fails by several percent.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T, backend)
            Δt = 1e-3
            forcing = OUForcing(g, band(4, 6); ε = 1, τ_f = 0.05, seed = 1, compressive_fraction = 0.5)
            model = Compressible(; eos = Isothermal(cs = 1), μ = 0, ζ = 0, κ = 0, nonlinear = false)
            problem = Problem(g, model, forcing)
            u = random_velocity(g; spectrum = k -> k^4 * exp(-k^2 / 8), E_kin = 0.5, seed = 2)
            X, Y = coordinates(g)
            ws0 = Workspace(g, State(u))
            state = State(merge(u, (rho = spectral(g, ws0, 1 .+ 0.1 .* cos.(X .+ Y)),)))
            ρ0 = Array(state.fields.rho)
            ws = Workspace(g, state; scheme, integrating_factor)
            v(name) = evaluate(:GLOBAL, name, state, problem, ws)
            E, Ω = v(:E_kin), v(:Omega)
            worst_E, worst_Ω = zero(T), zero(T)
            for _ in 1:100
                step!(state, problem, ws, Δt)
                E1, Ω1 = v(:E_kin), v(:Omega)
                worst_E = max(worst_E, abs(v(:inj_E_kin) - (E1 - E) / T(Δt)) / (eps(T) * E1 / T(Δt)))
                worst_Ω = max(worst_Ω, abs(v(:inj_Omega) - (Ω1 - Ω) / T(Δt)) / (eps(T) * Ω1 / T(Δt)))
                E, Ω = E1, Ω1
            end
            @test state.n == 100
            @test Array(state.fields.rho) == ρ0
            @test worst_E ≤ 10
            @test worst_Ω ≤ 10
        end

        @testset "nonlinear run through step!: mass and global budgets, IF $integrating_factor" for integrating_factor in (false, true)
            # Unforced nonlinear run with every dissipation term on. The specification (compressible model,
            # global budgets) gives
            #   dE_kin/dt = p_theta − dissS_E_kin − dissL_E_kin + kappa_E_kin,
            #   dE_int/dt = −p_theta + kappa_E_int,
            #   dOmega/dt = −dissS_Omega − dissL_Omega + kappa_Omega.
            # The trapezoidal rule integrates the rates sampled at every step; its error is O(Δt²). The truncation
            # and aliasing leave a floor that does not scale with Δt (specification, note 2): about 2e-4 of the
            # Omega scale at N = 32, below 3e-7 at N = 64. At N = 64 and Δt = 1e-3 the measured residuals are below
            # 1.3e-6 of their scale; halving Δt divides the E_kin and E_int residuals by about 4.
            # The mass mode has a zero tendency and a zero linear rate, so it does not change at all.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 64, Ny = 64, backend)
            model = Compressible(; eos = Isothermal(cs = 3), μ = 0.02, ζ = 0.01, κ = 0.01,
                                 large_scale = LargeScaleDissipation(α = 0.05, m = 0))
            problem = Problem(g, model, NoForcing())
            u = random_velocity(g; spectrum = k -> k^4 * exp(-k^2 / 2), E_kin = 0.5, seed = 5)
            X, Y = coordinates(g)
            ws0 = Workspace(g, State(u))
            state = State(merge(u, (rho = spectral(g, ws0, 1 .+ 0.05 .* cos.(X .+ Y) .+ 0.03 .* sin.(2Y)),)))
            ws = Workspace(g, state; scheme = RK4(), integrating_factor)
            v(name) = evaluate(:GLOBAL, name, state, problem, ws)
            rates() = (E_kin = v(:p_theta) - v(:dissS_E_kin) - v(:dissL_E_kin) + v(:kappa_E_kin),
                       E_int = -v(:p_theta) + v(:kappa_E_int),
                       Omega = -v(:dissS_Omega) - v(:dissL_Omega) + v(:kappa_Omega))
            sizes() = (E_kin = abs(v(:p_theta)) + v(:dissS_E_kin) + v(:dissL_E_kin) + abs(v(:kappa_E_kin)),
                       E_int = abs(v(:p_theta)) + abs(v(:kappa_E_int)),
                       Omega = abs(v(:dissS_Omega)) + abs(v(:dissL_Omega)) + abs(v(:kappa_Omega)))
            names = (:E_kin, :E_int, :Omega)
            start = map(v, NamedTuple{names}(names))
            mass0 = v(:mass)
            Δt, nsteps = 1e-3, 500
            integral, scale, r = map(zero, start), map(zero, start), rates()
            mass_change = 0.0
            for _ in 1:nsteps
                step!(state, problem, ws, Δt)
                r1 = rates()
                integral = map((i, a, b) -> i + Δt * (a + b) / 2, integral, r, r1)
                scale = map((s, a) -> s + Δt * a, scale, sizes())
                r = r1
                mass_change = max(mass_change, abs(v(:mass) - mass0))
            end
            @test mass_change == 0
            for name in names
                change = v(name) - start[name]
                @info "budget of $name" change integral[name] scale[name] relative_error = abs(change - integral[name]) / scale[name]
                @test abs(change - integral[name]) < 1e-5 * scale[name]
            end
        end
    end
end

end
