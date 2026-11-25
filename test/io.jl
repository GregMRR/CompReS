# Tests of the NetCDF output and the restart (WP-12, ADR-0009, ADR-0024): exact read-back of GLOBAL and SPECTRAL
# values, the restart error of a snapshot, the files' layout and parameters, and a split run against a straight run.
module IoTests

using Test
using KernelAbstractions
using NCDatasets: NCDataset
using TOML
using CompReS

const SPECTRUM = k -> k^4 * exp(-k^2 / 8)

# Grid, problem, state and workspace of a forced run of the incompressible model.
function setup(backend, T; compressible = false, N = 16)
    g = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N, T, backend)
    model = compressible ? Compressible(; eos = Isothermal(cs = 10), μ = 0.01, ζ = 0.005) : Incompressible(; ν = 0.01)
    forcing = OUForcing(g, band(2, 4); ε = 1.0, τ_f = 0.2, seed = 3)
    fields = random_velocity(g; spectrum = SPECTRUM, E_kin = 0.5, seed = 2)
    state = State(compressible ? merge(fields, uniform_density(g)) : fields)
    return g, Problem(g, model, forcing), state, Workspace(g, state)
end

function advance!(state, problem, ws, nsteps, Δt)
    for _ in 1:nsteps
        step!(state, problem, ws, Δt)
    end
    return state
end

relative_error(a, b) = sqrt(sum(abs2, Array(a) .- Array(b)) / sum(abs2, Array(b)))

# Relative 2-norm distance of two states over all their fields.
state_error(a, b) = sqrt(sum(name -> sum(abs2, Array(a.fields[name]) .- Array(b.fields[name])), keys(a.fields)) /
                         sum(name -> sum(abs2, Array(b.fields[name])), keys(b.fields)))

function run_tests(backend)
    @testset "io" begin
        @testset "Output and selection of diagnostics" begin
            @test Output(period = 2).period == 2.0
            @test_throws ArgumentError Output(period = 0)
            model = Incompressible(; ν = 0.01)
            names(o, class) = first.(selected_diagnostics(o, class, model))
            @test names(Output(period = 1), :GLOBAL) == collect(default_diagnostics(:GLOBAL, model))
            @test names(Output(period = 1, diagnostics = [:Z, :E_kin]), :GLOBAL) == [:Z, :E_kin]
            @test_throws ArgumentError selected_diagnostics(Output(period = 1, diagnostics = [:nope]), :GLOBAL, model)
        end

        @testset "GLOBAL and SPECTRAL files read back exactly: $T" for T in (Float64, Float32)
            g, problem, state, ws = setup(backend, T)
            global_d = selected_diagnostics(Output(period = 0.1), :GLOBAL, problem.model)
            spectral_d = selected_diagnostics(Output(period = 0.2), :SPECTRAL, problem.model)
            mktempdir() do dir
                gpath, spath = joinpath(dir, "global.nc"), joinpath(dir, "spectral.nc")
                gds = open_global(gpath, problem, global_d; period = 0.1, run = (; scheme = "RK4", Δt = 0.01))
                sds = open_spectral(spath, problem, spectral_d; period = 0.2)
                expected_global = Dict(name => T[] for (name, _) in global_d)
                expected_spectral = Dict(name => Vector{T}[] for (name, _) in spectral_d)
                times = T[]
                for sample in 1:3
                    sample > 1 && advance!(state, problem, ws, 2, 0.01)
                    write_global!(gds, state, problem, ws, global_d)
                    write_spectral!(sds, state, problem, ws, spectral_d)
                    push!(times, state.t)
                    for (name, d) in global_d
                        push!(expected_global[name], d(state, problem, ws))
                    end
                    for (name, d) in spectral_d
                        push!(expected_spectral[name], d(state, problem, ws))
                    end
                end
                close(gds)
                close(sds)
                NCDataset(gpath) do ds
                    @test ds.dim["time"] == 3
                    @test isequal(Array(ds["time"][:]), times)
                    @test Array(ds["step"][:]) == [0, 2, 4]
                    for (name, _) in global_d
                        @test eltype(ds[String(name)]) == T
                        @test isequal(Array(ds[String(name)][:]), expected_global[name])
                    end
                    # The injection rates are NaN at t₀, stored as plain values.
                    @test isnan(ds["inj_E_kin"][1]) && !isnan(ds["inj_E_kin"][3])
                    @test !haskey(ds["inj_E_kin"].attrib, "_FillValue")
                    @test ds.attrib["period"] == 0.1
                    @test ds.attrib["output_class"] == "GLOBAL"
                    @test ds.attrib["run_Δt"] == 0.01 && ds.attrib["run_scheme"] == "RK4"
                end
                NCDataset(spath) do ds
                    @test ds.dim["time"] == 3 && ds.dim["shell"] == length(g.shell_k)
                    @test isequal(Array(ds["shell_k"][:]), g.shell_k) && eltype(ds["shell_k"]) == T
                    @test isequal(Array(ds["shell_edge"][:]), g.shell_edge)
                    @test ds["shell_count"][:] == g.shell_count && ds["shell_count_all"][:] == g.shell_count_all
                    for (name, d) in spectral_d
                        @test eltype(ds[String(name)]) == T
                        @test isequal(Array(ds[String(name)][:, :]), reduce(hcat, expected_spectral[name]))
                        @test ds[String(name)].attrib["coordinates"] == String(d.coordinate)
                    end
                    @test ds["Pi_kin"].attrib["coordinates"] == "shell_edge"
                    @test ds["E_kin_inc_k"].attrib["coordinates"] == "shell_k"
                end
            end
        end

        @testset "parameters: attributes and params.toml: $T" for T in (Float64, Float32)
            for compressible in (false, true)
                g, problem, state, ws = setup(backend, T; compressible)
                mktempdir() do dir
                    path = write_parameters(joinpath(dir, "params.toml"), problem; run = (; Δt = 0.01))
                    p = TOML.parsefile(path)
                    @test p["float_type"] == string(T)
                    @test p["grid"]["Nx"] == 16 && isapprox(p["grid"]["Lx"], 2π; rtol = 1e-6)
                    @test p["forcing"]["type"] == "OUForcing" && p["forcing"]["seed"] == 3
                    @test isapprox(p["forcing"]["τ_f"], 0.2; rtol = 1e-6) && !haskey(p["forcing"], "z")
                    @test p["model"]["type"] == (compressible ? "Compressible" : "Incompressible")
                    @test p["run"]["Δt"] == 0.01
                    compressible && @test p["model"]["eos"]["cs"] == 10 && isapprox(p["model"]["μ"], 0.01; rtol = 1e-6)
                    snap = write_snapshot(dir, 7, state, problem, ws)
                    NCDataset(snap) do ds
                        @test ds.attrib["model_large_scale_m"] == 0
                        @test ds.attrib["forcing_seed"] == 3 && ds.attrib["grid_Ny"] == 16
                        @test ds.attrib["output_class"] == "FIELDS" && !haskey(ds.attrib, "period")
                        compressible && @test ds.attrib["model_eos_cs"] == 10
                    end
                end
            end
        end

        @testset "snapshot restart error: $(compressible ? "compressible" : "incompressible"), $T" for T in (Float64, Float32),
                                                                                                       compressible in (false, true)
            g, problem, state, ws = setup(backend, T; compressible)
            advance!(state, problem, ws, 5, T(0.005))
            mktempdir() do dir
                path = write_snapshot(dir, 12, state, problem, ws)
                @test basename(path) == "snap_000012.nc"
                NCDataset(path) do ds
                    @test ds.dim["x"] == 16 && ds.dim["y"] == 16
                    @test keys(state.fields) == (compressible ? (:ux, :uy, :rho) : (:ux, :uy))
                    @test all(name -> eltype(ds[String(name)]) == T && size(ds[String(name)]) == (16, 16), keys(state.fields))
                    @test eltype(ds["t"]) == T && ds["t"][] == state.t && ds["n"][] == 5
                    @test eltype(ds["forcing_z"]) == T && size(ds["forcing_z"]) == (2, 2 * mode_count_stored(problem.forcing))
                end
                _, problem2, state2, ws2 = setup(backend, T; compressible)
                read_snapshot!(state2, problem2, ws2, path)
                @test state2.t == state.t && state2.n == state.n
                @test state_error(state2, state) <= 10 * eps(T)
                for name in keys(state.fields)
                    @test relative_error(state2.fields[name], state.fields[name]) <= 10 * eps(T)
                end
                @test isequal(CompReS.forcing_state(problem2.forcing), CompReS.forcing_state(problem.forcing))
            end
        end

        @testset "restart rejects another run" begin
            g, problem, state, ws = setup(backend, Float64)
            mktempdir() do dir
                path = write_snapshot(dir, 0, state, problem, ws)
                _, other_seed, s2, w2 = setup(backend, Float64)
                other_seed = Problem(g, problem.model, OUForcing(g, band(2, 4); ε = 1.0, τ_f = 0.2, seed = 4))
                @test_throws ArgumentError read_snapshot!(s2, other_seed, w2, path)
                _, other_grid, s3, w3 = setup(backend, Float64; N = 32)
                @test_throws ArgumentError read_snapshot!(s3, other_grid, w3, path)
                _, compressible, s4, w4 = setup(backend, Float64; compressible = true)
                @test_throws ArgumentError read_snapshot!(s4, compressible, w4, path)
            end
        end

        @testset "2n steps equal n steps, restart, n steps: $(compressible ? "compressible" : "incompressible"), $T" for T in (Float64, Float32),
                                                                                                                      compressible in (false, true)
            n, Δt = 10, T(0.005)
            _, problem, straight, ws = setup(backend, T; compressible)
            advance!(straight, problem, ws, 2n, Δt)

            _, problem_a, first_half, ws_a = setup(backend, T; compressible)
            advance!(first_half, problem_a, ws_a, n, Δt)
            mktempdir() do dir
                path = write_snapshot(dir, 1, first_half, problem_a, ws_a)
                _, problem_b, resumed, ws_b = setup(backend, T; compressible)
                read_snapshot!(resumed, problem_b, ws_b, path)
                advance!(resumed, problem_b, ws_b, n, Δt)
                @test resumed.n == straight.n == 2n
                @test resumed.t == straight.t
                @test state_error(resumed, straight) <= 100 * eps(T)
                @test CompReS.forcing_state(problem_b.forcing).z ≈ CompReS.forcing_state(problem.forcing).z rtol = 100 * eps(T)
            end
        end

        @testset "upscaling: one Fourier mode, 8 → 16: $T" for T in (Float64, Float32)
            # The solenoidal mode k = (1, 2) and the kx = 0 pair (0, ±1); |k| < k_c = 8/3 on N = 8.
            mode_fields(N) = (ux = spectral_modes(N, T, [(1, 2) => Complex{T}(0.3, 0.1), (0, 1) => Complex{T}(0.2, -0.4)]),
                              uy = spectral_modes(N, T, [(1, 2) => Complex{T}(-0.15, -0.05)]))
            g8, problem8, state8, ws8 = unforced_setup(backend, T, 8, mode_fields(8))
            g16, problem16, state16, ws16 = unforced_setup(backend, T, 16, nothing)
            mktempdir() do dir
                path = write_snapshot(dir, 0, state8, problem8, ws8)
                read_snapshot!(state16, problem16, ws16, path; upscale = true)
                for name in (:ux, :uy)
                    expected = mode_fields(16)[name]
                    @test maximum(abs, Array(state16.fields[name]) .- expected) <= 10 * eps(T) * maximum(abs, expected)
                end
            end
        end

        @testset "upscaling: round trip 16 → 32 → 16 and invariants: $(compressible ? "compressible" : "incompressible"), $T" for T in (Float64, Float32),
                                                                                                                           compressible in (false, true)
            _, problem16, state16, ws16 = setup(backend, T; compressible)
            if compressible
                # A nonuniform density with the mean 1.
                rho = Array(state16.fields.rho)
                rho[2, 3] = Complex{T}(0.01, 0.02)
                copyto!(state16.fields.rho, rho)
                enforce_hermitian!(state16.fields.rho, problem16.grid)
            end
            _, problem32, state32, ws32 = setup(backend, T; compressible, N = 32)
            _, problem_back, state_back, ws_back = setup(backend, T; compressible)
            mktempdir() do dir
                path16 = write_snapshot(dir, 4, state16, problem16, ws16)
                @test_logs min_level = Base.CoreLogging.Warn read_snapshot!(state32, problem32, ws32, path16; upscale = true)
                @test CompReS.kinetic_energy(state32, problem32, ws32) ≈ CompReS.kinetic_energy(state16, problem16, ws16) rtol = 10 * eps(T)
                compressible && @test Array(state32.fields.rho)[1, 1] ≈ 1 rtol = 10 * eps(T)
                path32 = write_snapshot(mkpath(joinpath(dir, "fine")), 5, state32, problem32, ws32)
                read_snapshot!(state_back, problem_back, ws_back, path32; upscale = true)
                @test state_error(state_back, state16) <= 10 * eps(T)
                @test state_back.t == state16.t && state_back.n == state16.n
            end
        end

        @testset "upscaling: coarsening 32 → 16 truncates to the new disc" begin
            T = Float64
            # |(1, 1)| < k_c = 16/3 on N = 16; |(6, 0)| lies outside it and inside k_c = 32/3 on N = 32.
            fine = (ux = spectral_modes(32, T, [(1, 1) => 0.5 + 0.2im, (6, 0) => 0.1 + 0.0im]),
                    uy = spectral_modes(32, T, [(1, 1) => -0.5 - 0.2im]))
            _, problem32, state32, ws32 = unforced_setup(backend, T, 32, fine)
            _, problem16, state16, ws16 = unforced_setup(backend, T, 16, nothing)
            mktempdir() do dir
                path = write_snapshot(dir, 0, state32, problem32, ws32)
                read_snapshot!(state16, problem16, ws16, path; upscale = true)
                @test maximum(abs, Array(state16.fields.ux) .- spectral_modes(16, T, [(1, 1) => 0.5 + 0.2im])) <= 10 * eps(T)
                @test maximum(abs, Array(state16.fields.uy) .- spectral_modes(16, T, [(1, 1) => -0.5 - 0.2im])) <= 10 * eps(T)
            end
        end

        @testset "upscaling: guards and forcing mismatch" begin
            T = Float64
            g, problem, state, ws = setup(backend, T)
            mktempdir() do dir
                path = write_snapshot(dir, 0, state, problem, ws)
                saved_z = CompReS.forcing_state(problem.forcing).z
                # Without upscale, another grid size still throws.
                _, p32, s32, w32 = setup(backend, T; N = 32)
                @test_throws ArgumentError read_snapshot!(s32, p32, w32, path)
                # Another box throws even with upscale.
                big = Grid(; Lx = 4π, Ly = 4π, Nx = 32, Ny = 32, T, backend)
                sb = State(random_velocity(big; spectrum = SPECTRUM, E_kin = 0.5, seed = 2))
                pb = Problem(big, Incompressible(; ν = 0.01), NoForcing())
                @test_throws ArgumentError read_snapshot!(sb, pb, Workspace(big, sb), path; upscale = true)
                # Another seed: a warning, and the forcing state is restored.
                g32 = p32.grid
                p_seed = Problem(g32, p32.model, OUForcing(g32, band(2, 4); ε = 1.0, τ_f = 0.2, seed = 4))
                @test_logs (:warn, r"seed") read_snapshot!(s32, p_seed, w32, path; upscale = true)
                @test CompReS.forcing_state(p_seed.forcing).z == saved_z
                # Another forcing type: a warning, and the forcing keeps its initial state.
                p_type = Problem(g32, p32.model, RandomPhaseForcing(g32, band(2, 4); f_rms = 1.0, τ_r = 0.1, seed = 3))
                initial = CompReS.forcing_state(p_type.forcing)
                @test_logs (:warn, r"forcing type") read_snapshot!(s32, p_type, w32, path; upscale = true)
                @test CompReS.forcing_state(p_type.forcing) == initial
                # Another forced mode count: a warning, and the forcing keeps its initial state.
                p_band = Problem(g32, p32.model, OUForcing(g32, band(2, 6); ε = 1.0, τ_f = 0.2, seed = 3))
                initial = CompReS.forcing_state(p_band.forcing)
                @test_logs (:warn, r"forced modes") read_snapshot!(s32, p_band, w32, path; upscale = true)
                @test CompReS.forcing_state(p_band.forcing) == initial
            end
        end
    end
end

# Host spectral field of an N × N grid with the given modes (i, j) => value, i ≥ 0, signed j. A kx = 0 mode
# also sets its conjugate (0, −j).
function spectral_modes(N, T, values)
    f̂ = zeros(Complex{T}, N ÷ 2 + 1, N)
    row(j) = mod(j, N) + 1
    for ((i, j), v) in values
        f̂[i + 1, row(j)] = v
        i == 0 && (f̂[1, row(-j)] = conj(v))
    end
    return f̂
end

# Grid of side 2π, unforced incompressible problem, state and workspace; `fields` holds host spectral ux and
# uy, or is nothing for a state at rest.
function unforced_setup(backend, T, N, fields)
    g = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N, T, backend)
    host = isnothing(fields) ? (ux = spectral_modes(N, T, []), uy = spectral_modes(N, T, [])) : fields
    state = State(map(f -> CompReS.to_device(backend, f), host))
    return g, Problem(g, Incompressible(; ν = 0.01), NoForcing()), state, Workspace(g, state)
end

# Number of stored forced modes of a mode forcing.
mode_count_stored(f) = length(f.z) ÷ 2

end
