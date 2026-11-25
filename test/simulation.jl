# Tests of the simulation driver (WP-13, ADR-0009, ADR-0010, ADR-0027): output files and sampling times, the
# walltime stop, the NaN, density and CFL checks, restart, parameter conversion and the target user script.
module SimulationTests

using Test
using KernelAbstractions
using NCDatasets: NCDataset
using TOML
using CompReS

const SPECTRUM = k -> k^4 * exp(-k^2 / 8)

model_of(compressible) =
    compressible ? Compressible(; eos = Isothermal(cs = 10), μ = 0.01, ζ = 0.005) : Incompressible(; ν = 0.01)

# A forced simulation on an N × N grid of side 2π, with a random initial velocity.
function forced_simulation(backend, T, outdir; compressible = false, N = 16, Δt = 0.005, kwargs...)
    g = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N, T, backend)
    forcing = OUForcing(g, band(2, 4); ε = 1.0, τ_f = 0.2, seed = 3)
    sim = Simulation(g, model_of(compressible); forcing, Δt, outdir, kwargs...)
    set_initial_condition!(sim, random_velocity(g; spectrum = SPECTRUM, E_kin = 0.5, seed = 2))
    return sim
end

# Relative 2-norm distance of two states over all their fields.
state_error(a, b) = sqrt(sum(name -> sum(abs2, Array(a.fields[name]) .- Array(b.fields[name])), keys(a.fields)) /
                         sum(name -> sum(abs2, Array(b.fields[name])), keys(b.fields)))

snapshots(outdir) = sort!(filter(startswith("snap_"), readdir(joinpath(outdir, "fields"))))

# The target script of design/origin/ROADMAP-build.md at a tiny size: N = 32, Δt = 1e-3, t_end = 0.05, then a restart to 0.1.
function target_script(backend, outdir; compressible)
    grid    = Grid(; Lx = 2π, Ly = 2π, Nx = 32, Ny = 32, T = Float64, backend)
    model   = compressible ?
              Compressible(; eos = Isothermal(cs = 10.0), μ = 1e-4, ζ = 0.0, κ = 0.0,
                           large_scale = LargeScaleDissipation(α = 0.1, m = 0)) :
              Incompressible(; ν = 1e-4, large_scale = LargeScaleDissipation(α = 0.1, m = 0))
    forcing = OUForcing(grid, band(4, 6); ε = 1.0, τ_f = 0.05, seed = 1)
    new_simulation() =
        Simulation(grid, model; forcing, scheme = RK4(), integrating_factor = true, Δt = 1e-3,
                   outputs = (GLOBAL   = Output(period = 0.005),
                              SPECTRAL = Output(period = 0.01),
                              FIELDS   = Output(period = 0.05)),
                   outdir, max_walltime = 23.5 * 3600)
    sim = new_simulation()
    set_initial_condition!(sim, random_velocity(grid; spectrum = k -> k^4 * exp(-k^2 / 8), E_kin = 0.5, seed = 2))
    run!(sim; t_end = 0.05)
    sim = new_simulation()
    restart!(sim, joinpath(outdir, "fields", "snap_000001.nc"))
    run!(sim; t_end = 0.1)
    return sim
end

function run_tests(backend)
    @testset "simulation" begin
        @testset "a tiny forced run produces all files: $T" for T in (Float64, Float32)
            mktempdir() do dir
                Δt = 0.003
                outputs = (GLOBAL = Output(period = 0.01), SPECTRAL = Output(period = 0.05),
                           FIELDS = Output(period = 0.15))
                sim = forced_simulation(backend, T, dir; Δt, outputs)
                @test run!(sim; t_end = 0.3) === sim
                # The run ends at the step nearest to t_end; times follow the step count, without drift.
                @test sim.state.n == 100
                @test sim.state.t === T(100 * Float64(T(Δt)))
                @test sort!(readdir(dir)) == ["fields", "global.nc", "params.toml", "run.log", "spectral.nc"]
                @test snapshots(dir) == ["snap_000000.nc", "snap_000001.nc", "snap_000002.nc"]
                # Each class samples at t₀ and at the step nearest to each multiple of its period.
                expected_steps(period) = [round(Int, j * period / Δt) for j in 0:round(Int, 0.3 / period)]
                NCDataset(joinpath(dir, "global.nc")) do ds
                    @test ds["step"][:] == expected_steps(0.01)
                    @test ds["time"][:] == T.(expected_steps(0.01) .* Float64(T(Δt)))
                    @test eltype(ds["time"]) == T && eltype(ds["E_kin"]) == T
                    @test isnan(ds["inj_E_kin"][1]) && all(!isnan, ds["inj_E_kin"][2:end])
                    @test ds.attrib["run_Δt"] === T(Δt) && ds.attrib["model_ν"] isa T && ds.attrib["forcing_ε"] isa T
                    @test ds.attrib["run_scheme"] == "RK4" && ds.attrib["run_C_stop"] === one(T)
                end
                NCDataset(joinpath(dir, "spectral.nc")) do ds
                    @test ds["step"][:] == expected_steps(0.05)
                    @test eltype(ds["Pi_kin"]) == T
                end
                for (i, name) in enumerate(snapshots(dir))
                    NCDataset(joinpath(dir, "fields", name)) do ds
                        @test ds["n"][] == expected_steps(0.15)[i] && ds.attrib["snapshot_index"] == i - 1
                        @test eltype(ds["ux"]) == T && eltype(ds["t"]) == T
                    end
                end
                params = TOML.parsefile(joinpath(dir, "params.toml"))
                @test params["run"]["n_check"] == 10 && params["run"]["integrating_factor"] == true
                @test params["run"]["FIELDS_period"] == 0.15
                log = read(joinpath(dir, "run.log"), String)
                @test occursin("stop: t_end = 0.3 reached", log)
                @test occursin("step", log) && occursin("GLOBAL", log) && occursin("Pi_kin", log)
                @test !occursin("rhs!", log)
                @test count("t = ", log) >= 3
                # The state and the parameters keep the float type of the grid.
                @test sim.state.t isa T && eltype(sim.state.fields.ux) == Complex{T}
                @test sim.problem.model.ν isa T && sim.problem.forcing.ε isa T && sim.Δt isa T
            end
        end

        @testset "parameter conversion to the float type of the grid" begin
            g64 = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float64, backend)
            g32 = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float32, backend)
            forcing = OUForcing(g64, band(2, 4); ε = 1.0, τ_f = 0.2, seed = 3)
            mktempdir() do dir
                sim = Simulation(g32, Compressible(; eos = Isothermal(cs = 10.0), μ = 1e-2, ζ = 0.0); forcing,
                                 Δt = 1e-3, outdir = dir, C_stop = 0.5)
                @test sim.problem.model.eos.cs isa Float32 && sim.problem.model.μ isa Float32
                @test sim.problem.forcing.ε isa Float32 && eltype(sim.problem.forcing.z) == ComplexF32
                @test sim.C_stop === 0.5f0 && sim.Δt === 1f-3
                # Without set_initial_condition!, the compressible run starts from rest with ρ = 1.
                @test all(iszero, Array(sim.state.fields.ux)) && Array(sim.state.fields.rho)[1, 1] == 1
            end
        end

        @testset "walltime 0 writes one snapshot and stops" begin
            mktempdir() do dir
                outputs = (GLOBAL = Output(period = 0.01), FIELDS = Output(period = 1))
                sim = forced_simulation(backend, Float64, dir; outputs, max_walltime = 0)
                run!(sim; t_end = 1.0)
                @test sim.state.n == 0
                @test snapshots(dir) == ["snap_000000.nc"]
                NCDataset(ds -> (@test ds.dim["time"] == 1), joinpath(dir, "global.nc"))
                @test occursin("max_walltime", read(joinpath(dir, "run.log"), String))
            end
            mktempdir() do dir
                # Without a FIELDS output, the walltime stop still writes a snapshot.
                sim = forced_simulation(backend, Float64, dir; max_walltime = 0)
                run!(sim; t_end = 1.0)
                @test snapshots(dir) == ["snap_000000.nc"]
            end
        end

        @testset "an injected NaN writes a snapshot and raises an error" begin
            mktempdir() do dir
                sim = forced_simulation(backend, Float64, dir; n_check = 4, outputs = (FIELDS = Output(period = 1),))
                run!(sim; t_end = 0.025)                  # 5 steps
                ux = Array(sim.state.fields.ux)
                ux[2, 3] = NaN
                copyto!(sim.state.fields.ux, ux)
                @test_throws r"state holds NaN" run!(sim; t_end = 1.0)
                # The second run! checks the state at its start, step 5.
                @test sim.state.n == 5
                @test snapshots(dir) == ["snap_000000.nc", "snap_000001.nc"]
                NCDataset(joinpath(dir, "fields", "snap_000001.nc")) do ds
                    @test ds["n"][] == 5 && any(isnan, ds["ux"][:, :])
                end
                @test occursin("stop: error", read(joinpath(dir, "run.log"), String))
            end
        end

        @testset "the NaN check runs every n_check steps" begin
            mktempdir() do dir
                # Δt so large that the run overflows; the check every 2 steps stops it at an even step.
                g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float64, backend)
                sim = Simulation(g, Incompressible(; ν = 0); scheme = RK4(), Δt = 0.5, outdir = dir, n_check = 2,
                                 C_stop = Inf)
                set_initial_condition!(sim, random_velocity(g; spectrum = SPECTRUM, E_kin = 1e3, seed = 2))
                @test_throws r"state holds NaN" run!(sim; t_end = 1000.0)
                @test iseven(sim.state.n) && sim.state.n > 0
                @test !all(f -> all(isfinite, Array(f)), sim.state.fields)
            end
        end

        @testset "a non-positive density stops the run" begin
            mktempdir() do dir
                g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float64, backend)
                sim = Simulation(g, model_of(true); Δt = 1e-3, outdir = dir,
                                 outputs = (GLOBAL = Output(period = 0.01),))
                ρ = zeros(ComplexF64, 9, 16)
                ρ[1, 1] = 1
                ρ[2, 1] = 0.75                             # ρ = 1 + 1.5 cos(x), negative near x = π
                rho = copyto!(KernelAbstractions.allocate(backend, ComplexF64, 9, 16), ρ)
                set_initial_condition!(sim, merge(kolmogorov(g; A = 1, n = 1), (; rho)))
                @test_throws r"density" run!(sim; t_end = 0.1)
                @test snapshots(dir) == ["snap_000000.nc"]
                # The check runs before the sampling: GLOBAL holds no record of the failed state.
                NCDataset(ds -> (@test ds.dim["time"] == 0), joinpath(dir, "global.nc"))
            end
        end

        @testset "CFL number and CFL stop: $(compressible ? "compressible" : "incompressible")" for compressible in (false, true)
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float64, backend)
            model = compressible ? Compressible(; eos = Isothermal(cs = 2), μ = 0, ζ = 0) : Incompressible(; ν = 0)
            # u = sin(y) x̂: max|u_x| = 1 on the grid point y = π/2, u_y = 0, Δx = Δy = 2π/16.
            Δx = 2π / 16
            C(Δt) = compressible ? Δt * ((1 + 2) / Δx + 2 / Δx) : Δt / Δx
            mktempdir() do dir
                sim = Simulation(g, model; Δt = 0.1, outdir = dir)
                set_initial_condition!(sim, kolmogorov(g; A = 1, n = 1))
                @test CompReS.cfl_number(sim) ≈ C(0.1) rtol = 1e-12
            end
            Δt = 1.2 / C(1)                                  # C = 1.2
            mktempdir() do dir
                sim = Simulation(g, model; Δt, outdir = dir)
                set_initial_condition!(sim, kolmogorov(g; A = 1, n = 1))
                @test_throws r"CFL" run!(sim; t_end = 10Δt)
                @test sim.state.n == 0 && snapshots(dir) == ["snap_000000.nc"]
            end
            mktempdir() do dir
                sim = Simulation(g, model; Δt, outdir = dir, C_stop = 1.5)
                set_initial_condition!(sim, kolmogorov(g; A = 1, n = 1))
                run!(sim; t_end = 2Δt)
                @test sim.state.n == 2 && isapprox(sim.cfl, 1.2; rtol = 1e-6)
            end
        end

        @testset "restart! continues as an uninterrupted run: $(compressible ? "compressible" : "incompressible"), $T" for T in (Float64, Float32),
                                                                                                           compressible in (false, true)
            n, Δt = 10, 0.005
            outputs = (GLOBAL = Output(period = 2Δt), FIELDS = Output(period = n * Δt))
            mktempdir() do dir
                straight = forced_simulation(backend, T, joinpath(dir, "straight"); compressible, Δt, outputs)
                run!(straight; t_end = 2n * Δt)
                split_dir = joinpath(dir, "split")
                first_half = forced_simulation(backend, T, split_dir; compressible, Δt, outputs)
                run!(first_half; t_end = 2n * Δt)
                # Restart from the snapshot at step n; the files of split_dir hold the records up to step 2n.
                resumed = Simulation(first_half.grid, model_of(compressible);
                                     forcing = OUForcing(first_half.grid, band(2, 4); ε = 1.0, τ_f = 0.2, seed = 3),
                                     Δt, outdir = split_dir, outputs)
                restart!(resumed, joinpath(split_dir, "fields", "snap_000001.nc"))
                @test resumed.state.n == n
                run!(resumed; t_end = 2n * Δt)
                @test resumed.state.n == straight.state.n == 2n
                @test resumed.state.t ≈ straight.state.t
                @test state_error(resumed.state, straight.state) <= 100 * eps(T)
                # The restart keeps the records before step n, samples at step n again and continues.
                NCDataset(joinpath(split_dir, "global.nc")) do ds
                    @test ds["step"][:] == 0:2:2n
                    inj = ds["inj_E_kin"][:]
                    @test isnan(inj[1]) && isnan(inj[n ÷ 2 + 1]) && count(isnan, inj) == 2
                end
                @test snapshots(split_dir) == ["snap_000000.nc", "snap_000001.nc", "snap_000002.nc", "snap_000003.nc"]
                NCDataset(ds -> (@test ds["n"][] == 2n), joinpath(split_dir, "fields", "snap_000003.nc"))
            end
        end

        @testset "upscaling restart: $(compressible ? "compressible" : "incompressible")" for compressible in (false, true)
            T, n, Δt = Float64, 4, 0.005
            outputs = (GLOBAL = Output(period = Δt), SPECTRAL = Output(period = 2Δt), FIELDS = Output(period = n * Δt))
            mktempdir() do dir
                coarse_dir, fine_dir = joinpath(dir, "coarse"), joinpath(dir, "fine")
                coarse = forced_simulation(backend, T, coarse_dir; compressible, Δt, outputs)
                run!(coarse; t_end = n * Δt)
                source = joinpath(coarse_dir, "fields", "snap_000001.nc")
                fine_simulation(outdir) = forced_simulation(backend, T, outdir; compressible, N = 32, Δt, outputs)
                # The source outdir holds output of another grid size.
                @test_throws ArgumentError restart!(fine_simulation(coarse_dir), source; upscale = true)
                # Without upscale, another grid size throws.
                @test_throws ArgumentError restart!(fine_simulation(fine_dir), source)
                fine = restart!(fine_simulation(fine_dir), source; upscale = true)
                @test fine.state.n == n && fine.state.t ≈ n * Δt
                run!(fine; t_end = 2n * Δt)
                @test fine.state.n == 2n
                @test all(f -> all(isfinite, Array(f)), fine.state.fields)
                @test snapshots(fine_dir) == ["snap_000002.nc", "snap_000003.nc"]
                for file in ("global.nc", "spectral.nc", joinpath("fields", "snap_000003.nc"))
                    NCDataset(joinpath(fine_dir, file)) do ds
                        @test ds.attrib["run_upscaled_from"] == source
                        @test ds.attrib["run_upscaled_from_Nx"] == 16 && ds.attrib["run_upscaled_from_Ny"] == 16
                        @test ds.attrib["grid_Nx"] == 32
                    end
                end
                @test TOML.parsefile(joinpath(fine_dir, "params.toml"))["run"]["upscaled_from"] == source
                # A plain restart of the fine run keeps the provenance.
                resumed = restart!(fine_simulation(fine_dir), joinpath(fine_dir, "fields", "snap_000003.nc"))
                run!(resumed; t_end = 3n * Δt)
                params = TOML.parsefile(joinpath(fine_dir, "params.toml"))["run"]
                @test params["upscaled_from"] == source && params["upscaled_from_Nx"] == 16
                # A new run clears it.
                fresh = fine_simulation(joinpath(dir, "fresh"))
                run!(fresh; t_end = Δt)
                @test !haskey(TOML.parsefile(joinpath(dir, "fresh", "params.toml"))["run"], "upscaled_from")
            end
        end

        @testset "target user script, tiny: $(compressible ? "compressible" : "incompressible")" for compressible in (false, true)
            mktempdir() do dir
                sim = target_script(backend, dir; compressible)
                @test sim.state.n == 100 && sim.state.t ≈ 0.1
                @test snapshots(dir) == ["snap_000000.nc", "snap_000001.nc", "snap_000002.nc", "snap_000003.nc"]
                NCDataset(ds -> (@test ds["step"][:] == 0:5:100), joinpath(dir, "global.nc"))
                NCDataset(ds -> (@test ds["step"][:] == 0:10:100), joinpath(dir, "spectral.nc"))
                log = read(joinpath(dir, "run.log"), String)
                @test count("stop: t_end", log) == 2
            end
        end

        @testset "a second run! continues the run" begin
            mktempdir() do dir
                sim = forced_simulation(backend, Float64, dir; outputs = (GLOBAL = Output(period = 0.01),))
                run!(sim; t_end = 0.02)
                run!(sim; t_end = 0.04)
                NCDataset(joinpath(dir, "global.nc")) do ds
                    @test ds["step"][:] == 0:2:8
                    # Each segment samples at its t₀ with NaN injection rates.
                    @test isequal(isnan.(ds["inj_E_kin"][:]), [true, false, true, false, false])
                end
                @test count("== run segment", read(joinpath(dir, "run.log"), String)) == 2
            end
        end

        @testset "profile mode times the tasks of step!" begin
            mktempdir() do dir
                sim = forced_simulation(backend, Float64, dir; profile = true, scheme = RK4(), integrating_factor = false)
                run!(sim; t_end = 0.02)
                log = read(joinpath(dir, "run.log"), String)
                @test all(label -> occursin(label, log),
                          ("update_forcing!", "record_injection!", "rhs!", "linear terms", "stage combination"))
            end
        end

        @testset "scheme checks" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float32, backend)
            forcing = OUForcing(g, band(2, 4); ε = 1.0, τ_f = 0.2, seed = 3)
            mktempdir() do dir
                @test_logs (:warn, r"Euler") Simulation(g, Incompressible(; ν = 0.01); forcing, scheme = Euler(),
                                                        Δt = 0.01, outdir = dir)
                @test_logs (:warn, r"RK2") Simulation(g, Incompressible(; ν = 0.01); forcing, scheme = RK2(),
                                                      Δt = 0.01, outdir = dir)
                @test_logs Simulation(g, Incompressible(; ν = 0.01); scheme = Euler(), Δt = 0.01, outdir = dir)
                # SSPRK3 with the integrating factor applies exp(σΔt/2). The rates cover the grid corner, k² = 128:
                # σΔt/2 = 10·128·Δt/2 is 640 for Δt = 1 and 64 for Δt = 0.1; log(floatmax(Float32)) ≈ 88.7.
                @test_throws r"SSPRK3" Simulation(g, Incompressible(; ν = 10); scheme = SSPRK3(), Δt = 1, outdir = dir)
                @test Simulation(g, Incompressible(; ν = 10); scheme = SSPRK3(), Δt = 0.1, outdir = dir) isa Simulation
                @test Simulation(g, Incompressible(; ν = 10); scheme = RK4(), Δt = 1, outdir = dir) isa Simulation
                @test Simulation(g, Incompressible(; ν = 10); scheme = SSPRK3(), integrating_factor = false, Δt = 1,
                                 outdir = dir) isa Simulation
            end
        end

        @testset "argument checks" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16, T = Float64, backend)
            model = Incompressible(; ν = 0.01)
            mktempdir() do dir
                @test_throws ArgumentError Simulation(g, model; Δt = 0, outdir = dir)
                @test_throws ArgumentError Simulation(g, model; Δt = 0.01, outdir = dir, n_check = 0)
                @test_throws ArgumentError Simulation(g, model; Δt = 0.01, outdir = dir, C_stop = 0)
                @test_throws ArgumentError Simulation(g, model; Δt = 0.01, outdir = dir, outputs = (FOO = Output(period = 1),))
                @test_throws ArgumentError Simulation(g, model; Δt = 0.01, outdir = dir,
                                                      outputs = (GLOBAL = Output(period = 1, diagnostics = [:nope]),))
                sim = Simulation(g, model; Δt = 0.01, outdir = dir)
                @test_throws ArgumentError set_initial_condition!(sim, (; ux = sim.state.fields.ux))
                @test_throws ArgumentError set_initial_condition!(sim, merge(kolmogorov(g; A = 1, n = 1), uniform_density(g)))
                @test_throws ArgumentError run!(sim; t_end = -1.0)
            end
            mktempdir() do dir
                # A new run does not overwrite the output of another run, unless overwrite = true.
                sim = forced_simulation(backend, Float64, dir; outputs = (GLOBAL = Output(period = 0.01),))
                run!(sim; t_end = 0.02)
                again = forced_simulation(backend, Float64, dir; outputs = (GLOBAL = Output(period = 0.01),))
                @test_throws ArgumentError run!(again; t_end = 0.02)
                again = forced_simulation(backend, Float64, dir; overwrite = true, outputs = (GLOBAL = Output(period = 0.01),))
                run!(again; t_end = 0.01)
                NCDataset(ds -> (@test ds["step"][:] == [0, 2]), joinpath(dir, "global.nc"))
                @test count("stop: t_end", read(joinpath(dir, "run.log"), String)) == 1
            end
        end
    end
end

end
