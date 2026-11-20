# Tests of the plotting extensions (WP-15, ADR-0012, ADR-0029): each plotting function returns a figure in
# CompReSMakieExt (through CairoMakie) and CompReSPlotsExt (through Plots with the headless GR backend), and
# animate_snapshots writes a file. Plotting runs on the CPU, so the tests run on the CPU backend only.
module PlottingTests

ENV["GKSwstype"] = "100"   # GR without a display

using Test
using CompReS
import CairoMakie
import Plots

# A tiny incompressible and a tiny compressible crossed Kolmogorov run with all three output classes.
function tiny_run(dir; compressible)
    g = Grid(; Lx = 2π, Ly = 2π, Nx = 16, Ny = 16)
    model = compressible ? Compressible(; eos = Isothermal(cs = 5), μ = 0.05, ζ = 0) : Incompressible(; ν = 0.05)
    sim = Simulation(g, model; Δt = 0.01, outdir = dir,
                     outputs = (GLOBAL = Output(period = 0.02), SPECTRAL = Output(period = 0.05),
                                FIELDS = Output(period = 0.05)))
    ic = crossed_kolmogorov(g; A1 = 1, k1 = 1, A2 = 1, k2 = 2)
    set_initial_condition!(sim, compressible ? merge(ic, uniform_density(g)) : ic)
    run!(sim; t_end = 0.15)
    return dir
end

function run_tests(backend)
    Main.is_cpu(backend) || return
    @testset "plotting" begin
        @testset "vorticity of a snapshot" begin
            mktempdir() do dir
                tiny_run(dir; compressible = false)
                # The crossed Kolmogorov flow has ω = -A₁k₁ cos(k₁y) + A₂k₂ cos(k₂x); after the start time the
                # nonlinear term changes it, so the check uses snapshot 0.
                (; x, y, values) = CompReS.snapshot_field(joinpath(dir, "fields", "snap_000000.nc"), :omega)
                @test values ≈ [-cos(yj) + 2cos(2xi) for xi in x, yj in y] atol = 1e-12
            end
        end
        @testset "both extensions loaded: library is required" begin
            mktempdir() do dir
                tiny_run(dir; compressible = false)
                @test_throws ArgumentError plot_field(dir)
            end
        end
        @testset "$library, compressible = $compressible" for library in (:Makie, :Plots), compressible in (false, true)
            figure_type = library === :Makie ? CairoMakie.Figure : Plots.Plot
            mktempdir() do dir
                tiny_run(dir; compressible)
                @test plot_field(dir; library) isa figure_type
                @test plot_field(joinpath(dir, "fields", "snap_000001.nc"); field = :ux, library) isa figure_type
                compressible && (@test plot_field(dir; field = :rho, library) isa figure_type)
                @test plot_spectrum(dir; library) isa figure_type
                @test plot_spectrum(joinpath(dir, "spectral.nc"); names = (:E_kin_inc_k,), time = 0.05, library) isa
                      figure_type
                @test plot_global(dir; library) isa figure_type
                @test plot_global(joinpath(dir, "global.nc"); names = (:E_kin, :Z), library) isa figure_type
                @test_throws ArgumentError plot_spectrum(dir; names = (:Pi_kin,), library)
                @test_throws ArgumentError plot_field(dir; field = :pressure, library)
                for name in ("anim.mp4", "anim.gif")
                    output = joinpath(dir, name)
                    @test animate_snapshots(dir, output; framerate = 5, library) == output
                    @test filesize(output) > 0
                end
            end
        end
    end
end

end
