# User acceptance run (ADR-0011, ADR-0029): the unforced, incompressible crossed Kolmogorov flow
# u = A₁ sin(k₁y) x̂ + A₂ sin(k₂x) ŷ, animated as the vorticity and the tail of E_kin_inc(k).
#
# Run: julia --project=test test/acceptance/crossed_kolmogorov.jl [outdir]
# The run uses CUDA when it is functional, else the CPU threads. The output files go to `outdir`
# (default runs/acceptance_crossed_kolmogorov, ignored by git); the animation goes to
# test/acceptance/crossed_kolmogorov.mp4.
#
# ν is fixed once by the rule of ADR-0029, which amends ADR-0011: the smallest value of a scan whose vorticity
# part with |k| > 0.9 k_max stays below 1% of max|ω|, with Re = A₁/(k₁ν) ≥ 100. Here k_max/k_ν = 4.3 at the
# palinstrophy peak, k_ν = (2νP/ν³)^{1/6}, and Re = 1250. test/acceptance/RESULTS.md records the measured values.
using CompReS
using CUDA
import CairoMakie

const ν = 4e-4
const N = 512
const T_END = 30.0

backend = CUDA.functional() ? CUDABackend() : CPU()
outdir = isempty(ARGS) ? joinpath(@__DIR__, "..", "..", "runs", "acceptance_crossed_kolmogorov") : ARGS[1]
mkpath(outdir)

grid = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N, T = Float64, backend)
sim = Simulation(grid, Incompressible(; ν); scheme = RK4(), integrating_factor = true, Δt = 1e-3,
                 outputs = (GLOBAL = Output(period = 0.05), SPECTRAL = Output(period = 0.25),
                            FIELDS = Output(period = 0.25)),
                 outdir, overwrite = true)
set_initial_condition!(sim, crossed_kolmogorov(grid; A1 = 1, k1 = 2, A2 = 1, k2 = 3))
run!(sim; t_end = T_END)

animate_snapshots(outdir, joinpath(@__DIR__, "crossed_kolmogorov.mp4"); field = :omega, spectrum = :E_kin_inc_k,
                  framerate = 10, decades = 16)
