# Entry point of the CompReS package: each source file is included in dependency order, followed by its exports.
"""
    CompReS

Pseudospectral solver for forced-dissipated 2D turbulence, compressible (low Mach) and
incompressible, on CPU threads, NVIDIA GPUs and AMD GPUs.
"""
module CompReS

using KernelAbstractions
using AbstractFFTs: plan_rfft, plan_brfft
using FFTW: FFTW
using LinearAlgebra: mul!
using Random: Xoshiro, randn!
using NCDatasets: NCDataset, defDim, defVar, dimnames, sync
using TimerOutputs: TimerOutput, @timeit, print_timer, reset_timer!
using TOML: TOML

export CPU

include("grid.jl")
export Grid

include("transforms.jl")
export Transforms, to_spectral!, to_physical!, truncate!, ∂x!, ∂y!, helmholtz_split!, project_solenoidal!,
    enforce_hermitian!

include("eos.jl")
export Isothermal, Polytropic, convert_float, pressure, enthalpy, internal_energy, sound_speed²,
    density_from_enthalpy

include("state.jl")
export State, Workspace

include("initial_conditions.jl")
export random_velocity, kolmogorov, crossed_kolmogorov, uniform_density

include("linear.jl")
export LargeScaleDissipation, Incompressible, Compressible, linear_rates, linear_tendency!, propagate!

include("nonlinear_incompressible.jl")
export Problem, NoForcing, rhs!

include("forcing.jl")
export band, modes, support, KolmogorovForcing, OUForcing, RandomPhaseForcing, mode_count

include("nonlinear_compressible.jl")

include("diagnostics/registry.jl")
export Diagnostic, register_diagnostic!, select_diagnostics, default_diagnostics, record_injection!

include("diagnostics/incompressible.jl")
include("diagnostics/compressible.jl")

include("timestepping.jl")
export RungeKuttaScheme, Euler, RK2, SSPRK3, RK4, step!

include("balanced_density.jl")
export balanced_density

include("io.jl")
export Output, selected_diagnostics, parameters, write_parameters, open_global, write_global!, open_spectral,
    write_spectral!, write_snapshot, read_snapshot!

include("simulation.jl")
export Simulation, set_initial_condition!, run!, restart!

end
