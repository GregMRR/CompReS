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
using Random: Xoshiro

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

end
