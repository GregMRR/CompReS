# Entry point of the NS2D package: each source file is included in dependency order, followed by its exports.
"""
    NS2D

Pseudospectral solver for forced-dissipated 2D turbulence, compressible (low Mach) and
incompressible, on CPU threads, NVIDIA GPUs and AMD GPUs.
"""
module NS2D

using KernelAbstractions
using AbstractFFTs: plan_rfft, plan_brfft
using FFTW: FFTW
using LinearAlgebra: mul!

export CPU

include("grid.jl")
export Grid

include("transforms.jl")
export Transforms, to_spectral!, to_physical!, truncate!, ∂x!, ∂y!, helmholtz_split!, project_solenoidal!,
    enforce_hermitian!

include("eos.jl")
export Isothermal, Polytropic, convert_float, pressure, enthalpy, internal_energy, sound_speed²,
    density_from_enthalpy

end
