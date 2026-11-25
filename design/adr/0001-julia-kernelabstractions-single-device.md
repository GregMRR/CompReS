# Julia, KernelAbstractions.jl, one device per run

Date: 2026-11-19. Status: accepted.

The solver is a Julia package. Every custom kernel uses KernelAbstractions.jl, so one source runs on CPU threads, NVIDIA (CUDA.jl) and AMD (AMDGPU.jl). FFTs go through the AbstractFFTs.jl interface: FFTW on CPU, cuFFT and rocFFT on GPU. CUDA.jl and AMDGPU.jl are weak dependencies (package extensions), so a CPU run loads no GPU stack.

A run uses one device: one CPU node with threads, or one GPU. A single 80 GB GPU holds a 2D grid of about 16384², which covers the planned production runs. There is no MPI decomposition.

## Consequences

- Keep the state in plain device arrays and keep all FFT calls behind the transform module, so a later pencil-decomposition backend replaces one module only.
- CPU and CUDA tests decide pass/reject on the development machine. AMD tests are written and stored but not run (ADR-0011).
