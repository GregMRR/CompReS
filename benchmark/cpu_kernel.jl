# Time of one CompReS kernel under the KernelAbstractions CPU() backend against a plain Threads.@threads loop
# (the POCL backend check of ADR-0030, WP-16).
#
#   julia --project=benchmark -t auto benchmark/cpu_kernel.jl [N]
#
# The kernel is `propagate!` of the incompressible model: `out = exp(-στ)·û` on both velocity components,
# a pointwise kernel that the integrating-factor stages call 10 times per RK4 step.
using CompReS, Random, Statistics

N = isempty(ARGS) ? 2048 : parse(Int, ARGS[1])
grid = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N)
rates = linear_rates(grid, Incompressible(; ν = 0.01)).velocity
û = (ux = randn(Complex{Float64}, size(rates)), uy = randn(Complex{Float64}, size(rates)))
out = (ux = similar(û.ux), uy = similar(û.uy))
τ = 1e-4

function threaded_propagate!(out, û, σ, τ)
    Threads.@threads :static for I in eachindex(σ)
        @inbounds begin
            e = exp(-σ[I] * τ)
            out.ux[I] = e * û.ux[I]
            out.uy[I] = e * û.uy[I]
        end
    end
    return out
end

median_time(f; reps = 200) = median([(@elapsed f()) for _ in 1:reps])

ka_call() = propagate!(out, û, CompReS.IncompressibleRates(rates), τ)
threaded_call() = threaded_propagate!(out, û, rates, τ)
ka_call(); threaded_call()
reference = deepcopy(out)
ka_call()
@assert out.ux ≈ reference.ux && out.uy ≈ reference.uy "the two versions disagree"
t_ka, t_threads = median_time(ka_call), median_time(threaded_call)
println("N = $N, $(Threads.nthreads()) threads, propagate! on $(length(rates)) modes")
println("KernelAbstractions CPU(): ", round(1e3 * t_ka; sigdigits = 4), " ms")
println("Threads.@threads loop:    ", round(1e3 * t_threads; sigdigits = 4), " ms")
println("ratio KA / threads:       ", round(t_ka / t_threads; sigdigits = 3))
