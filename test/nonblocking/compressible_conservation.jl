# Non-blocking test of the inviscid invariants of the compressible model (ADR-0011, specification §4): unforced,
# μ = ζ = κ = 0, RK4 without integrating factor. E_tot and Omega drift for two reasons:
#   - the time scheme, whose error falls as Δt^4 or faster;
#   - the aliasing floor: the products with ρ, 1/ρ and √ρ are not band-limited, so even the exact-in-time
#     truncated equations do not conserve E_tot and Omega. The floor falls with N while the flow is resolved.
# The file runs on the CPU, prints the measured values and writes the figure compressible_conservation.png
# next to this file for visual inspection.
module CompressibleConservationNonblocking

using Test
using CompReS
import CairoMakie as CM

const INVARIANTS = Dict(select_diagnostics(:GLOBAL, (:E_tot, :Omega, :mass)))
const T_CHECK = 1.0     # time of the measurements; the flow is resolved on every grid up to about t = 1.2
const T_HISTORY = 2.0   # end of the drift histories; after t ≈ 1.5 the inviscid flow forms grid-scale fronts
const CS = 5.0          # sound speed; Ma ≈ 0.2

# Smooth, band-limited velocity (|k| ≤ 3), the same function on every grid.
ux(x, y) = 0.7sin(3y) + 0.5sin(2x) + 0.3cos(x + y)
uy(x, y) = 0.4cos(2x - y) + 0.6sin(x)

# Inviscid compressible problem on an N × N grid of [0, 2π)², with the balanced density of the velocity.
function setup(N)
    g = Grid(; Lx = 2π, Ly = 2π, Nx = N, Ny = N)
    x = [2π * (i - 1) / N for i in 1:N]
    ws0 = Workspace(g, State(kolmogorov(g; A = 1, n = 1)))
    spectral(f) = to_spectral!(similar(ws0.forcing.ux), f.(x, x'), ws0.transforms)
    eos = Isothermal(; cs = CS)
    u = (ux = spectral(ux), uy = spectral(uy))
    state = State((ux = u.ux, uy = u.uy, rho = balanced_density(u, eos, g).rho))
    problem = Problem(g, Compressible(; eos, μ = 0, ζ = 0), NoForcing())
    return state, problem, Workspace(g, state; scheme = RK4(), integrating_factor = false)
end

invariants(state, problem, ws) = [INVARIANTS[name](state, problem, ws) for name in (:E_tot, :Omega, :mass)]

"""
Run to `t_end` with the step `Δt`. Return the times and the relative changes `(q(t) − q(0))/q(0)` of
E_tot, Omega and mass (rows), sampled every `every`.
"""
function drift_history(N, Δt, t_end; every)
    state, problem, ws = setup(N)
    q0 = invariants(state, problem, ws)
    m = round(Int, every / Δt)
    times, drift = [0.0], zeros(3, 1)
    for n in 1:round(Int, t_end / Δt)
        step!(state, problem, ws, Δt)
        if n % m == 0
            push!(times, n * Δt)
            drift = hcat(drift, (invariants(state, problem, ws) .- q0) ./ q0)
        end
    end
    return times, drift
end

final_drift(N, Δt) = drift_history(N, Δt, T_CHECK; every = T_CHECK)[2][:, end]

# Least-squares slope of log(err) against log(Δt).
function loglog_slope(x, y)
    lx, ly = log.(x), log.(y)
    return sum((lx .- sum(lx) / length(lx)) .* ly) / sum((lx .- sum(lx) / length(lx)) .^ 2)
end

# Four panels: time-scheme error against Δt, maximum drift up to T_CHECK against N, drift histories of E_tot and Omega.
function plot_conservation(Δts, time_error, Ns, floors, histories)
    fig = CM.Figure(; size = (1200, 900), fontsize = 16)
    clip(v) = max.(abs.(v), 1e-17)
    names = ("E_tot", "Omega")

    ax = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "Δt", ylabel = "|q(Δt) − q(Δt_ref)| / q(0)",
                 title = "RK4 time-scheme error, N = 64, t = $T_CHECK")
    for (i, name) in enumerate(names)
        CM.scatterlines!(ax, Δts, clip(time_error[i, :]); label = name)
    end
    CM.lines!(ax, Δts, time_error[1, 1] .* (Δts ./ Δts[1]) .^ 4; color = :gray, linestyle = :dash, label = "Δt⁴")
    CM.axislegend(ax; position = :rb)

    ax = CM.Axis(fig[1, 2]; xscale = log2, yscale = log10, xlabel = "N", ylabel = "max |q(t) − q(0)| / q(0)",
                 xticks = (Ns, string.(Ns)), title = "Maximum drift for t ≤ $T_CHECK, Δt = 5·10⁻⁴ (aliasing floor)")
    for (i, name) in enumerate(names)
        CM.scatterlines!(ax, Ns, clip(floors[i, :]); label = name)
    end
    CM.axislegend(ax; position = :rt)

    for (i, name) in enumerate(names)
        ax = CM.Axis(fig[2, i]; yscale = log10, xlabel = "t", ylabel = "|q(t) − q(0)| / q(0)",
                     title = "Drift of $name, Δt = 5·10⁻⁴")
        for (N, (times, drift)) in zip(Ns, histories)
            CM.lines!(ax, times[2:end], clip(drift[i, 2:end]); label = "N = $N")
        end
        CM.vlines!(ax, [T_CHECK]; color = :gray, linestyle = :dash)
        CM.axislegend(ax; position = :rb)
    end
    return fig
end

@testset "non-blocking: compressible inviscid invariants" begin
    # Time-scheme error at N = 64: the difference to a run with a much smaller step removes the aliasing floor.
    Δts = [8e-3, 4e-3, 2e-3, 1e-3]
    reference = final_drift(64, 1.25e-4)
    time_error = reduce(hcat, [abs.(final_drift(64, Δt) .- reference) for Δt in Δts])
    slopes = map(1:2) do i
        resolved = time_error[i, :] .> 1e-13   # points above round-off
        loglog_slope(Δts[resolved], time_error[i, resolved])
    end
    println("RK4 time-scheme error at N = 64, t = $T_CHECK, Δt = $Δts")
    println("  E_tot: ", time_error[1, :], "  slope ", slopes[1])
    println("  Omega: ", time_error[2, :], "  slope ", slopes[2])
    @testset "time-scheme error of $name falls at least as Δt^3.8" for (i, name) in enumerate(("E_tot", "Omega"))
        @test slopes[i] > 3.8
    end

    # Aliasing floor against N: Δt = 5e-4 keeps the time-scheme error near 1e-12 at t = 1.
    Ns = [32, 64, 128, 256]
    histories = [drift_history(N, 5e-4, T_HISTORY; every = 0.05) for N in Ns]
    i_check = findfirst(≈(T_CHECK), histories[1][1])
    # The maximum over 0 ≤ t ≤ T_CHECK: the drift at one time can pass through zero.
    floors = reduce(hcat, [vec(maximum(abs, d[:, 1:i_check]; dims = 2)) for (_, d) in histories])
    println("Maximum relative drift for t ≤ $T_CHECK, Δt = 5e-4, N = $Ns")
    println("  E_tot: ", floors[1, :])
    println("  Omega: ", floors[2, :])
    println("  mass:  ", floors[3, :])
    @test all(iszero, floors[3, :])
    @testset "maximum drift of $name for t ≤ $T_CHECK falls with N" for (i, name) in enumerate(("E_tot", "Omega"))
        @test issorted(floors[i, :]; rev = true)
        @test floors[i, end] < 1e-9
    end

    figure = plot_conservation(Δts, time_error, Ns, floors, histories)
    path = joinpath(@__DIR__, "compressible_conservation.png")
    CM.save(path, figure)
    println("figure: ", path)
end

end
