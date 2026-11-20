# Plotting functions: arguments, library selection and the acceptance run

Date: 2026-11-20. Status: accepted.

WP-15 fixes the parts of ADR-0012 that leave the plotting interface open, and the parameters of the user acceptance run of ADR-0011. Changes:

- **Files, not the simulation.** The plotting functions read the output files of a run, so a user plots a finished run without the `Simulation` in memory. Each takes a file or an output directory: `plot_field(path; field)` takes a snapshot or the output directory (then the snapshot of the highest index), `plot_spectrum(path; names, time)` takes `spectral.nc` or the output directory, `plot_global(path; names)` takes `global.nc` or the output directory, and `animate_snapshots(outdir, file; field, spectrum, framerate, decades)` takes the output directory, because it reads the snapshots and `spectral.nc`.
- **Vorticity from the snapshot.** A snapshot holds `ux`, `uy` (and `rho`), not `ω`. `field = :omega` computes `ω = ∂x u_y − ∂y u_x` spectrally on the CPU with a `Grid` and `Transforms` built from the file, so the derivative is the one of the solver.
- **Library selection.** Both extensions add methods to the same generic functions. Two methods with the same signature would overwrite each other when CairoMakie and Plots are both loaded. The extensions therefore add methods for `Val(:Makie)` and `Val(:Plots)` as the first argument, and the public method dispatches on the keyword `library`. Its default is the only loaded extension; with none or both loaded the call throws an `ArgumentError`. The extensions return a `Makie.Figure` or a `Plots.Plot`, so the user changes the figure with the library's own functions.
- **Shared readers in the core.** The readers of the files (snapshot field, spectra, GLOBAL diagnostics, colour range, frame data) are internal functions of `src/plotting.jl`. They use only NCDatasets and the solver, so the two extensions hold only drawing code and stay symmetric.
- **Spectra only.** `plot_spectrum` draws variables with the coordinate `shell_k` on log–log axes and leaves out the shell `k = 0` and every non-positive value. A flux (coordinate `shell_edge`) takes both signs and is rejected (`ArgumentError`); a flux plot is in `docs/BACKLOG.md`.
- **Animation frames.** One frame per snapshot. The spectrum of a frame is the SPECTRAL record nearest in time to the snapshot; a run that sets the SPECTRAL period equal to the FIELDS period gets matching times. All frames share the colour range (extrema over all snapshots, symmetric about zero when the field takes both signs) and the spectrum axis (`decades` decades below the largest value), so the evolution of the spectrum tail is visible. A vertical line marks `k_max`. The extension of `file` selects the format.
- **Tests.** `test/plotting.jl` is in the blocking tier, CPU backend only. It loads CairoMakie and Plots (GR with `GKSwstype = 100`, no display). Loading and compiling both libraries adds about 80 s to `Pkg.test`, below the 3 min limit of WP-15.
- **Acceptance run.** `test/acceptance/crossed_kolmogorov.jl`: incompressible, unforced, `A₁ = A₂ = 1`, `k₁ = 2`, `k₂ = 3`, `N = 512`, `Lx = Ly = 2π` (`k_max = 170`), Float64, RK4 with the integrating factor, `Δt = 10⁻³`, to `t = 30`; FIELDS and SPECTRAL every 0.25, GLOBAL every 0.05. The script fixes `ν = 4·10⁻⁴`.
- **Choice of `ν` (proposed amendment of ADR-0011; awaits the user's decision).** A scan at `N = 512` to `t = 10` measured, at each `ν`, the palinstrophy peak, `k_max/k_ν` with `k_ν = (2νP/ν³)^{1/6}`, and the ringing: the largest grid-point value of the part of `ω` with `|k| > 0.9 k_max`, relative to `max|ω|`, over the snapshots.

  | `ν` | `k_max/k_ν` | `Re = A₁/(k₁ν)` | ringing |
  |-----|-------------|-----------------|---------|
  | 6·10⁻⁴ | 5.1 | 833 | 0.05% |
  | 4·10⁻⁴ | 4.3 | 1250 | 0.5% |
  | 3·10⁻⁴ | 3.8 | 1667 | 1.3% |
  | 2·10⁻⁴ | 3.2 | 2500 | 2.9% |
  | 10⁻⁴ | 2.3 | 5000 | 7.3% |
  | 5·10⁻⁵ | 1.7 | 10⁴ | 13% |
  | 4·10⁻⁵ | 1.50 | 1.25·10⁴ | — |

  The rule of ADR-0011 (`k_max/k_ν ≥ 1.5`) selects `ν = 5·10⁻⁵`. That run shows grid-scale stripes in `ω` at the fronts of the vorticity sheets (Gibbs oscillations of the truncation), visible in the animation, although its shell spectrum at `k = 160` stays 8 decades below the peak. The energy spectrum hides them because the stripes are weighted by `k²` in `ω` and add coherently over the modes of a shell. The rule becomes: the smallest `ν` of the scan with ringing below 1% and `Re ≥ 100`. This gives `ν = 4·10⁻⁴`: `k_max/k_ν = 4.3`, `Re = 1250`, the spectrum at `k = 160` at most `4·10⁻¹¹` of its peak over the run, and the nonlinear roll-up of the shear layers into vortices.
