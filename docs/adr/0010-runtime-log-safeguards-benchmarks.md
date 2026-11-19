# Runtime log, run safeguards, benchmarks

Date: 2026-11-19. Status: accepted.

Every run writes `run.log` with TimerOutputs.jl. By default it times the full step, the I/O per output class and the diagnostics per output class; it writes the table periodically and at the end. On a GPU an accurate per-kernel timing needs a device synchronisation after each kernel, which slows the run, so per-kernel timings (FFTs, nonlinear term, truncation, linear step) come only with `profile = true`.

Safeguards: `max_walltime` writes a final snapshot and stops before a cluster job limit. Every `n_check` steps the driver checks for NaN and for CFL > 1; on failure it writes the last snapshot and stops with an error.

The benchmark suite measures time per step and memory for each `(model, scheme, N, device, T)`. It fits `t_step = a·N² log N + b·N²` and predicts the wall time and memory of a target run from `N`, the end time and the CFL-limited `Δt`. It prints the prediction and writes it as TOML. A fit made on one machine does not transfer to another: rerun the short benchmark on a cluster node before extrapolating.
