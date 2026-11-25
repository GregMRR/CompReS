# r4-probe-hratio: ROADMAP:192 states "h(1) is about 10³ times its variation" for the WP-14 set-up
# (cellular flow, u_rms = 1, Ma = 0.1, Polytropic γ = 1.4). Computes h(1)/(max h − min h) on N = 64.
# Run: julia --project=<env with FFTW> r4-probe-hratio.jl
include("r4-probe-balanced.jl")   # reuses setup(); also prints that probe's table
S = setup(64); cs2 = 100.0
println("Polytropic h(1) = ", cs2 / 0.4, ";  range of h = ", maximum(S.h) - minimum(S.h),
        ";  ratio = ", cs2 / 0.4 / (maximum(S.h) - minimum(S.h)), ";  Isothermal h(1) = 0")
