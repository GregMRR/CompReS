# Non-blocking test runner: runs the tests of this directory, reports failures and always exits 0.
using Test

const TEST_NAMES = sort!([f for f in readdir(@__DIR__) if endswith(f, ".jl") && f != "runtests.jl"])

results = Pair{String,Bool}[]
for file in TEST_NAMES
    ok = try
        include(joinpath(@__DIR__, file))
        true
    catch err
        @warn "non-blocking test failed" file exception = err
        false
    end
    push!(results, file => ok)
end

println("Non-blocking tier: $(count(last, results)) of $(length(results)) files passed.")
foreach(r -> println("  ", last(r) ? "pass  " : "FAIL  ", first(r)), results)
exit(0)
