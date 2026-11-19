# Blocking test runner: runs every test file on the CPU and, when functional, on CUDA.
#
# Test file convention: `test/<name>.jl` defines a module `<Name>Tests` with a function
# `run_tests(backend)` that contains its own `@testset`. This file includes the test files
# and calls `run_tests` once per backend.
#
# Name filter: `julia --project=test test/runtests.jl <name>` runs only `test/<name>.jl`.
using Test
using KernelAbstractions
using CUDA

backends() = CUDA.functional() ? [CPU(), CUDABackend()] : [CPU()]

is_cpu(backend) = backend isa CPU

const TEST_NAMES = sort!([splitext(f)[1] for f in readdir(@__DIR__)
                          if endswith(f, ".jl") && f != "runtests.jl"])

selected = isempty(ARGS) ? TEST_NAMES : filter(in(ARGS), TEST_NAMES)
isempty(selected) && error("no test file matches $(ARGS); available: $TEST_NAMES")

for name in selected
    include(joinpath(@__DIR__, "$name.jl"))
    mod = getfield(Main, Symbol(uppercasefirst(name), "Tests"))
    for backend in backends()
        @testset "$name on $(nameof(typeof(backend)))" begin
            Base.invokelatest(mod.run_tests, backend)
        end
    end
end
