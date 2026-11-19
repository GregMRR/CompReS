# AMD runner: runs the blocking test files on `ROCBackend()`; needs AMD hardware.
using Test
using AMDGPU
using AMDGPU: ROCBackend

const TEST_DIR = joinpath(@__DIR__, "..")
const TEST_NAMES = sort!([splitext(f)[1] for f in readdir(TEST_DIR)
                          if endswith(f, ".jl") && f != "runtests.jl"])

is_cpu(backend) = false

for name in TEST_NAMES
    include(joinpath(TEST_DIR, "$name.jl"))
    mod = getfield(Main, Symbol(uppercasefirst(name), "Tests"))
    @testset "$name on ROCBackend" begin
        Base.invokelatest(mod.run_tests, ROCBackend())
    end
end
