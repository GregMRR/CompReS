# Trivial test that checks the harness: a KernelAbstractions kernel fills an array on each backend.
module SkeletonTests

using Test
using KernelAbstractions

@kernel function fill_index!(a)
    i = @index(Global)
    a[i] = i
end

function run_tests(backend)
    @testset "kernel fills an array" begin
        n = 16
        a = KernelAbstractions.zeros(backend, Float64, n)
        fill_index!(backend)(a; ndrange = n)
        KernelAbstractions.synchronize(backend)
        @test Array(a) == collect(1.0:n)
    end
end

end
