# r4-probe-int128-cuda: does the amended shell rule (ROADMAP:105, Int128 products) compile and run in a
# KernelAbstractions kernel on CUDA? Compares the device shell_index with the CPU one on N = 2048 (square).
# Run: julia --project=<env with CUDA, KernelAbstractions> r4-probe-int128-cuda.jl
using CUDA, KernelAbstractions
@kernel function shells!(s, Nx, Ny)
    I = @index(Global, Cartesian); a = I[1] - 1; b = I[2] - 1; j = b < Ny ÷ 2 ? b : b - Ny
    Nm = Int128(max(Nx, Ny)); A = 4 * Nm^2 * (Int128(a)^2 * Ny^2 + Int128(j)^2 * Nx^2); D = Int128(Nx)^2 * Ny^2
    n = 0
    while A >= Int128(2n + 1)^2 * D; n += 1; end
    s[I] = n
end
N = 2048
for be in (CPU(), CUDABackend())
    s = KernelAbstractions.zeros(be, Int32, N ÷ 2 + 1, N)
    shells!(be)(s, N, N; ndrange = size(s)); KernelAbstractions.synchronize(be)
    global res = be isa CPU ? Array(s) : (println("CUDA == CPU: ", Array(s) == res); nothing)
end
