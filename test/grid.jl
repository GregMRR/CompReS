# Tests of the grid (WP-1): wavenumbers, circular truncation and shells, checked on hand-counted examples.
module GridTests

using Test
using KernelAbstractions
using CompReS

# Array position of the mode with integer indices (i, j) in the real-FFT layout.
pos(g, i, j) = (i + 1, j < 0 ? g.Ny + j + 1 : j + 1)

function run_tests(backend)
    @testset "grid" begin
        @testset "wavenumbers, truncation and shells for N = 8" begin
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 8, backend)
            @test size(g.kx) == (5, 8)
            # Lx = 2π gives k = integer; ky follows the FFT order 0, 1, 2, 3, -4, -3, -2, -1.
            @test Array(g.kx)[:, 1] == [0, 1, 2, 3, 4]
            @test Array(g.ky)[1, :] == [0, 1, 2, 3, -4, -3, -2, -1]
            @test Array(g.k2)[pos(g, 2, -1)...] ≈ 5
            # k_c = 8/3, so |k|² < 7.11 keeps |k|² in {0, 1, 2, 4, 5}; the largest is 5.
            @test g.k_max ≈ sqrt(5)
            @test g.Δk ≈ 1
            # Full-plane counts: 21 modes, 1 + 4 + 4 + 4 + 8.
            # The grid corner (4, -4) has |k| = 5.66, in shell 6, so the shell vectors have 7 entries.
            @test g.shell_count == [1, 8, 12, 0, 0, 0, 0]
            @test g.shell_k ≈ 0:6
            @test g.shell_edge ≈ (0:6) .+ 0.5
            # All 64 grid modes count once in the full plane: the kx = 0 and Nyquist columns have weight 1.
            @test sum(g.shell_count_all) == 64
            # Only the corner mode (4, -4) lies in shell 6; the Nyquist column does not double it.
            @test g.shell_count_all[7] == 1
            @test Array(g.shell_index)[pos(g, 4, -4)...] == 7
            # Half-layout count: kx = 0 column has 5 modes, kx = 1 has 5, kx = 2 has 3 (21 = 5 + 2·(5 + 3)).
            @test sum(Array(g.truncation)) == 13
            @test Array(g.truncation)[pos(g, 2, 1)...]
            @test !Array(g.truncation)[pos(g, 2, 2)...]
            # Shell index is n + 1 for every mode; mode (2, 2) has |k| = 2.83, shell 3. Mode (2, 1) has |k| = 2.24, shell 2. Discarded modes keep their shell.
            @test Array(g.shell_index)[pos(g, 2, 1)...] == 3
            @test Array(g.shell_index)[pos(g, 0, 0)...] == 1
            @test Array(g.shell_index)[pos(g, 2, 2)...] == 4
        end

        @testset "N = 12 keeps no mode at |k| = 4" begin
            # k_c = 4 and the bound is strict, so (4, 0), (0, 4) and (0, -4) are out.
            g = Grid(; Lx = 2π, Ly = 2π, Nx = 12, Ny = 12, backend)
            tr = Array(g.truncation)
            @test !tr[pos(g, 4, 0)...]
            @test !tr[pos(g, 0, 4)...]
            @test !tr[pos(g, 0, -4)...]
            @test tr[pos(g, 3, 2)...]
            @test g.k_max ≈ sqrt(13)
        end

        @testset "rectangle Lx = 4π, Ly = 2π" begin
            # Equal spacing needs Nx = 16, Ny = 8; Δk = 1/2 and k_c = 8/3.
            g = Grid(; Lx = 4π, Ly = 2π, Nx = 16, Ny = 8, backend)
            @test Array(g.kx)[:, 1] ≈ (0:8) ./ 2
            # Retained: i² + 4j² ≤ 28, 43 modes counted by hand; the largest has i² + 4j² = 25.
            @test g.k_max ≈ 2.5
            @test g.Δk ≈ 0.5
            @test g.shell_count[1:6] == [1, 2, 8, 6, 20, 6]
            @test sum(g.shell_count) == 43
            @test sum(g.shell_count_all) == 16 * 8
        end

        @testset "mode on a shell edge goes to the upper shell" begin
            # Mode (6, -14): |k| = sqrt(1.2² + 3.5²) = 3.7 = 18.5 Δk with Δk = 0.2, so shell 19.
            g = Grid(; Lx = 5 * 2π, Ly = 4 * 2π, Nx = 60, Ny = 48, backend)
            @test Array(g.truncation)[pos(g, 6, -14)...]
            @test Array(g.shell_index)[pos(g, 6, -14)...] == 19 + 1
            @test g.shell_edge[19] ≈ 3.7
        end

        @testset "CPU is re-exported" begin
            @test CompReS.CPU === KernelAbstractions.CPU
        end

        @testset "construction checks" begin
            # Unequal spacing: 2π/8 against 2π/10.
            @test_throws ArgumentError Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 10, backend)
            # Odd N is rejected.
            @test_throws ArgumentError Grid(; Lx = 2π, Ly = 2π, Nx = 9, Ny = 9, backend)
            @test_throws ArgumentError Grid(; Lx = 2π, Ly = 2π, Nx = 8, Ny = 9, backend)
            # Spacing 0.01 on both axes passes for both float types, whatever the rounding of Lx / Nx.
            for T in (Float64, Float32)
                g = Grid(; Lx = 0.3, Ly = 0.1, Nx = 30, Ny = 10, T, backend)
                @test eltype(g.kx) == T
                @test eltype(g.k2) == T
                @test g.k_max isa T
            end
        end

        @testset "no mode on the disc boundary when 3 divides N" begin
            # Mode (0, 30) has |k| = k_c exactly: 9·30²·30² = 30²·90², so the strict test drops it.
            g = Grid(; Lx = 0.3, Ly = 0.9, Nx = 30, Ny = 90, backend)
            @test !Array(g.truncation)[pos(g, 0, 30)...]
            @test Array(g.truncation)[pos(g, 0, 29)...]
        end
    end
end

end
