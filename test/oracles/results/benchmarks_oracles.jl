# NC Polynomial Benchmarks Oracle Results
# ========================================
# Generated from NCTSSOS with MosekTools
# Run: cd ~/NCTSSOS && julia --project path/to/nctssos_benchmarks.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: {Benchmark}_{sparsity}_n{n}_d{order}
#
# Benchmarks included:
#   - Broyden Banded (n=4, degree=6, d=3)
#   - Broyden Tridiagonal (n=6, degree=4, d=2)
#   - Chained Singular (n=8, degree=4, d=2)
#   - Chained Wood (n=8, degree=4, d=2)
#
# Note: Rosenbrock oracles are in rosenbrock_oracles.jl

# Generated 2026-01-03 from NCTSSOS with MosekTools
const BENCHMARKS_ORACLES = Dict(
    "BroydenBanded_CS_TS_n4_d3" => (opt=-1.7480268537791176e-9, sides=[9, 9, 7, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=134),
    "BroydenBanded_CS_n4_d3" => (opt=-2.0628417282628906e-8, sides=[85], nuniq=2815),
    "BroydenBanded_Dense_n4_d3" => (opt=-2.0628417282628906e-8, sides=[85], nuniq=2815),
    "BroydenTridiagonal_CS_TS_n6_d2" => (opt=-1.3851129121677345e-7, sides=[5, 4, 4, 3, 3, 3, 3, 2, 2, 5, 4, 4, 3, 3, 3, 3, 2, 2, 5, 4, 4, 3, 3, 3, 3, 2, 2, 5, 4, 4, 3, 3, 3, 3, 2, 2], nuniq=62),
    "BroydenTridiagonal_CS_n6_d2" => (opt=-7.683336739295744e-9, sides=[13, 13, 13, 13], nuniq=226),
    "BroydenTridiagonal_Dense_n6_d2" => (opt=6.580588236100922e-11, sides=[43], nuniq=820),
    "ChainedSingular_CS_TS_n8_d2" => (opt=-2.4177846987055753e-7, sides=[4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1], nuniq=83),
    "ChainedSingular_CS_n8_d2" => (opt=-1.2251940085171628e-7, sides=[13, 13, 13, 13, 13, 13], nuniq=328),
    "ChainedSingular_Dense_n8_d2" => (opt=-6.428097254420157e-8, sides=[73], nuniq=2413),
    "ChainedWood_CS_TS_n8_d2" => (opt=0.9999998815901131, sides=[3, 2, 2, 2, 1, 3, 2, 2, 2, 2, 3, 2, 2, 2, 1, 3, 2, 2, 2, 2, 3, 2, 2, 2, 1, 3, 2, 2, 2, 2, 3, 2, 2, 2, 1], nuniq=46),
    "ChainedWood_CS_n8_d2" => (opt=0.9999998300594172, sides=[7, 7, 7, 7, 7, 7, 7], nuniq=124),
    "ChainedWood_Dense_n8_d2" => (opt=0.9999998013387236, sides=[73], nuniq=2413),
)
