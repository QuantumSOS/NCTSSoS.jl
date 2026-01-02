# Generalized Rosenbrock Oracle Results
# ======================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/projects/NCTSSOS && julia --project path/to/nctssos_rosenbrock.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: Rosenbrock_{sparsity}_n{vars}_d{order}
#
# Note: Expected optimal is 1.0 at xᵢ = 0 for all i
# The problem has chain-like sparsity: each term involves (xᵢ₋₁, xᵢ)

# Generated 2026-01-02 from NCTSSOS on a800 with MosekTools
const ROSENBROCK_ORACLES = Dict(
    "Rosenbrock_Dense_n6_d2" => (opt=0.999999995930163, sides=[43], nuniq=820),
    "Rosenbrock_CS_n6_d2" => (opt=0.999999973478842, sides=[7, 7, 7, 7, 7], nuniq=90),
    "Rosenbrock_CS_TS_n6_d2" => (opt=0.9999997821660428, sides=[3, 2, 2, 2, 1, 3, 2, 2, 2, 1, 3, 2, 2, 2, 1, 3, 2, 2, 2, 1, 3, 2, 2, 1], nuniq=33),
)
