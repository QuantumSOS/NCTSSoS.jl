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

# TODO: Run nctssos_rosenbrock.jl on a800 to generate actual values
const ROSENBROCK_ORACLES = Dict(
    "Rosenbrock_Dense_n6_d2" => (opt=1.0, sides=Int[], nuniq=0),
    "Rosenbrock_CS_n6_d2" => (opt=1.0, sides=Int[], nuniq=0),
    "Rosenbrock_CS_TS_n6_d2" => (opt=1.0, sides=Int[], nuniq=0),
)
