# CHSH Bell Inequality Oracle Results
# =====================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/projects/NCTSSOS && julia --project path/to/nctssos_chsh.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (maximizes -f → opt ≈ -2√2 for correct result)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: CHSH_{sparsity}_d{order}
#   Each variant tested at its specified order from CHSH_VARIANTS
#
# NOTE: CS_TS gives opt=-4.0 (not -2.828) - this is NCTSSOS behavior.
#       The combined CS+TS relaxation does NOT converge to -2.828 even at order 6.
#       The cliques share only the identity (nuniq=5 constant), so cross-clique
#       moment constraints are missing. This is a structural limitation for CHSH.

const CHSH_ORACLES = Dict(
    "CHSH_Dense_d1" => (opt=-2.828427132162319, sides=[5], nuniq=11),
    "CHSH_CS_d1" => (opt=-2.8284271247170487, sides=[4, 4], nuniq=10),
    "CHSH_TS_d1" => (opt=-2.828427124732117, sides=[3, 3, 1], nuniq=6),
    "CHSH_CS_TS_d2" => (opt=-3.999999999803648, sides=[3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 2, 2], nuniq=5),
)
