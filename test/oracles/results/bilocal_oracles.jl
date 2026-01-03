# Bilocal Network Oracle Results (Example 8.1.1)
# ===============================================
# Generated from NCTSSOS with MosekTools
# Run: cd ~/NCTSSOS && julia --project path/to/nctssos_bilocal.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: Bilocal_8_1_1_{sparsity}_d{order}
#
# Problem: Bilocal network quantum bound from Example 8.1.1
# Variables: x[1:2], y[1:2], z[1:2] (6 vars total)
#   - Alice (x): 2 measurement operators
#   - Bob (y): 2 measurement operators
#   - Charlie (z): 2 measurement operators
# Constraint: unipotent (U²=I)
# Expected quantum bound: -4.0

# Generated 2026-01-03 from NCTSSOS with MosekTools
# Note: All sparsity variants give block sizes of all 1s due to the
# state polynomial structure with product terms.
const BILOCAL_ORACLES = Dict(
    # Dense: 3380 blocks of size 1, 3380 constraints
    "Bilocal_8_1_1_Dense_d3" => (opt=-3.9999999896685883, nblocks=3380, nuniq=3380),
    # CS: Same structure
    "Bilocal_8_1_1_CS_d3" => (opt=-3.9999999896685883, nblocks=3380, nuniq=3380),
    # TS: Reduced constraints
    "Bilocal_8_1_1_TS_d3" => (opt=-3.9999999994438804, nblocks=1596, nuniq=263),
    # CS+TS: Same as TS
    "Bilocal_8_1_1_CS_TS_d3" => (opt=-3.9999999994438804, nblocks=1596, nuniq=263),
)
