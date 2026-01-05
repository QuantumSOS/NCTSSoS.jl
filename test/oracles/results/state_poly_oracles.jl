# State Polynomial Oracle Results
# ================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/projects/NCTSSOS && julia --project path/to/nctssos_state_poly.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: State_{section}_{sparsity}_d{order}
#
# Note: State polynomial optimization uses stateopt_first in NCTSSOS
# which handles the state expectation structure (ς operator).

# Generated 2026-01-02 from NCTSSOS on a800 with MosekTools
# Note: Only pure state expectation problems can use pstateopt_first.
# Problems 7.2.1, 7.2.2, 7.2.3 require mixword API for product terms and are not tested.
const STATE_POLY_ORACLES = Dict(
    # 7.2.0: CHSH State Polynomial (expected: -2√2 ≈ -2.8284)
    "State_7_2_0_Dense_d1" => (opt=-2.828427124746231, sides=[11], nuniq=34),
    "State_7_2_0_TS_d1" => (opt=-2.828427124732117, sides=[3, 3, 1, 1, 1, 1, 1, 1, 1], nuniq=12),
    
    # 7.2.1, 7.2.2, 7.2.3: Require mixword API - not testable with current script
)
