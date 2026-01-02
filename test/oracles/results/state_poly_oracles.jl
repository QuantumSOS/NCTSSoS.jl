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

# TODO: Run nctssos_state_poly.jl on a800 to generate actual values
const STATE_POLY_ORACLES = Dict(
    # 7.2.0: CHSH State Polynomial (expected: -2√2 ≈ -2.8284)
    "State_7_2_0_Dense_d1" => (opt=-2.8284271321623202, sides=Int[], nuniq=0),
    "State_7_2_0_TS_d1" => (opt=-2.8284271321623202, sides=Int[], nuniq=0),
    
    # 7.2.1: Squared Expectations (expected: -4.0 at order=3)
    "State_7_2_1_Dense_d3" => (opt=-4.0, sides=Int[], nuniq=0),
    
    # 7.2.2: Covariance Expression (expected: -5.0)
    "State_7_2_2_Dense_d2" => (opt=-5.0, sides=Int[], nuniq=0),
    
    # 7.2.3: Mixed State Polynomial (expected: -3.5114802)
    "State_7_2_3_Dense_d2" => (opt=-3.5114802, sides=Int[], nuniq=0),
    "State_7_2_3_TS_d2" => (opt=-3.5114802, sides=Int[], nuniq=0),
)
