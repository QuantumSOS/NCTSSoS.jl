# State Polynomial Extended Oracle Results (7.2.1, 7.2.2, 7.2.3)
# ================================================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/NCTSSOS && julia --project path/to/nctssos_state_poly_extended.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: State_{section}_{sparsity}_d{order}
#
# NOTE: These problems use product terms ς(w1)·ς(w2), requiring the
# multi-word support format in pstateopt_first.

# PLACEHOLDER VALUES - Run nctssos_state_poly_extended.jl on server to generate actual values
# Expected optimal values from NCTSSoS.jl tests:
#   7.2.1: -4.0 (at order=3)
#   7.2.2: -5.0 (at order=2)
#   7.2.3: -3.5114802 (at order=2)

const STATE_POLY_EXTENDED_ORACLES = Dict(
    # 7.2.1: Squared Expectations (expected: -4.0)
    # Objective: -(ς(x₁y₂) + ς(x₂y₁))² - (ς(x₁y₁) - ς(x₂y₂))²
    "State_7_2_1_Dense_d3" => (opt=-4.0, sides=[Int[]], nuniq=0),  # PLACEHOLDER
    "State_7_2_1_TS_d3" => (opt=-4.0, sides=[Int[]], nuniq=0),     # PLACEHOLDER

    # 7.2.2: Covariance (expected: -5.0)
    # Objective: Σ±cov(i,j) where cov(a,b) = ς(xₐyᵦ) - ς(xₐ)ς(yᵦ)
    "State_7_2_2_Dense_d2" => (opt=-5.0, sides=[Int[]], nuniq=0),  # PLACEHOLDER
    "State_7_2_2_TS_d2" => (opt=-5.0, sides=[Int[]], nuniq=0),     # PLACEHOLDER

    # 7.2.3: Mixed State Polynomial (expected: -3.5114802)
    # Objective: Mixed linear, product, and squared expectations
    "State_7_2_3_Dense_d2" => (opt=-3.5114802, sides=[Int[]], nuniq=0),  # PLACEHOLDER
    "State_7_2_3_TS_d2" => (opt=-3.5114802, sides=[Int[]], nuniq=0),     # PLACEHOLDER
)
