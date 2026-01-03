# State Polynomial Extended Oracle Results (7.2.1, 7.2.2, 7.2.3)
# ================================================================
# Generated from NCTSSOS with MosekTools
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

# Generated 2026-01-03 from NCTSSOS with MosekTools
const STATE_POLY_EXTENDED_ORACLES = Dict(
    # 7.2.1: Squared Expectations (expected: -4.0)
    # Objective: -(ς(x₁y₂) + ς(x₂y₁))² - (ς(x₁y₁) - ς(x₂y₂))²
    # Note: TS gives -8.0 (looser bound due to sparsity)
    "State_7_2_1_Dense_d3" => (opt=-3.9999999914666895, sides=[209], nuniq=1887),
    "State_7_2_1_TS_d3" => (opt=-7.9999999998214975, sides=[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=75),

    # 7.2.2: Covariance (expected: -5.0)
    # Objective: Σ±cov(i,j) where cov(a,b) = ς(xₐyᵦ) - ς(xₐ)ς(yᵦ)
    "State_7_2_2_Dense_d2" => (opt=-4.999999999824081, sides=[106], nuniq=1098),
    "State_7_2_2_TS_d2" => (opt=-4.99999999745226, sides=[9, 9, 9, 9, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=93),

    # 7.2.3: Mixed State Polynomial (expected: -3.5114802)
    # Objective: Mixed linear, product, and squared expectations
    # Note: TS gives slightly worse bound (-3.582)
    "State_7_2_3_Dense_d2" => (opt=-3.511480225797076, sides=[49], nuniq=233),
    "State_7_2_3_TS_d2" => (opt=-3.582132180463948, sides=[7, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1], nuniq=41),
)
