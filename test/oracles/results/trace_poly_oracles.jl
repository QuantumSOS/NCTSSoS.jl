# Trace Polynomial Oracle Results (6.1, 6.2.0, 6.2.1, 6.2.2)
# ===========================================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/NCTSSOS && julia --project path/to/nctssos_trace_poly.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: Trace_{section}_{sparsity}_d{order}
#
# NOTE: Trace polynomials use the tr() operator and ptraceopt_first API.
# Support format: [[w1], [w2]] means tr(w1)·tr(w2)

# PLACEHOLDER VALUES - Run nctssos_trace_poly.jl on server to generate actual values
# Expected optimal values from NCTSSoS.jl tests:
#   6.1:   ≈ -0.0467 (d=2), ≈ -0.03125 (d=3)
#   6.2.0: -2√2 ≈ -2.8284
#   6.2.1: -4.0
#   6.2.2: -5.0

const TRACE_POLY_ORACLES = Dict(
    # 6.1: Projector Algebra
    # Objective: tr(x₁x₂x₃) + tr(x₁x₂)·tr(x₃)
    "Trace_6_1_Dense_d2" => (opt=-0.046717378455438933, sides=[Int[]], nuniq=0),  # PLACEHOLDER
    "Trace_6_1_Dense_d3" => (opt=-0.03124998978001017, sides=[Int[]], nuniq=0),   # PLACEHOLDER
    "Trace_6_1_TS_d2" => (opt=-0.046717378455438933, sides=[Int[]], nuniq=0),     # PLACEHOLDER
    "Trace_6_1_TS_d3" => (opt=-0.03124998978001017, sides=[Int[]], nuniq=0),      # PLACEHOLDER

    # 6.2.0: CHSH Trace Polynomial (expected: -2√2 ≈ -2.8284)
    # Objective: -tr(x₁y₁) - tr(x₁y₂) - tr(x₂y₁) + tr(x₂y₂)
    "Trace_6_2_0_Dense_d1" => (opt=-2.8284271247461903, sides=[Int[]], nuniq=0),  # PLACEHOLDER
    "Trace_6_2_0_TS_d1" => (opt=-2.8284271247461903, sides=[Int[]], nuniq=0),     # PLACEHOLDER

    # 6.2.1: Squared Trace Expressions (expected: -4.0)
    # Objective: -(tr(x₁y₂) + tr(x₂y₁))² - (tr(x₁y₁) - tr(x₂y₂))²
    "Trace_6_2_1_Dense_d2" => (opt=-4.0, sides=[Int[]], nuniq=0),  # PLACEHOLDER
    "Trace_6_2_1_TS_d2" => (opt=-4.0, sides=[Int[]], nuniq=0),     # PLACEHOLDER

    # 6.2.2: Covariance Trace Polynomial (expected: -5.0)
    # Objective: Σ±cov(i,j) where cov(a,b) = tr(xₐyᵦ) - tr(xₐ)tr(yᵦ)
    "Trace_6_2_2_Dense_d2" => (opt=-5.0, sides=[Int[]], nuniq=0),  # PLACEHOLDER
    "Trace_6_2_2_TS_d2" => (opt=-5.0, sides=[Int[]], nuniq=0),     # PLACEHOLDER
)
