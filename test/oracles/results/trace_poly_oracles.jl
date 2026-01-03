# Trace Polynomial Oracle Results (6.1, 6.2.0, 6.2.1, 6.2.2)
# ===========================================================
# Generated from NCTSSOS with MosekTools
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

# Generated 2026-01-03 from NCTSSOS with MosekTools
const TRACE_POLY_ORACLES = Dict(
    # 6.1: Projector Algebra
    # Objective: tr(x₁x₂x₃) + tr(x₁x₂)·tr(x₃)
    # Note: TS_d2 gives 0.0 (too sparse), TS_d3 gives worse bound
    "Trace_6_1_Dense_d2" => (opt=-0.04671737845552321, sides=[31], nuniq=81),
    "Trace_6_1_Dense_d3" => (opt=-0.031249989780027937, sides=[108], nuniq=395),
    "Trace_6_1_TS_d2" => (opt=0.0, sides=[5, 5, 5, 5, 4, 4, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=39),
    "Trace_6_1_TS_d3" => (opt=-0.27196633115451874, sides=[7, 7, 7, 7, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=94),

    # 6.2.0: CHSH Trace Polynomial (expected: -2√2 ≈ -2.8284)
    # Objective: -tr(x₁y₁) - tr(x₁y₂) - tr(x₂y₁) + tr(x₂y₂)
    "Trace_6_2_0_Dense_d1" => (opt=-2.8284271247462307, sides=[9], nuniq=21),
    "Trace_6_2_0_TS_d1" => (opt=-2.8284271247321175, sides=[3, 3, 1, 1, 1, 1, 1], nuniq=10),

    # 6.2.1: Squared Trace Expressions (expected: -4.0)
    # Objective: -(tr(x₁y₂) + tr(x₂y₁))² - (tr(x₁y₁) - tr(x₂y₂))²
    # Note: TS gives -8.0 (looser bound due to sparsity structure)
    "Trace_6_2_1_Dense_d2" => (opt=-4.000000007251562, sides=[53], nuniq=222),
    "Trace_6_2_1_TS_d2" => (opt=-8.00000000078321, sides=[3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=23),

    # 6.2.2: Covariance Trace Polynomial (expected: -5.0)
    # Objective: Σ±cov(i,j) where cov(a,b) = tr(xₐyᵦ) - tr(xₐ)tr(yᵦ)
    "Trace_6_2_2_Dense_d2" => (opt=-4.999999997079296, sides=[115], nuniq=1010),
    "Trace_6_2_2_TS_d2" => (opt=-4.999999995608242, sides=[25, 19, 19, 19, 17, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 9, 9, 8, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=199),
)
