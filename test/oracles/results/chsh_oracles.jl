# CHSH Bell Inequality Oracle Results (All Formulations)
# =======================================================
# Generated from NCTSSOS with MosekTools
# Run: cd ~/projects/NCTSSOS && julia --project path/to/nctssos_chsh.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization → opt ≈ -2√2 for correct result)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# This file contains oracles for THREE formulations:
# 1. NC Polynomial: CHSH_{Dense|CS|TS|CS_TS}_d{order}
# 2. State Polynomial: CHSH_State_{Dense|TS}_d{order}
# 3. Trace Polynomial: CHSH_Trace_{Dense|TS}_d{order}
#
# NOTE: CS_TS gives opt=-4.0 (not -2.828) - this is NCTSSOS behavior.
#       The combined CS+TS relaxation does NOT converge to -2.828 even at order 6.
#       The cliques share only the identity (nuniq=5 constant), so cross-clique
#       moment constraints are missing. This is a structural limitation for CHSH.

const CHSH_ORACLES = Dict(
    # =========================================================================
    # NC Polynomial Formulation (nctssos_first / cs_nctssos_first)
    # =========================================================================
    "CHSH_Dense_d1" => (opt=-2.8284271321623193, sides=[5], nuniq=11),
    "CHSH_CS_d1" => (opt=-2.8284271247170496, sides=[4, 4], nuniq=10),
    "CHSH_TS_d1" => (opt=-2.8284271247321175, sides=[3, 3, 1], nuniq=6),
    "CHSH_CS_TS_d2" => (opt=-3.999999999803662, sides=[3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 2, 2], nuniq=5),

    # =========================================================================
    # State Polynomial Formulation (pstateopt_first)
    # Objective: -ς(x₁y₁) - ς(x₁y₂) - ς(x₂y₁) + ς(x₂y₂)
    # =========================================================================
    "CHSH_State_Dense_d1" => (opt=-2.8284271247462325, sides=[11], nuniq=34),
    "CHSH_State_TS_d1" => (opt=-2.8284271247321175, sides=[3, 3, 1, 1, 1, 1, 1, 1, 1], nuniq=12),

    # =========================================================================
    # Trace Polynomial Formulation (ptraceopt_first)
    # Objective: -tr(x₁y₁) - tr(x₁y₂) - tr(x₂y₁) + tr(x₂y₂)
    # =========================================================================
    "CHSH_Trace_Dense_d1" => (opt=-2.8284271247462307, sides=[9], nuniq=21),
    "CHSH_Trace_TS_d1" => (opt=-2.8284271247321175, sides=[3, 3, 1, 1, 1, 1, 1], nuniq=10),
)
