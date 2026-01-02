# NCTSSOS Oracle Script: State Polynomial Optimization
# =====================================================
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_state_poly.jl
#
# Problems from Section 7.2 of the NCTSSOS paper.
# State polynomials use the ς() operator to represent quantum state expectations.
# API: pstateopt_first for pure state polynomial optimization
#      Uses coefficient/support format for state expectations

include("oracle_utils.jl")

# =============================================================================
# 7.2.0: CHSH State Polynomial (basic, linear in expectations)
# =============================================================================
# Variables: x[1:2], y[1:2] (unipotent, U²=I)
# Objective: sp = -ς(x₁y₁) - ς(x₁y₂) - ς(x₂y₁) + ς(x₂y₂)
# Expected: -2√2 ≈ -2.8284

const STATE_7_2_0_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=1),
    (name="TS", cs=false, ts="MD", order=1),
]

const STATE_7_2_0_NCTSSOS_PARAMS = (
    n = 4,
    vargroup = [2, 4],
    constraint = "unipotent",
)

# =============================================================================
# 7.2.1: Squared Expectations (requires order=3 for tight bound)
# =============================================================================
# Variables: x[1:2], y[1:2] (unipotent)
# Objective: sp = -(ς(x₁y₂) + ς(x₂y₁))² - (ς(x₁y₁) - ς(x₂y₂))²
# Expected: -4.0 (at order=3)
# NOTE: Requires mixword API for product terms - not supported by pstateopt_first

const STATE_7_2_1_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=3),
]

# =============================================================================
# 7.2.2: Covariance Expression
# =============================================================================
# Variables: x[1:3], y[1:3] (unipotent)
# Objective: sp = Σᵢⱼ cᵢⱼ * cov(i,j) where cov(a,b) = ς(xₐyᵦ) - ς(xₐ)ς(yᵦ)
# Expected: -5.0
# NOTE: Requires mixword API for product terms

const STATE_7_2_2_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
]

# =============================================================================
# 7.2.3: Mixed State Polynomial (squared expectations + products)
# =============================================================================
# Variables: x[1:2], y[1:2] (unipotent)
# Objective: Complex mix of linear, product, and squared expectations
# Expected: -3.5114802
# NOTE: Requires mixword API for product terms

const STATE_7_2_3_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="TS", cs=false, ts="MD", order=2),
]

print_header("State Polynomial Optimization")

# =============================================================================
# 7.2.0: CHSH State Polynomial
# =============================================================================
println("# 7.2.0: CHSH State Polynomial")
println("# Variables: x[1:2], y[1:2] (x=1,2, y=3,4)")
println("# Expected: -2√2 ≈ -2.8284")
println()

# sp = -ς(x₁y₁) - ς(x₁y₂) - ς(x₂y₁) + ς(x₂y₂)
coe_7_2_0 = [-1.0, -1.0, -1.0, 1.0]
supp_7_2_0 = [[[1,3]], [[1,4]], [[2,3]], [[2,4]]]

results_7_2_0 = map(STATE_7_2_0_VARIANTS) do v
    key = "State_7_2_0_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    ts_param = v.ts == false ? false : v.ts
    opt, data = pstateopt_first(supp_7_2_0, coe_7_2_0,
        STATE_7_2_0_NCTSSOS_PARAMS.n, v.order;
        vargroup=STATE_7_2_0_NCTSSOS_PARAMS.vargroup,
        TS=ts_param,
        constraint=STATE_7_2_0_NCTSSOS_PARAMS.constraint)
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# =============================================================================
# 7.2.1, 7.2.2, 7.2.3: Product terms (require mixword API)
# =============================================================================
println("\n# Note: Problems 7.2.1, 7.2.2, 7.2.3 require mixword API for product terms")
println("# Only pure state expectation problems (7.2.0 CHSH) can use pstateopt_first")
println()

print_summary("STATE_POLY", results_7_2_0)
