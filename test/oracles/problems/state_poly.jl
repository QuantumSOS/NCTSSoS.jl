# State Polynomial Optimization Problem Definitions
# ==================================================
# Problems from Section 7.2 of the NCTSSOS paper.
# State polynomials use the ς() operator to represent quantum state expectations.
#
# Sign convention:
#   NCTSSOS maximizes → solve max(-f)
#   NCTSSoS minimizes → solve min(-f)
#
# Note: NCTSSOS uses `stateopt_first` for state polynomial optimization,
# which handles the state expectation structure automatically.

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

const STATE_7_2_0_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:2] y[1:2]",
    # NCTSSOS stateopt format: coe and supp arrays
    coe_expr = "[-1, -1, -1, 1]",
    supp_expr = "[[[1;3]], [[1;4]], [[2;3]], [[2;4]]]",
    vars = "[x;y]",
    partition = 2,
    constraint = "unipotent",
)

const STATE_7_2_0_NCTSS = """
reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
sp = -ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
spop = polyopt(sp * one(typeof(x[1])), reg)
"""

# =============================================================================
# 7.2.1: Squared Expectations (requires order=3 for tight bound)
# =============================================================================
# Variables: x[1:2], y[1:2] (unipotent)
# Objective: sp = -(ς(x₁y₂) + ς(x₂y₁))² - (ς(x₁y₁) - ς(x₂y₂))²
# Expected: -4.0 (at order=3)

const STATE_7_2_1_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=3),
]

const STATE_7_2_1_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:2] y[1:2]",
    # State polynomial with squared expectations
    coe_expr = "[-1, -1, -2, -1, -1, 2]",
    supp_expr = "[[[1;4], [1;4]], [[2;3], [2;3]], [[1;4], [2;3]], [[1;3], [1;3]], [[2;4], [2;4]], [[1;3], [2;4]]]",
    vars = "[x;y]",
    partition = 2,
    constraint = "unipotent",
)

const STATE_7_2_1_NCTSS = """
reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2
spop = polyopt(sp * one(typeof(x[1])), reg)
"""

# =============================================================================
# 7.2.2: Covariance Expression
# =============================================================================
# Variables: x[1:3], y[1:3] (unipotent)
# Objective: sp = Σᵢⱼ cᵢⱼ * cov(i,j) where cov(a,b) = ς(xₐyᵦ) - ς(xₐ)ς(yᵦ)
# Expected: -5.0

const STATE_7_2_2_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
]

const STATE_7_2_2_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:3] y[1:3]",
    # Covariance terms: cov(i,j) = <xy> - <x><y>
    coe_expr = "[1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1]",
    supp_expr = """
        [[[1;4]], [[1], [4]], [[1;5]], [[1], [5]], [[1;6]], [[1], [6]],
         [[2;4]], [[2], [4]], [[2;5]], [[2], [5]], [[2;6]], [[2], [6]],
         [[3;4]], [[3], [4]], [[3;5]], [[3], [5]]]
    """,
    vars = "[x;y]",
    partition = 3,
    constraint = "unipotent",
)

const STATE_7_2_2_NCTSS = """
reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
cov(a, b) = 1.0 * ς(x[a] * y[b]) - 1.0 * ς(x[a]) * ς(y[b])
sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)
spop = polyopt(sp * one(typeof(x[1])), reg)
"""

# =============================================================================
# 7.2.3: Mixed State Polynomial (squared expectations + products)
# =============================================================================
# Variables: x[1:2], y[1:2] (unipotent)
# Objective: Complex mix of linear, product, and squared expectations
# Expected: -3.5114802

const STATE_7_2_3_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="TS", cs=false, ts="MD", order=2),
]

const STATE_7_2_3_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:2] y[1:2]",
    # From NCTSSOS: coe = -[1; 1; 1; -1; 1; 1; 1; -1; -1; -1; -1; -1]
    coe_expr = "[-1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1]",
    supp_expr = """
        [[[2]], [[3]], [[4]], [[1;3]], [[2;3]], [[1;4]], [[2;4]], 
         [[1], [3]], [[2], [3]], [[2], [4]], [[1], [1]], [[4], [4]]]
    """,
    vars = "[x;y]",
    partition = 2,
    constraint = "unipotent",
)

const STATE_7_2_3_NCTSS = """
reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
sp = -1.0 * ς(x[2]) - 1.0 * ς(y[1]) - 1.0 * ς(y[2]) +
     1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[2] * y[1]) -
     1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[2]) +
     1.0 * ς(x[1]) * ς(y[1]) + 1.0 * ς(x[2]) * ς(y[1]) +
     1.0 * ς(x[2]) * ς(y[2]) +
     1.0 * ς(x[1]) * ς(x[1]) + 1.0 * ς(y[2]) * ς(y[2])
spop = polyopt(sp * one(typeof(x[1])), reg)
"""
