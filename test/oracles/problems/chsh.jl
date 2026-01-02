# CHSH Bell Inequality Problem Definition
# =========================================
# Problem: CHSH Bell inequality with unipotent constraint (U²=I)
# Variables: x[1:2] (Alice), y[1:2] (Bob)
# Objective: f = x₁y₁ + x₁y₂ + x₂y₁ - x₂y₂
# Expected optimal: 2√2 ≈ 2.8284 (Tsirelson bound)
#
# Sign convention:
#   NCTSSOS maximizes → solve max(-f) → opt ≈ -2.8284
#   NCTSSoS minimizes → solve min(-f) → opt ≈ -2.8284

# Sparsity variants: order specifies relaxation order for oracle tests
const CHSH_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=1),
    (name="CS", cs="MF", ts=false, order=1),
    (name="TS", cs=false, ts="MD", order=1),
    (name="CS_TS", cs="MF", ts="MD", order=2),
]

# NCTSSOS setup
const CHSH_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:2] y[1:2]",
    objective_expr = "-(x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2])",
    vars = "[x;y]",
    partition = 2,
    constraint = "unipotent",
)

# NCTSSoS setup
const CHSH_NCTSS = """
reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
f = x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
pop = polyopt(-f, reg)
"""
