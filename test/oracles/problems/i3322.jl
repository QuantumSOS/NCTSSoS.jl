# I_3322 Bell Inequality Problem Definition
# ==========================================
# Problem: I_3322 Bell inequality with projector constraint (P²=P)
# Variables: x[1:3] (Alice), y[1:3] (Bob)
# Objective: f = x₁(y₁+y₂+y₃) + x₂(y₁+y₂-y₃) + x₃(y₁-y₂) - x₁ - 2y₁ - y₂

# Sparsity variants: order matches test cases in test/solvers/
const I3322_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="CS", cs="MF", ts=false, order=2),
    (name="TS", cs=false, ts="MD", order=3),
    (name="CS_TS", cs="MF", ts="MD", order=3),
]

# NCTSSOS setup
const I3322_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:3] y[1:3]",
    objective_expr = "-(x[1]*(y[1]+y[2]+y[3]) + x[2]*(y[1]+y[2]-y[3]) + x[3]*(y[1]-y[2]) - x[1] - 2.0*y[1] - y[2])",
    vars = "[x;y]",
    partition = 3,
    constraint = "projector",
)

# NCTSSoS setup
const I3322_NCTSS = """
reg, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
f = x[1]*(y[1]+y[2]+y[3]) + x[2]*(y[1]+y[2]-y[3]) + x[3]*(y[1]-y[2]) - x[1] - 2.0*y[1] - y[2]
pop = polyopt(-f, reg)
"""
