# Correlative Sparsity Example
# =============================
# Problem: NC polynomial with ball and linear constraints
# Variables: x[1:3]
# Objective: same as Example 1
# Inequality: 1 - xᵢ² ≥ 0 for i=1:3
# Inequality: xᵢ - 1/3 ≥ 0 for i=1:3

# Sparsity variants: order specifies relaxation order for oracle tests
const CORR_SPARSITY_VARIANTS = [
    (name="CS", cs="MF", ts=false, order=3),
    (name="TS", cs=false, ts="MD", order=3),
]

# NCTSSOS setup
const CORR_SPARSITY_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:3]",
    objective_expr = """
        x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3*x[2]^2 - 2*x[1]*x[2]*x[1] +
        2*x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] + 6*x[3]^2 +
        9*x[2]^2*x[3] + 9*x[3]*x[2]^2 - 54*x[3]*x[2]*x[3] + 142*x[3]*x[2]^2*x[3]
    """,
    vars = "x",
    partition = 3,
    constraint = false,
)

# NCTSSoS setup
const CORR_SPARSITY_NCTSS = """
reg, (x,) = create_noncommutative_variables([("x", 1:3)])
f = 1.0*x[1]^2 - 1.0*x[1]*x[2] - 1.0*x[2]*x[1] + 3.0*x[2]^2 - 2.0*x[1]*x[2]*x[1] +
    2.0*x[1]*x[2]^2*x[1] - 1.0*x[2]*x[3] - 1.0*x[3]*x[2] + 6.0*x[3]^2 +
    9.0*x[2]^2*x[3] + 9.0*x[3]*x[2]^2 - 54.0*x[3]*x[2]*x[3] + 142.0*x[3]*x[2]^2*x[3]
cons = vcat([1.0 - 1.0*x[i]^2 for i = 1:3], [1.0*x[i] - 1.0/3 for i = 1:3])
pop = polyopt(f, reg; ineq_constraints=cons)
"""
