# Example 2: 2-variable NC Polynomial with Constraints
# =====================================================
# Problem: Constrained NC polynomial optimization
# Variables: x[1:2]
# Objective: f = 2 - x₁² + x₁x₂²x₁ - x₂²
# Inequality: g = 4 - x₁² - x₂² ≥ 0
# Equality: h = x₁x₂ + x₂x₁ - 2 = 0

# Sparsity variants: order specifies relaxation order for oracle tests
const EXAMPLE2_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="CS_TS", cs="MF", ts="MD", order=2),
]

# NCTSSOS setup
const EXAMPLE2_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:2]",
    objective_expr = "2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2",
    ineq_expr = "4 - x[1]^2 - x[2]^2",
    eq_expr = "x[1]*x[2] + x[2]*x[1] - 2",
    vars = "x",
    partition = 2,
    constraint = false,
)

# NCTSSoS setup
const EXAMPLE2_NCTSS = """
reg, (x,) = create_noncommutative_variables([("x", 1:2)])
f = 2.0 - 1.0*x[1]^2 + 1.0*x[1]*x[2]^2*x[1] - 1.0*x[2]^2
g = 4.0 - 1.0*x[1]^2 - 1.0*x[2]^2
h = 1.0*x[1]*x[2] + 1.0*x[2]*x[1] - 2.0
pop = polyopt(f, reg; eq_constraints=[h], ineq_constraints=[g])
"""
