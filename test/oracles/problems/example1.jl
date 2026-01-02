# Example 1: 3-variable NC Polynomial
# ====================================
# Problem: Unconstrained NC polynomial optimization
# Variables: x[1:3]
# Objective: f = x₁² - x₁x₂ - x₂x₁ + 3x₂² - 2x₁x₂x₁ + 2x₁x₂²x₁ 
#              - x₂x₃ - x₃x₂ + 6x₃² + 9x₂²x₃ + 9x₃x₂² - 54x₃x₂x₃ + 142x₃x₂²x₃

# Sparsity variants: order specifies relaxation order for oracle tests
const EXAMPLE1_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="TS", cs=false, ts="MD", order=2),
]

# NCTSSOS setup
const EXAMPLE1_NCTSSOS = (
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
const EXAMPLE1_NCTSS = """
reg, (x,) = create_noncommutative_variables([("x", 1:3)])
f = 1.0*x[1]^2 - 1.0*x[1]*x[2] - 1.0*x[2]*x[1] + 3.0*x[2]^2 - 2.0*x[1]*x[2]*x[1] +
    2.0*x[1]*x[2]^2*x[1] - 1.0*x[2]*x[3] - 1.0*x[3]*x[2] + 6.0*x[3]^2 +
    9.0*x[2]^2*x[3] + 9.0*x[3]*x[2]^2 - 54.0*x[3]*x[2]*x[3] + 142.0*x[3]*x[2]^2*x[3]
pop = polyopt(f, reg)
"""
