# Generalized Rosenbrock Problem Definition
# ==========================================
# Problem: NC version of the generalized Rosenbrock function
# Variables: x[1:n]
# Objective: f = n + Σᵢ₌₂ⁿ [100xᵢ₋₁⁴ - 200xᵢ₋₁²xᵢ - 2xᵢ + 101xᵢ²]
# Constraints: none (unconstrained)
# Expected optimal: 1.0 at xᵢ = 0 for all i
#
# The problem has a chain-like sparsity structure: each term
# involves only consecutive variables (xᵢ₋₁, xᵢ).
# Degree: 4, so order d=2 is sufficient.

# Sparsity variants
const ROSENBROCK_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2, n=6),
    (name="CS", cs="MF", ts=false, order=2, n=6),
    (name="CS_TS", cs="MF", ts="MD", order=2, n=6),
]

# NCTSSOS setup
const ROSENBROCK_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:6]",
    objective_expr = """
        n = 6
        f = n * 1.0
        for i in 2:n
            f += 100 * x[i-1]^4 - 200 * x[i-1]^2 * x[i] - 2 * x[i] + 101 * x[i]^2
        end
    """,
    vars = "x",
    partition = 6,
    constraint = false,
)

# NCTSSoS setup
const ROSENBROCK_NCTSS = """
n = 6
registry, (x,) = create_noncommutative_variables([("x", 1:n)])
f = Float64(n) * one(typeof(x[1]))
for i in 2:n
    f = f + 100.0 * x[i-1]^4 - 200.0 * x[i-1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
end
pop = polyopt(f, registry)
"""
