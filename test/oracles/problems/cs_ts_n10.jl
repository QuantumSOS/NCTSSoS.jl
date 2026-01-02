# CS TS Example (n=10) Problem Definition
# ========================================
# Problem: Large NC polynomial with correlative and term sparsity structure
# Variables: x[1:10]
# Objective: Sum of local terms with overlapping support (banded structure)
# 
# Each term has support within i±5 range, creating natural sparsity.
# Used to validate combined CS+TS exploitation.

# Sparsity variants: order specifies relaxation order for oracle tests
const CS_TS_N10_VARIANTS = [
    (name="CS_TS", cs="MF", ts="MD", order=3),
]

# NCTSSOS setup
# Note: NCTSSOS minimizes by default via nctssos_first
const CS_TS_N10_NCTSSOS = (
    vars_expr = "@ncpolyvar x[1:10]",
    objective_expr = """
        n = 10
        f = 0 * x[1]
        for i = 1:n
            jset = max(1, i - 5):min(n, i + 1)
            jset = setdiff(jset, i)
            f += (2x[i] + 5x[i]^3 + 1)^2
            f -= sum([
                4x[i]*x[j] + 10x[i]^3*x[j] + 2x[j] +
                4x[i]*x[j]^2 + 10x[i]^3*x[j]^2 + 2x[j]^2 for j in jset
            ])
            f += sum([
                x[j]*x[k] + 2x[j]^2*x[k] + x[j]^2*x[k]^2 for j in jset for k in jset
            ])
        end
    """,
    vars = "x",
    partition = 10,
    constraint = false,
)

# NCTSSoS setup
const CS_TS_N10_NCTSS = """
n = 10
registry, (x,) = create_noncommutative_variables([("x", 1:n)])
f = 1.0 * x[1] - 1.0 * x[1]  # zero polynomial
for i = 1:n
    jset = max(1, i - 5):min(n, i + 1)
    jset = setdiff(jset, i)
    f += (2.0 * x[i] + 5.0 * x[i]^3 + 1)^2
    f -= sum([
        4.0 * x[i] * x[j] + 10.0 * x[i]^3 * x[j] + 2.0 * x[j] +
        4.0 * x[i] * x[j]^2 + 10.0 * x[i]^3 * x[j]^2 + 2.0 * x[j]^2 for j in jset
    ])
    f += sum([
        1.0 * x[j] * x[k] + 2.0 * x[j]^2 * x[k] + 1.0 * x[j]^2 * x[k]^2 for j in jset for k in jset
    ])
end
cons = vcat([(1.0 - x[i]^2) for i = 1:n], [(1.0 * x[i] - 1.0 / 3) for i = 1:n])
pop = polyopt(f, registry; ineq_constraints=cons)
"""
