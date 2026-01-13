# NCTSSOS Oracle Script: Bilocal Network (Example 8.1.1)
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_bilocal.jl
#
# Problem: Bilocal network quantum bound from Example 8.1.1
# Variables: x[1:2], y[1:2], z[1:2] (6 vars, unipotent U²=I)
# The bilocal scenario has three parties:
#   Alice (x), Bob (y), Charlie (z)
# with two independent sources connecting Alice-Bob and Bob-Charlie.
#
# Expected: Quantum bound of -4.0

include("oracle_utils.jl")

# Sparsity variants (CS not supported for state polynomials)
const BILOCAL_VARIANTS = [
    (name="Dense", ts=false, order=3),
    (name="TS", ts="MD", order=3),
]

print_header("Bilocal Network (Example 8.1.1)")

println("# Variables: x[1:2]={1,2}, y[1:2]={3,4}, z[1:2]={5,6}")
println("# Constraint: unipotent (U²=I)")
println("# Expected: -4.0")
println()

# Support structure from Example 8.1.1
# Each element is either:
#   - [[a;b;c]] → single expectation ς(xyz)
#   - [[a;b;c], [d;e;f]] → product of expectations ς(xyz) * ς(xyz')
# Variable mapping: x1=1, x2=2, y1=3, y2=4, z1=5, z2=6

# Build support for state polynomial optimization
# The pstateopt_first API uses vargroup to specify party structure
supp = [
    # Squared terms ς(·)ς(·) with same argument (8 terms)
    [[1, 3, 5], [1, 3, 5]], [[1, 3, 6], [1, 3, 6]], [[2, 3, 5], [2, 3, 5]],
    [[2, 3, 6], [2, 3, 6]], [[1, 4, 5], [1, 4, 5]], [[1, 4, 6], [1, 4, 6]],
    [[2, 4, 5], [2, 4, 5]], [[2, 4, 6], [2, 4, 6]],
    # Cross terms ς(·)ς(·') with different arguments (28 terms)
    [[1, 3, 5], [1, 3, 6]], [[1, 3, 5], [2, 3, 5]], [[1, 3, 5], [2, 3, 6]],
    [[1, 3, 5], [1, 4, 5]], [[1, 3, 5], [1, 4, 6]], [[1, 3, 5], [2, 4, 5]],
    [[1, 3, 5], [2, 4, 6]], [[1, 3, 6], [2, 3, 5]], [[1, 3, 6], [2, 3, 6]],
    [[1, 3, 6], [1, 4, 5]], [[1, 3, 6], [1, 4, 6]], [[1, 3, 6], [2, 4, 5]],
    [[1, 3, 6], [2, 4, 6]], [[2, 3, 5], [2, 3, 6]], [[2, 3, 5], [1, 4, 5]],
    [[2, 3, 5], [1, 4, 6]], [[2, 3, 5], [2, 4, 5]], [[2, 3, 5], [2, 4, 6]],
    [[2, 3, 6], [1, 4, 5]], [[2, 3, 6], [1, 4, 6]], [[2, 3, 6], [2, 4, 5]],
    [[2, 3, 6], [2, 4, 6]], [[1, 4, 5], [1, 4, 6]], [[1, 4, 5], [2, 4, 5]],
    [[1, 4, 5], [2, 4, 6]], [[1, 4, 6], [2, 4, 5]], [[1, 4, 6], [2, 4, 6]],
    [[2, 4, 5], [2, 4, 6]],
    # Linear terms ς(·) (8 terms)
    [[1, 3, 5]], [[1, 3, 6]], [[2, 3, 5]], [[2, 3, 6]],
    [[1, 4, 5]], [[1, 4, 6]], [[2, 4, 5]], [[2, 4, 6]]
]

# Coefficients from the paper
# Quadratic terms: 1/8 * [...] (36 terms)
coe_quad = (1 / 8) .* [
    1, 1, 1, 1, 1, 1, 1, 1,           # Squared terms (8)
    2, 2, 2, -2, 2, 2, -2,             # Cross terms row 1 (7)
    2, 2, -2, 2, 2, -2,                # Cross terms row 2 (6)
    2, -2, 2, 2, -2,                   # Cross terms row 3 (5)
    -2, 2, 2, -2,                      # Cross terms row 4 (4)
    -2, -2, 2, 2, -2, -2               # Cross terms row 5 (6)
]
# Linear terms: -[1, 1, 1, 1, 1, -1, -1, 1] (8 terms)
coe_linear = Float64[-1, -1, -1, -1, -1, 1, 1, -1]

# Combine coefficients
coe = vcat(coe_quad, coe_linear)

results = map(BILOCAL_VARIANTS) do v
    key = "Bilocal_8_1_1_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")

    # pstateopt_first(supp, coe, n, d; vargroup, TS, constraint)
    # n=6 total variables, vargroup=[2,2,2] for Alice, Bob, Charlie
    opt, data = pstateopt_first(supp, coe, 6, v.order;
        vargroup=[2, 2, 2],
        TS=v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

print_summary("BILOCAL", results)
