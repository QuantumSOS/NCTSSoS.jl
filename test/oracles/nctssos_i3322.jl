# NCTSSOS Oracle Script: I_3322 Bell Inequality
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_i3322.jl
#
# Problem: I_3322 Bell inequality with projection constraint (P²=P)
# Variables: x[1:3] (Alice), y[1:3] (Bob)
# Objective: f = x₁(y₁+y₂+y₃) + x₂(y₁+y₂-y₃) + x₃(y₁-y₂) - x₁ - 2y₁ - y₂
# Expected optimal: ≈ -0.2509 (quantum bound)

include("oracle_utils.jl")

# Sparsity variants - generates all cases needed by test file
# NOTE: For Bell inequalities with partition parameter, correlative sparsity (CS)
# is NOT applicable via cs_nctssos_first - partition already exploits bipartite structure.
# All variants use nctssos_first with partition parameter.
const I3322_VARIANTS = [
    # Order 2: Dense
    (name="Dense", ts=false, order=2),
    # Order 3: Dense and Term Sparsity
    (name="Dense", ts=false, order=3),
    (name="TS", ts="MD", order=3),
]

# NCTSSOS parameters
const I3322_NCTSSOS_PARAMS = (
    partition = 3,
    constraint = "projection",
)

print_header("I_3322 Bell Inequality")

# Problem setup
@ncpolyvar x[1:3] y[1:3]
obj = -(x[1]*(y[1]+y[2]+y[3]) + x[2]*(y[1]+y[2]-y[3]) + x[3]*(y[1]-y[2]) - x[1] - 2.0*y[1] - y[2])
vars = [x; y]

println("Objective: -(x₁(y₁+y₂+y₃) + x₂(y₁+y₂-y₃) + x₃(y₁-y₂) - x₁ - 2y₁ - y₂)")
println("Constraint: $(I3322_NCTSSOS_PARAMS.constraint)")
println("Expected: ≈ -0.2509")
println()

results = map(I3322_VARIANTS) do v
    key = "I3322_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")

    opt, data = nctssos_first([obj], vars, v.order;
        TS=v.ts, partition=I3322_NCTSSOS_PARAMS.partition,
        constraint=I3322_NCTSSOS_PARAMS.constraint)
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

print_summary("I3322", results)
