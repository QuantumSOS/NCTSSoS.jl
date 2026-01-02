# NCTSSOS Oracle Script: I_3322 Bell Inequality
# ===============================================
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_i3322.jl
#
# Problem: I_3322 Bell inequality with projector constraint (P²=P)
# Variables: x[1:3] (Alice), y[1:3] (Bob)
# Objective: f = x₁(y₁+y₂+y₃) + x₂(y₁+y₂-y₃) + x₃(y₁-y₂) - x₁ - 2y₁ - y₂
# Expected optimal: ≈ -0.2509 (quantum bound)
#
# NOTE on CS_TS: Does NOT converge to -0.2509 even at high order.
# Combined sparsity loses cross-clique constraints (similar to CHSH).

include("oracle_utils.jl")

# Sparsity variants
const I3322_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="CS", cs="MF", ts=false, order=2),
    (name="TS", cs=false, ts="MD", order=3),
    (name="CS_TS", cs="MF", ts="MD", order=3),
    (name="CS_TS", cs="MF", ts="MD", order=4),
]

# NCTSSOS parameters
const I3322_NCTSSOS_PARAMS = (
    partition = 3,
    constraint = "projector",
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
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    if v.cs == false
        opt, data = nctssos_first([obj], vars, v.order;
            TS=v.ts, partition=I3322_NCTSSOS_PARAMS.partition,
            constraint=I3322_NCTSSOS_PARAMS.constraint)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first([obj], vars, v.order;
            TS=v.ts, CS=v.cs, partition=I3322_NCTSSOS_PARAMS.partition,
            constraint=I3322_NCTSSOS_PARAMS.constraint)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

print_summary("I3322", results)
