# NCTSSOS Oracle Script: CS TS Example (n=10)
# =============================================
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_cs_ts_n10.jl
#
# Problem: Large NC polynomial with correlative and term sparsity structure
# Variables: x[1:10]
# Objective: Sum of local terms with overlapping support (banded structure)
#
# Each term has support within i±5 range, creating natural sparsity.
# Used to validate combined CS+TS exploitation.

include("oracle_utils.jl")

# Sparsity variants
const CS_TS_N10_VARIANTS = [
    (name="CS_TS", cs="MF", ts="MD", order=3),
]

# NCTSSOS parameters
const CS_TS_N10_NCTSSOS_PARAMS = (
    constraint = false,
    numeq = 0,
)

print_header("CS TS Example (n=10)")

# Problem setup
n = 10
@ncpolyvar x[1:n]

# Build objective
obj = 0 * x[1]
for i = 1:n
    global obj
    jset = max(1, i - 5):min(n, i + 1)
    jset = setdiff(jset, i)
    obj += (2x[i] + 5x[i]^3 + 1)^2
    obj -= sum([
        4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] +
        4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset
    ])
    obj += sum([
        x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
    ])
end

# Constraints: ball and linear
cons = vcat([1 - x[i]^2 for i = 1:n], [x[i] - 1 / 3 for i = 1:n])
pop = vcat([obj], cons)
vars = x

println("Variables: x[1:$n]")
println("Objective: sum of local terms with banded support")
println("Inequality: 1 - xᵢ² ≥ 0, xᵢ - 1/3 ≥ 0 for i=1:$n")
println()

results = map(CS_TS_N10_VARIANTS) do v
    key = "CS_TS_N10_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    opt, data = cs_nctssos_first(pop, vars, v.order; TS=v.ts, CS=v.cs, numeq=0)
    result = extract_oracle(key, opt, data; use_cs=true)
    print_oracle(result)
    println()
    result
end

print_summary("CS_TS_N10", results)
