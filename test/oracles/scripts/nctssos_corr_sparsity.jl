# NCTSSOS Oracle Script: Correlative Sparsity Example
# =====================================================
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_corr_sparsity.jl
#
# Problem: NC polynomial with ball and linear constraints
# Variables: x[1:3]
# Objective: same as Example 1
# Inequality: 1 - xᵢ² ≥ 0 for i=1:3
# Inequality: xᵢ - 1/3 ≥ 0 for i=1:3

include("oracle_utils.jl")

# Sparsity variants
const CORR_SPARSITY_VARIANTS = [
    (name="CS", cs="MF", ts=false, order=3),
    (name="TS", cs=false, ts="MD", order=3),
]

# NCTSSOS parameters
const CORR_SPARSITY_NCTSSOS_PARAMS = (
    constraint = false,
    numeq = 0,
)

print_header("Correlative Sparsity Example")

# Problem setup
@ncpolyvar x[1:3]
obj = x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3*x[2]^2 - 2*x[1]*x[2]*x[1] +
    2*x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] + 6*x[3]^2 +
    9*x[2]^2*x[3] + 9*x[3]*x[2]^2 - 54*x[3]*x[2]*x[3] + 142*x[3]*x[2]^2*x[3]
ineqs = vcat([1 - x[i]^2 for i = 1:3], [x[i] - 1/3 for i = 1:3])
pop = vcat([obj], ineqs)
vars = x

println("Variables: x[1:3]")
println("Inequality: 1 - xᵢ² ≥ 0, xᵢ - 1/3 ≥ 0 for i=1:3")
println()

results = map(CORR_SPARSITY_VARIANTS) do v
    key = "CorrSparsity_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    if v.cs == false
        opt, data = nctssos_first(pop, vars, v.order; TS=v.ts, numeq=0)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first(pop, vars, v.order; TS=v.ts, CS=v.cs, numeq=0)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

print_summary("CORR_SPARSITY", results)
