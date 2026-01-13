# NCTSSOS Oracle Script: Example 2
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_example2.jl
#
# Problem: Constrained NC polynomial optimization
# Variables: x[1:2]
# Objective: f = 2 - x₁² + x₁x₂²x₁ - x₂²
# Inequality: g = 4 - x₁² - x₂² ≥ 0

include("oracle_utils.jl")

# Sparsity variants
const EXAMPLE2_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="CS_TS", cs="MF", ts="MD", order=2),
]

# NCTSSOS parameters
const EXAMPLE2_NCTSSOS_PARAMS = (
    constraint = false,
    numeq = 0,
)

print_header("Example 2 - 2-variable NC Polynomial with Constraints")

# Problem setup
@ncpolyvar x[1:2]
obj = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
g = 4 - x[1]^2 - x[2]^2  # inequality constraint
vars = x
pop = [obj; g]

println("Variables: x[1:2]")
println("Inequality: g = 4 - x₁² - x₂² ≥ 0")
println()

results = map(EXAMPLE2_VARIANTS) do v
    key = "Example2_$(v.name)_d$(v.order)"
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

print_summary("EXAMPLE2", results)
