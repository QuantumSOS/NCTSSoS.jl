# NCTSSOS Oracle Script: Generalized Rosenbrock
# ==============================================
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_rosenbrock.jl
#
# Problem: NC version of the generalized Rosenbrock function
# Variables: x[1:n]
# Objective: f = n + Σᵢ₌₂ⁿ [100xᵢ₋₁⁴ - 200xᵢ₋₁²xᵢ - 2xᵢ + 101xᵢ²]
# Constraints: none (unconstrained)
# Expected optimal: 1.0 at xᵢ = 0 for all i
#
# The problem has a chain-like sparsity structure: each term
# involves only consecutive variables (xᵢ₋₁, xᵢ).
# Degree: 4, so order d=2 is sufficient.

include("oracle_utils.jl")

# Sparsity variants
const ROSENBROCK_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2, n=6),
    (name="CS", cs="MF", ts=false, order=2, n=6),
    (name="CS_TS", cs="MF", ts="MD", order=2, n=6),
]

# NCTSSOS parameters
const ROSENBROCK_NCTSSOS_PARAMS = (
    constraint = false,
)

print_header("Generalized Rosenbrock")

println("Constraint: none (unconstrained NC)")
println("Expected: 1.0")
println()

results = map(ROSENBROCK_VARIANTS) do v
    n = v.n
    @ncpolyvar x[1:n]
    
    # Build Rosenbrock objective
    # f = n + Σᵢ₌₂ⁿ [100xᵢ₋₁⁴ - 200xᵢ₋₁²xᵢ - 2xᵢ + 101xᵢ²]
    obj = n * 1.0
    for i in 2:n
        obj += 100 * x[i-1]^4 - 200 * x[i-1]^2 * x[i] - 2 * x[i] + 101 * x[i]^2
    end
    
    key = "Rosenbrock_$(v.name)_n$(v.n)_d$(v.order)"
    println("# $(v.name) (n=$(v.n), order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    if v.cs == false
        opt, data = nctssos_first([obj], x, v.order; TS=v.ts)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first([obj], x, v.order; TS=v.ts, CS=v.cs)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

print_summary("ROSENBROCK", results)
