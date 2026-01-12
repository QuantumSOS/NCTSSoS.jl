# NCTSSOS Oracle Script: Example 1
# ==================================
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_example1.jl
#
# Problem: Unconstrained NC polynomial optimization
# Variables: x[1:3]
# Objective: f = x₁² - x₁x₂ - x₂x₁ + 3x₂² - 2x₁x₂x₁ + 2x₁x₂²x₁
#              - x₂x₃ - x₃x₂ + 6x₃² + 9x₂²x₃ + 9x₃x₂² - 54x₃x₂x₃ + 142x₃x₂²x₃

include("oracle_utils.jl")

# Sparsity variants
const EXAMPLE1_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="TS", cs=false, ts="MD", order=2),
]

# NCTSSOS parameters (no constraint for unconstrained)
const EXAMPLE1_NCTSSOS_PARAMS = (
    constraint = false,
)

print_header("Example 1 - 3-variable NC Polynomial")

# Problem setup
@ncpolyvar x[1:3]
obj = x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3*x[2]^2 - 2*x[1]*x[2]*x[1] +
    2*x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] + 6*x[3]^2 +
    9*x[2]^2*x[3] + 9*x[3]*x[2]^2 - 54*x[3]*x[2]*x[3] + 142*x[3]*x[2]^2*x[3]
vars = x

println("Variables: x[1:3]")
println("Constraint: none (unconstrained NC)")
println()

results = map(EXAMPLE1_VARIANTS) do v
    key = "Example1_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")

    opt, data = nctssos_first(obj, vars; order=v.order, TS=v.ts)
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

print_summary("EXAMPLE1", results)
