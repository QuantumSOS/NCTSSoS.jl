# NCTSSOS Oracle Script: Example 1
# ==================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_example1.jl
#
# Output: Oracle values for Example 1 problem

using NCTSSOS, DynamicPolynomials
using MosekTools

# Variant definitions (mirrors EXAMPLE1_VARIANTS in example1.jl)
const VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="TS", cs=false, ts="MD", order=2),
]

println("=" ^ 70)
println("NCTSSOS Oracle: Example 1")
println("=" ^ 70)
println()

# Problem setup
@ncpolyvar x[1:3]
f = x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3*x[2]^2 - 2*x[1]*x[2]*x[1] +
    2*x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] + 6*x[3]^2 +
    9*x[2]^2*x[3] + 9*x[3]*x[2]^2 - 54*x[3]*x[2]*x[3] + 142*x[3]*x[2]^2*x[3]

println("Problem: Example 1 - 3-variable NC polynomial")
println("Variables: x[1:3]")
println("Constraint: none (unconstrained NC)")
println()

# Store results
results = Dict{String, NamedTuple}()

# Helper to extract oracle info
function extract_nctssos(name, opt, data)
    sides = [size(M, 1) for M in data.moment]
    nuniq = length(data.ksupp)
    results[name] = (opt=opt, sides=sides, nuniq=nuniq)
    println("\"$name\" => (opt=$opt, sides=$sides, nuniq=$nuniq),")
end

# Run each variant at its specified order
for v in VARIANTS
    key = "Example1_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    opt, data = nctssos_first([f], x, v.order; TS=v.ts)
    extract_nctssos(key, opt, data)
    println()
end

# =============================================================================
# Summary
# =============================================================================
println("=" ^ 70)
println("EXAMPLE1 ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const EXAMPLE1_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
