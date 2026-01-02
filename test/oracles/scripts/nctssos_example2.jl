# NCTSSOS Oracle Script: Example 2
# ==================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_example2.jl
#
# Output: Oracle values for Example 2 problem

using NCTSSOS, DynamicPolynomials
using MosekTools

# Variant definitions (mirrors EXAMPLE2_VARIANTS in example2.jl)
const VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="CS_TS", cs="MF", ts="MD", order=2),
]

println("=" ^ 70)
println("NCTSSOS Oracle: Example 2")
println("=" ^ 70)
println()

# Problem setup
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
g = 4 - x[1]^2 - x[2]^2  # inequality

println("Problem: Example 2 - 2-variable NC polynomial with constraints")
println("Variables: x[1:2]")
println("Inequality: g = 4 - x₁² - x₂² ≥ 0")
println()

# Store results
results = Dict{String, NamedTuple}()

# Helper to extract oracle info from nctssos_first (CS=false)
function extract_nctssos(name, opt, data)
    sides = [size(M, 1) for M in data.moment]
    nuniq = length(data.ksupp)
    results[name] = (opt=opt, sides=sides, nuniq=nuniq)
    println("\"$name\" => (opt=$opt, sides=$sides, nuniq=$nuniq),")
end

# Helper to extract oracle info from cs_nctssos_first (CS="MF")
function extract_cs_nctssos(name, opt, data)
    sides = [size(M, 1) for clique in data.moment for M in clique]
    nuniq = length(data.ksupp)
    results[name] = (opt=opt, sides=sides, nuniq=nuniq)
    println("\"$name\" => (opt=$opt, sides=$sides, nuniq=$nuniq),")
end

# Run each variant at its specified order
for v in VARIANTS
    key = "Example2_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    if v.cs == false
        opt, data = nctssos_first([f; g], x, v.order; TS=v.ts, numeq=0)
        extract_nctssos(key, opt, data)
    else
        opt, data = cs_nctssos_first([f; g], x, v.order; TS=v.ts, CS=v.cs, numeq=0)
        extract_cs_nctssos(key, opt, data)
    end
    println()
end

# =============================================================================
# Summary
# =============================================================================
println("=" ^ 70)
println("EXAMPLE2 ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const EXAMPLE2_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
