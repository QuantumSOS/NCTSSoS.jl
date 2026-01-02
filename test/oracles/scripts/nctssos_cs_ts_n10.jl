# NCTSSOS Oracle Script: CS TS Example (n=10)
# =============================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_cs_ts_n10.jl
#
# Output: Oracle values for CS/TS n=10 problem
#
# API:
#   cs_nctssos_first for combined CS+TS

using NCTSSOS, DynamicPolynomials
using MosekTools

# Variant definitions (mirrors CS_TS_N10_VARIANTS in cs_ts_n10.jl)
const VARIANTS = [
    (name="CS_TS", cs="MF", ts="MD", order=3),
]

println("=" ^ 70)
println("NCTSSOS Oracle: CS TS Example (n=10)")
println("=" ^ 70)
println()

# Problem setup
n = 10
@ncpolyvar x[1:n]

# Build objective
f = 0 * x[1]
for i = 1:n
    jset = max(1, i - 5):min(n, i + 1)
    jset = setdiff(jset, i)
    f += (2x[i] + 5x[i]^3 + 1)^2
    f -= sum([
        4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] +
        4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset
    ])
    f += sum([
        x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
    ])
end

# Constraints: ball and linear
cons = vcat([1 - x[i]^2 for i = 1:n], [x[i] - 1 / 3 for i = 1:n])
pop = vcat([f], cons)

println("Problem: CS/TS Example (n=10)")
println("Variables: x[1:$n]")
println("Objective: sum of local terms with banded support")
println("Inequality: 1 - xᵢ² ≥ 0, xᵢ - 1/3 ≥ 0 for i=1:$n")
println()

# Store results
results = Dict{String, NamedTuple}()

# Helper to extract oracle info from cs_nctssos_first (CS="MF")
function extract_cs_nctssos(name, opt, data)
    sides = [size(M, 1) for clique in data.moment for M in clique]
    nuniq = length(data.ksupp)
    results[name] = (opt=opt, sides=sides, nuniq=nuniq)
    println("\"$name\" => (opt=$opt, sides=$sides, nuniq=$nuniq),")
end

# Run each variant
for v in VARIANTS
    key = "CS_TS_N10_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    opt, data = cs_nctssos_first(pop, x, v.order; TS=v.ts, CS=v.cs, numeq=0)
    extract_cs_nctssos(key, opt, data)
    println()
end

# =============================================================================
# Summary
# =============================================================================
println("=" ^ 70)
println("CS_TS_N10 ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const CS_TS_N10_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
