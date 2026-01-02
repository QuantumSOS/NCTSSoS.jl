# NCTSSOS Oracle Script: Correlative Sparsity Example
# =====================================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_corr_sparsity.jl
#
# Output: Oracle values for correlative sparsity problem

using NCTSSOS, DynamicPolynomials
using MosekTools

# Variant definitions (mirrors CORR_SPARSITY_VARIANTS in corr_sparsity.jl)
const VARIANTS = [
    (name="CS", cs="MF", ts=false, order=3),
    (name="TS", cs=false, ts="MD", order=3),
]

println("=" ^ 70)
println("NCTSSOS Oracle: Correlative Sparsity Example")
println("=" ^ 70)
println()

# Problem setup
@ncpolyvar x[1:3]
f = x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3*x[2]^2 - 2*x[1]*x[2]*x[1] +
    2*x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] + 6*x[3]^2 +
    9*x[2]^2*x[3] + 9*x[3]*x[2]^2 - 54*x[3]*x[2]*x[3] + 142*x[3]*x[2]^2*x[3]
ineqs = vcat([1 - x[i]^2 for i = 1:3], [x[i] - 1/3 for i = 1:3])
pop = vcat([f], ineqs)

println("Problem: Correlative Sparsity - 3-variable NC polynomial")
println("Variables: x[1:3]")
println("Inequality: 1 - xᵢ² ≥ 0, xᵢ - 1/3 ≥ 0 for i=1:3")
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
    key = "CorrSparsity_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    if v.cs == false
        opt, data = nctssos_first(pop, x, v.order; TS=v.ts, numeq=0)
        extract_nctssos(key, opt, data)
    else
        opt, data = cs_nctssos_first(pop, x, v.order; TS=v.ts, CS=v.cs, numeq=0)
        extract_cs_nctssos(key, opt, data)
    end
    println()
end

# =============================================================================
# Summary
# =============================================================================
println("=" ^ 70)
println("CORR_SPARSITY ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const CORR_SPARSITY_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
