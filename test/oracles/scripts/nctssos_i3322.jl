# NCTSSOS Oracle Script: I_3322 Bell Inequality
# ===============================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_i3322.jl
#
# Output: Oracle values for I_3322 problem with sparsity variants
#
# API:
#   CS=false → nctssos_first (dense or TS only)
#   CS="MF"  → cs_nctssos_first (CS or CS+TS)

using NCTSSOS, DynamicPolynomials
using MosekTools

# Variant definitions (mirrors I3322_VARIANTS in i3322.jl)
const VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2),
    (name="CS", cs="MF", ts=false, order=2),
    (name="TS", cs=false, ts="MD", order=3),
    (name="CS_TS", cs="MF", ts="MD", order=3),
]

println("=" ^ 70)
println("NCTSSOS Oracle: I_3322 Bell Inequality")
println("=" ^ 70)
println()

# Problem setup
@ncpolyvar x[1:3] y[1:3]
f = x[1]*(y[1]+y[2]+y[3]) + x[2]*(y[1]+y[2]-y[3]) + x[3]*(y[1]-y[2]) - x[1] - 2.0*y[1] - y[2]

println("Problem: I_3322 Bell inequality")
println("Variables: x[1:3] (Alice), y[1:3] (Bob)")
println("Objective: f = x₁(y₁+y₂+y₃) + x₂(y₁+y₂-y₃) + x₃(y₁-y₂) - x₁ - 2y₁ - y₂")
println("Constraint: projector (P²=P)")
println("Solving: max(-f)")
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
    key = "I3322_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    if v.cs == false
        # Dense or TS only → nctssos_first
        opt, data = nctssos_first([-f], [x;y], v.order;
                                  TS=v.ts, partition=3, constraint="projector")
        extract_nctssos(key, opt, data)
    else
        # CS or CS+TS → cs_nctssos_first
        opt, data = cs_nctssos_first([-f], [x;y], v.order;
                                     TS=v.ts, CS=v.cs, partition=3, constraint="projector")
        extract_cs_nctssos(key, opt, data)
    end
    println()
end

# =============================================================================
# Summary
# =============================================================================
println("=" ^ 70)
println("I3322 ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const I3322_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
