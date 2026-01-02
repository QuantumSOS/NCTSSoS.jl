# NCTSSOS Oracle Script: CHSH Bell Inequality
# =============================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_chsh.jl
#
# Output: Oracle values for CHSH problem with sparsity variants
#
# API:
#   CS=false → nctssos_first (dense or TS only)
#   CS="MF"  → cs_nctssos_first (CS or CS+TS)

using NCTSSOS, DynamicPolynomials
using MosekTools

# Variant definitions (mirrors CHSH_VARIANTS in chsh.jl)
const VARIANTS = [
    (name="Dense", cs=false, ts=false, order=1),
    (name="CS", cs="MF", ts=false, order=1),
    (name="TS", cs=false, ts="MD", order=1),
    (name="CS_TS", cs="MF", ts="MD", order=2),
]

println("=" ^ 70)
println("NCTSSOS Oracle: CHSH Bell Inequality")
println("=" ^ 70)
println()

# Problem setup
@ncpolyvar x[1:2] y[1:2]
f = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]

println("Problem: CHSH Bell inequality")
println("Variables: x[1:2] (Alice), y[1:2] (Bob)")
println("Objective: f = x₁y₁ + x₁y₂ + x₂y₁ - x₂y₂")
println("Constraint: unipotent (U²=I)")
println("Expected optimal: 2√2 ≈ 2.8284 (Tsirelson bound)")
println("Solving: max(-f) → opt ≈ -2.8284")
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
    key = "CHSH_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    if v.cs == false
        # Dense or TS only → nctssos_first
        opt, data = nctssos_first([-f], [x;y], v.order;
                                  TS=v.ts, partition=2, constraint="unipotent")
        extract_nctssos(key, opt, data)
    else
        # CS or CS+TS → cs_nctssos_first
        opt, data = cs_nctssos_first([-f], [x;y], v.order;
                                     TS=v.ts, CS=v.cs, partition=2, constraint="unipotent")
        extract_cs_nctssos(key, opt, data)
    end
    println()
end

# =============================================================================
# Summary
# =============================================================================
println("=" ^ 70)
println("CHSH ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const CHSH_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
