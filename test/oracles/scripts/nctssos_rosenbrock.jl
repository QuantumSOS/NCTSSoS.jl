# NCTSSOS Oracle Script: Generalized Rosenbrock
# ==============================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_rosenbrock.jl
#
# Output: Oracle values for Rosenbrock problem with sparsity variants

using NCTSSOS, DynamicPolynomials
using MosekTools

# Variant definitions (mirrors ROSENBROCK_VARIANTS in rosenbrock.jl)
const VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2, n=6),
    (name="CS", cs="MF", ts=false, order=2, n=6),
    (name="CS_TS", cs="MF", ts="MD", order=2, n=6),
]

println("=" ^ 70)
println("NCTSSOS Oracle: Generalized Rosenbrock")
println("=" ^ 70)
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

# Run each variant
for v in VARIANTS
    n = v.n
    @ncpolyvar x[1:n]
    
    # Build Rosenbrock objective
    # f = n + Σᵢ₌₂ⁿ [100xᵢ₋₁⁴ - 200xᵢ₋₁²xᵢ - 2xᵢ + 101xᵢ²]
    f = n * 1.0
    for i in 2:n
        f += 100 * x[i-1]^4 - 200 * x[i-1]^2 * x[i] - 2 * x[i] + 101 * x[i]^2
    end
    
    key = "Rosenbrock_$(v.name)_n$(v.n)_d$(v.order)"
    println("# $(v.name) (n=$(v.n), order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    
    if v.cs == false
        opt, data = nctssos_first([f], x, v.order; TS=v.ts)
        extract_nctssos(key, opt, data)
    else
        opt, data = cs_nctssos_first([f], x, v.order; TS=v.ts, CS=v.cs)
        extract_cs_nctssos(key, opt, data)
    end
    println()
end

# =============================================================================
# Summary
# =============================================================================
println("=" ^ 70)
println("ROSENBROCK ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const ROSENBROCK_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
