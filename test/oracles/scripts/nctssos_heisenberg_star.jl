# NCTSSOS Oracle Script: Heisenberg Star Graph
# =============================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_heisenberg_star.jl
#
# Output: Oracle values for Heisenberg star graph problem
#
# API:
#   CS=false → nctssos_first (dense or TS only)
#   CS="MF"  → cs_nctssos_first (CS or CS+TS)

using NCTSSOS, DynamicPolynomials
using MosekTools

# Variant definitions
const VARIANTS = [
    (name="Dense", cs=false, ts=false, order=1, num_sites=10),
    (name="CS", cs="MF", ts=false, order=1, num_sites=10),
    (name="Dense_n8", cs=false, ts=false, order=1, num_sites=8),
]

println("=" ^ 70)
println("NCTSSOS Oracle: Heisenberg Star Graph")
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
    num_sites = v.num_sites
    vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]
    n_vars = length(vec_idx2ij)
    
    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)
    
    # Define variables
    @ncpolyvar p[1:n_vars]
    
    # Objective: sum over star graph edges (center=1)
    obj = sum(p[findvaridx(1, k)] for k in 2:num_sites)
    
    # Triangle consistency constraints (equality constraints)
    gs = unique!([
        (
            p[findvaridx(sort([i, j])...)] * p[findvaridx(sort([j, k])...)] +
            p[findvaridx(sort([j, k])...)] * p[findvaridx(sort([i, j])...)] -
            p[findvaridx(sort([i, j])...)] - p[findvaridx(sort([j, k])...)] -
            p[findvaridx(sort([i, k])...)] + 1.0
        ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
        (i != j && j != k && i != k)
    ])
    
    pop = vcat([obj], gs)
    numeq = length(gs)
    
    key = "HeisenbergStar_$(v.name)_n$(v.num_sites)_d$(v.order)"
    println("# $(v.name) (n=$(v.num_sites), order=$(v.order), CS=$(v.cs), TS=$(v.ts))")
    println("# Variables: $n_vars, Equality constraints: $numeq")
    
    if v.cs == false
        opt, data = nctssos_first(pop, p, v.order;
                                  TS=v.ts, numeq=numeq, constraint="unipotent")
        extract_nctssos(key, opt, data)
    else
        opt, data = cs_nctssos_first(pop, p, v.order;
                                     TS=v.ts, CS=v.cs, numeq=numeq, constraint="unipotent")
        extract_cs_nctssos(key, opt, data)
    end
    println()
end

# =============================================================================
# Summary
# =============================================================================
println("=" ^ 70)
println("HEISENBERG_STAR ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const HEISENBERG_STAR_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
