# NCTSSOS Oracle Script: Heisenberg Star Graph
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_heisenberg_star.jl
#
# Problem: Heisenberg model on star graph with unipotent constraint (U²=I)
# Variables: pij for each edge (i,j) in the star graph
# Objective: sum of pij over star edges (center=1)
# Constraints: Jordan-Wigner type triangle consistency relations
#
# The star graph has n sites with site 1 at center.
# Edges: (1,2), (1,3), ..., (1,n)

include("oracle_utils.jl")

# Sparsity variants
const HEISENBERG_STAR_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=1, num_sites=10),
    (name="CS", cs="MF", ts=false, order=1, num_sites=10),
    (name="Dense_n8", cs=false, ts=false, order=1, num_sites=8),
]

# NCTSSOS parameters
const HEISENBERG_STAR_NCTSSOS_PARAMS = (
    constraint = "unipotent",
)

print_header("Heisenberg Star Graph")

println("Constraint: $(HEISENBERG_STAR_NCTSSOS_PARAMS.constraint)")
println()

results = map(HEISENBERG_STAR_VARIANTS) do v
    num_sites = v.num_sites
    vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]
    n_vars = length(vec_idx2ij)
    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)
    
    # Define variables
    @ncpolyvar p[1:n_vars]
    
    # Objective: sum over star graph edges (center=1)
    obj = sum(p[findvaridx(1, k)] for k in 2:num_sites)
    
    # Triangle consistency constraints (equality)
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
            TS=v.ts, numeq=numeq,
            constraint=HEISENBERG_STAR_NCTSSOS_PARAMS.constraint)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first(pop, p, v.order;
            TS=v.ts, CS=v.cs, numeq=numeq,
            constraint=HEISENBERG_STAR_NCTSSOS_PARAMS.constraint)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

print_summary("HEISENBERG_STAR", results)
