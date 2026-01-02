# Heisenberg Star Graph Problem Definition
# =========================================
# Problem: Heisenberg model on star graph with unipotent constraint (U²=I)
# Variables: pij for each edge (i,j) in the star graph
# Objective: sum of pij over edges
# Constraints: Jordan-Wigner type triangle consistency relations
#
# The star graph has n sites with site 1 at center.
# Edge (1,2), (1,3), ..., (1,n)
#
# Expected optimal: -1.0 (tight for star graph)

# Sparsity variants
const HEISENBERG_STAR_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=1, num_sites=10),
    (name="CS", cs="MF", ts=false, order=1, num_sites=10),
    (name="Dense_n8", cs=false, ts=false, order=1, num_sites=8),
]

# NCTSSOS setup
# Note: NCTSSOS uses unipotent constraint
const HEISENBERG_STAR_NCTSSOS = (
    vars_expr = "num_sites = 10; vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]; @ncpolyvar p[1:length(vec_idx2ij)]",
    objective_expr = """
        # Build star graph edges and objective
        findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)
        edges_in_star = [(1, k) for k in 2:num_sites]
        obj = sum(p[findvaridx(i, j)] for (i, j) in edges_in_star)
    """,
    vars = "p",
    partition = :num_edges,
    constraint = "unipotent",
)

# NCTSSoS setup
const HEISENBERG_STAR_NCTSS = """
num_sites = 10
vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]
findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

registry, (pij,) = create_unipotent_variables([("pij", 1:length(vec_idx2ij))])

# Objective: sum over star graph edges (center=1)
objective = sum(1.0 * pij[[findvaridx(1, k) for k in 2:num_sites]])

# Triangle consistency constraints for Heisenberg model
gs = unique!([
    (
        pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
        pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
        pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
        pij[findvaridx(sort([i, k])...)] + 1.0
    ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
    (i != j && j != k && i != k)
])

pop = polyopt(objective, registry; eq_constraints=gs)
"""
