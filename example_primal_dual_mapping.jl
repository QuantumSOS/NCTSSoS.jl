"""
Practical Example: Mapping Primal JuMP Variables to Dual Constraints

This example shows how to access the correspondence between:
- Primal moment variables (in MomentProblem)
- Dual polynomial equality constraints (in SOSProblem)
"""

using NCTSSoS
using NCTSSoS: sos_dualize, moment_relax, sorted_unique, canonicalize

# Create a simple problem
@ncpolyvar x y

f = 2.0 - x^2 - y^2
g = 1.0 - x^2 - y^2

pop = polyopt(f; ineq_constraints=[g])

println("="^70)
println("Original Problem")
println("="^70)
println("Objective: $f")
println("Constraint: $g ≥ 0")

# Create moment problem (assuming you have the sparsity structures)
# For this example, we'll use a simple no-sparsity setup
using NCTSSoS: CorrelativeSparsity, TermSparsity, compute_basis

# Simple correlative sparsity (no exploitation)
corr_sparsity = CorrelativeSparsity(
    clq_cons = [[1]],  # One clique with constraint 1
    global_cons = Int[],
    cons = [g]
)

# Simple term sparsity (order 1 basis)
order = 1
basis = compute_basis([x, y], order)

# Create term sparsity for objective and constraint
term_sparsity = TermSparsity(
    term_sparse_graph_supp = basis,
    block_bases = [basis]  # One block with full basis
)

cliques_term_sparsities = [[term_sparsity, term_sparsity]]  # For objective and constraint

# Create moment problem
moment_problem = moment_relax(pop, corr_sparsity, cliques_term_sparsities)

println("\n" * "="^70)
println("Moment Problem (Primal)")
println("="^70)

println("\nMonomap (Monomial → JuMP Variable):")
for (mono, var) in sort(collect(moment_problem.monomap), by=x->string(x[1]))
    println("  $mono → $var")
end

# Dualize to SOS problem
sos_problem = sos_dualize(moment_problem)
dual_model = sos_problem.model

println("\n" * "="^70)
println("SOS Problem (Dual)")
println("="^70)

println("\nDual model objective:")
println("  ", objective_function(dual_model))

println("\nDual constraints stored as :coef_cons")

# ============================================================================
# BUILD THE MAPPING
# ============================================================================

sa = moment_problem.sa
unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))
symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))
symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

println("\n" * "="^70)
println("Primal-Dual Mapping")
println("="^70)

println("\nTotal unique monomials: $(length(symmetric_basis))")
println("Total primal variables: $(length(symmetric_variables))")
println("Total dual constraints: $(length(dual_model[:coef_cons]))")

println("\n" * "-"^70)
println("Complete Mapping Table")
println("-"^70)
println(rpad("Index", 8) * rpad("Monomial", 20) * rpad("Primal Variable", 20) * "Dual Constraint")
println("-"^70)

for (i, (mono, pvar)) in enumerate(zip(symmetric_basis, symmetric_variables))
    dcons = dual_model[:coef_cons][i]
    println(rpad(string(i), 8) *
            rpad(string(mono), 20) *
            rpad(string(pvar), 20) *
            string(dcons))
end

# ============================================================================
# EXAMPLE QUERIES
# ============================================================================

println("\n" * "="^70)
println("Example Queries")
println("="^70)

# Query 1: Find dual constraint for a specific primal variable
println("\n--- Query 1: Find dual constraint for primal variable y[3] ---")

if length(symmetric_variables) >= 3
    target_var = symmetric_variables[3]
    constraint_idx = findfirst(==(target_var), symmetric_variables)

    println("Primal variable: $target_var")
    println("Represents monomial: $(symmetric_basis[constraint_idx])")
    println("Dual constraint index: $constraint_idx")
    println("Dual constraint: $(dual_model[:coef_cons][constraint_idx])")
end

# Query 2: Find primal variable and dual constraint for a monomial
println("\n--- Query 2: Find mappings for monomial x^2 ---")

target_mono = x^2
canonical_mono = canonicalize(target_mono, sa)
constraint_idx = searchsortedfirst(symmetric_basis, canonical_mono)

if constraint_idx <= length(symmetric_basis) && symmetric_basis[constraint_idx] == canonical_mono
    println("Original monomial: $target_mono")
    println("Canonical form: $canonical_mono")
    println("Constraint index: $constraint_idx")
    println("Primal variable: $(symmetric_variables[constraint_idx])")
    println("Dual constraint: $(dual_model[:coef_cons][constraint_idx])")
else
    println("Monomial $target_mono not found in basis")
end

# Query 3: Build a complete mapping dictionary
println("\n--- Query 3: Build mapping dictionary ---")

mapping = Dict{Any, NamedTuple}()
for (i, (mono, pvar)) in enumerate(zip(symmetric_basis, symmetric_variables))
    mapping[pvar] = (
        constraint_index = i,
        monomial = mono,
        dual_constraint = dual_model[:coef_cons][i]
    )
end

println("Mapping dictionary created with $(length(mapping)) entries")
println("\nExample access:")
sample_var = symmetric_variables[1]
println("  mapping[$sample_var] = $(mapping[sample_var])")

# ============================================================================
# VERIFY DUALITY
# ============================================================================

println("\n" * "="^70)
println("Verify Strong Duality")
println("="^70)

println("\nThis mapping ensures that:")
println("  1. Each primal variable y_α corresponds to exactly one dual constraint")
println("  2. The dual constraint enforces the polynomial equality for monomial α")
println("  3. The correspondence is maintained through index alignment")
println("\nWhen both problems are solved:")
println("  primal_optimal = dual_optimal (by SDP strong duality)")

println("\n" * "="^70)
println("Helper Functions")
println("="^70)

# Define helper functions
function get_dual_constraint_index(primal_var, symmetric_variables)
    return findfirst(==(primal_var), symmetric_variables)
end

function get_monomial_for_variable(primal_var, symmetric_basis, symmetric_variables)
    idx = get_dual_constraint_index(primal_var, symmetric_variables)
    return idx !== nothing ? symmetric_basis[idx] : nothing
end

println("""
Available helper functions:

1. get_dual_constraint_index(primal_var, symmetric_variables)
   Returns the index i such that:
   - primal_var = symmetric_variables[i]
   - dual_constraint = dual_model[:coef_cons][i]

2. get_monomial_for_variable(primal_var, symmetric_basis, symmetric_variables)
   Returns the monomial corresponding to primal_var

Example usage:
  idx = get_dual_constraint_index(y[2], symmetric_variables)
  mono = get_monomial_for_variable(y[2], symmetric_basis, symmetric_variables)
""")

println("\n" * "="^70)
println("Summary")
println("="^70)
println("""
The mapping works as follows:

  Index i | Monomial           | Primal Variable      | Dual Constraint
  --------|-------------------|---------------------|------------------
  1       | symmetric_basis[1] | symmetric_variables[1] | dual_model[:coef_cons][1]
  2       | symmetric_basis[2] | symmetric_variables[2] | dual_model[:coef_cons][2]
  ...     | ...               | ...                 | ...

Key arrays:
- symmetric_basis: Vector of unique monomials (after canonicalization)
- symmetric_variables: Vector of primal JuMP variables
- dual_model[:coef_cons]: Vector of dual constraint references

All three arrays use the SAME indexing scheme!
""")
