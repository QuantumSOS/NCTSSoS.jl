# Reconstructing Moment Matrices from Dual Constraints

## Your Approach (with corrections)

### Building the Index Mapping

You're correct about the approach to build `moment_matrices` as indices:

```julia
# Loop structure:
for clique_idx in 1:num_cliques
    for block_idx in 1:num_blocks_in_clique
        # Get the moment matrix constraint
        moment_matrix_constraint = model[:mom_mtx_clique_$(clique_idx)_block_$(block_idx)]

        # Get the matrix of JuMP variables
        moment_matrix_vars = ... # Extract from constraint

        # Build index matrix
        index_matrix = similar(moment_matrix_vars, Int)
        for i in 1:size(moment_matrix_vars, 1)
            for j in 1:size(moment_matrix_vars, 2)
                a = moment_matrix_vars[i, j]  # JuMP variable

                # Map to constraint index:
                unsym_idx = unsymmetrized_basis_vals_dict[a]
                monomial = unsymmetrized_basis[unsym_idx]
                constraint_idx = symmetrized_α2cons_dict[monomial]

                index_matrix[i, j] = constraint_idx
            end
        end

        moment_matrices[clique_idx][block_idx] = index_matrix
    end
end
```

### Simplified Mapping

Your one-liner is correct:
```julia
constraint_idx = symmetrized_α2cons_dict[unsymmetrized_basis[unsymmetrized_basis_vals_dict[a]]]
```

This can be simplified to:
```julia
# Given JuMP variable a:
unsym_idx = unsymmetrized_basis_vals_dict[a]
monomial = unsymmetrized_basis[unsym_idx]
constraint_idx = symmetrized_α2cons_dict[monomial]
```

## Important Corrections

### 1. Reconstructing Values - It Depends on What You Solved!

**If you solved the PRIMAL problem** (`dualize=false`):
```julia
# The moment variables have values:
moment_value = value(a)  # Direct value of the JuMP variable

# Reconstruct moment matrix:
reconstructed_moment_matrix[i,j] = value(moment_matrix_vars[i,j])
```

**If you solved the DUAL problem** (`dualize=true`):
```julia
# The dual constraints are EQUALITIES, not variables!
# You CANNOT do: value(dual_model[:coef_cons][idx])  ← This doesn't work!

# Instead, the dual constraints are of the form: f_α - b - Σ... = 0
# At optimality, this equals 0 (satisfied)

# To get moment values from the dual solution, you need to extract them
# from the dual variables G_j, which is more complex.
```

### 2. What the Dual Constraints Represent

The dual constraints `dual_model[:coef_cons]` are **constraint references**, not variable references:

```julia
# This is WRONG:
value(dual_model[:coef_cons][idx])  # Error! Constraints don't have values

# This is CORRECT (for primal problem):
primal_var = symmetric_variables[idx]
moment_value = value(primal_var)
```

### 3. Accessing Constraint Dual Values

If you solved the **primal** problem and want the **dual values** (shadow prices):

```julia
# For a constraint reference:
dual_value = dual(constraint_ref)

# For the polynomial constraints in the dual formulation,
# the dual values would be the optimal moment values!
```

## Correct Reconstruction Approaches

### Approach 1: Solve Primal Problem

```julia
# Solve primal
moment_problem = moment_relax(pop, corr_sparsity, term_sparsities)
optimize!(moment_problem.model)

# Reconstruct moment matrices directly
for clique_idx in 1:num_cliques
    for block_idx in 1:num_blocks
        constraint_name = Symbol("mom_mtx_clique_$(clique_idx)_block_$(block_idx)")

        # Get the constraint
        constraint_ref = moment_problem.model[constraint_name]

        # Extract the moment matrix variables from the constraint
        constraint_obj = constraint_object(constraint_ref)
        dim = get_dim(constraint_obj)

        # Reconstruct values
        moment_matrix_values = zeros(dim, dim)
        cis = CartesianIndices((dim, dim))

        for (ci, cur_expr) in zip(cis, constraint_obj.func)
            i, j = ci.I[1], ci.I[2]
            # cur_expr is the affine expression at position [i,j]
            moment_matrix_values[i,j] = value(cur_expr)
        end

        moment_matrices[clique_idx][block_idx] = moment_matrix_values
    end
end
```

### Approach 2: Build Index Mapping (for later use)

If you want to store the **structure** (indices) for later reconstruction:

```julia
function build_moment_matrix_index_structure(
    moment_problem::MomentProblem,
    corr_sparsity::CorrelativeSparsity,
    cliques_term_sparsities
)
    # Reconstruct the mapping dictionaries
    sa = moment_problem.sa
    unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))
    symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))

    symmetrized_α2cons_dict = Dict(
        zip(unsymmetrized_basis,
            map(x -> searchsortedfirst(symmetric_basis, canonicalize(x, sa)),
                unsymmetrized_basis))
    )

    unsymmetrized_basis_vals_dict = Dict(
        zip(getindex.(Ref(moment_problem.monomap), unsymmetrized_basis),
            1:length(unsymmetrized_basis))
    )

    # Storage for index matrices
    moment_matrices_indices = Vector{Vector{Matrix{Int}}}()

    # Loop over cliques and blocks
    for (clq_idx, term_sparsities) in enumerate(cliques_term_sparsities)
        clique_matrices = Vector{Matrix{Int}}()

        for (blk_idx, term_sparsity) in enumerate(term_sparsities)
            # Get the moment matrix constraint name
            constraint_name = Symbol("mom_mtx_clique_$(clq_idx)_block_$(blk_idx)")

            # Check if this constraint exists
            if !haskey(moment_problem.model.obj_dict, constraint_name)
                continue  # Skip if not a moment matrix (might be localizing)
            end

            # Get the constraint
            constraint_ref = moment_problem.model[constraint_name]
            constraint_obj = constraint_object(constraint_ref)
            dim = get_dim(constraint_obj)

            # Build index matrix
            index_matrix = zeros(Int, dim, dim)
            cis = CartesianIndices((dim, dim))

            for (ci, cur_expr) in zip(cis, constraint_obj.func)
                i, j = ci.I[1], ci.I[2]

                # Get the variable (assuming single variable per entry for moment matrices)
                # For moment matrices, each entry should be a single variable
                if length(cur_expr.terms) == 1
                    a = first(keys(cur_expr.terms))  # Get the JuMP variable

                    # Map to constraint index
                    unsym_idx = unsymmetrized_basis_vals_dict[a]
                    monomial = unsymmetrized_basis[unsym_idx]
                    constraint_idx = symmetrized_α2cons_dict[monomial]

                    index_matrix[i, j] = constraint_idx
                else
                    # Entry is a linear combination (e.g., in localizing matrices)
                    # Store -1 or handle specially
                    index_matrix[i, j] = -1
                end
            end

            push!(clique_matrices, index_matrix)
        end

        push!(moment_matrices_indices, clique_matrices)
    end

    return moment_matrices_indices
end

# Usage:
# moment_matrices_indices = build_moment_matrix_index_structure(
#     moment_problem, corr_sparsity, cliques_term_sparsities
# )
#
# Later, to get values from primal solution:
# value_matrix[i,j] = value(symmetric_variables[moment_matrices_indices[clq][blk][i,j]])
```

### Approach 3: Extract Moment Values from Dual Solution

If you solved the dual problem, you need to extract moment values differently:

```julia
# Solve dual
sos_problem = sos_dualize(moment_problem)
optimize!(sos_problem.model)

# The dual of the polynomial constraints gives the moment values!
# For each constraint in dual_model[:coef_cons], its dual value
# is the optimal value of the corresponding primal variable

sa = moment_problem.sa
unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))
symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))
symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

# Extract moment values from dual of constraints
moment_values = Dict()
for (i, var) in enumerate(symmetric_variables)
    # The dual of the i-th polynomial constraint is the optimal value of y_i
    dual_constraint = sos_problem.model[:coef_cons][i]
    moment_values[var] = dual(dual_constraint)
end

# Now reconstruct moment matrices using these values
# ... (similar to Approach 1, but using moment_values[var] instead of value(var))
```

## Complete Example Code

```julia
using NCTSSoS
using NCTSSoS: sos_dualize, moment_relax, sorted_unique, canonicalize, constraint_object, get_dim

# Create problem
@ncpolyvar x y
f = 2.0 - x^2 - y^2
g = 1.0 - x^2 - y^2
pop = polyopt(f; ineq_constraints=[g])

# Compute sparsity (simplified)
# ... (setup corr_sparsity and cliques_term_sparsities)

# Create moment problem
moment_problem = moment_relax(pop, corr_sparsity, cliques_term_sparsities)

# Option 1: Solve primal and get values directly
using Clarabel
set_optimizer(moment_problem.model, Clarabel.Optimizer)
optimize!(moment_problem.model)

# Extract moment matrix values
constraint_name = :mom_mtx_clique_1_block_1
constraint_ref = moment_problem.model[constraint_name]
constraint_obj = constraint_object(constraint_ref)
dim = get_dim(constraint_obj)

moment_matrix_values = zeros(dim, dim)
cis = CartesianIndices((dim, dim))
for (ci, cur_expr) in zip(cis, constraint_obj.func)
    i, j = ci.I[1], ci.I[2]
    moment_matrix_values[i,j] = value(cur_expr)
end

println("Moment matrix values:")
display(moment_matrix_values)

# Option 2: Build index structure then reconstruct
sa = moment_problem.sa
unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))
symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))
symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

symmetrized_α2cons_dict = Dict(
    zip(unsymmetrized_basis,
        map(x -> searchsortedfirst(symmetric_basis, canonicalize(x, sa)),
            unsymmetrized_basis))
)

unsymmetrized_basis_vals_dict = Dict(
    zip(getindex.(Ref(moment_problem.monomap), unsymmetrized_basis),
        1:length(unsymmetrized_basis))
)

# Build index matrix
index_matrix = zeros(Int, dim, dim)
for (ci, cur_expr) in zip(cis, constraint_obj.func)
    i, j = ci.I[1], ci.I[2]

    if length(cur_expr.terms) == 1
        a = first(keys(cur_expr.terms))

        # Your mapping:
        constraint_idx = symmetrized_α2cons_dict[
            unsymmetrized_basis[unsymmetrized_basis_vals_dict[a]]
        ]

        index_matrix[i, j] = constraint_idx
    end
end

println("\nIndex matrix (maps to symmetric_variables):")
display(index_matrix)

# Reconstruct using indices
reconstructed = zeros(dim, dim)
for i in 1:dim
    for j in 1:dim
        idx = index_matrix[i, j]
        reconstructed[i, j] = value(symmetric_variables[idx])
    end
end

println("\nReconstructed moment matrix:")
display(reconstructed)

println("\nVerification (should be same):")
println("Direct:        ", moment_matrix_values)
println("Reconstructed: ", reconstructed)
println("Match: ", isapprox(moment_matrix_values, reconstructed))
```

## Summary

Your approach is **correct for building the index structure**:
```julia
constraint_idx = symmetrized_α2cons_dict[unsymmetrized_basis[unsymmetrized_basis_vals_dict[a]]]
```

**However**, to reconstruct moment matrix **values**:

1. **If you solved the primal**: Use `value(a)` directly on the JuMP variable
2. **If you solved the dual**: Use `dual(dual_model[:coef_cons][idx])` to get moment values from constraint duals
3. **You cannot**: Use `value(dual_model[:coef_cons][idx])` - constraints don't have values!

The correct statement should be:
```julia
# Build index matrix:
constraint_idx = symmetrized_α2cons_dict[unsymmetrized_basis[unsymmetrized_basis_vals_dict[a]]]
moment_matrices_indices[clique][block][i,j] = constraint_idx

# Reconstruct values (if primal was solved):
reconstructed_value = value(symmetric_variables[constraint_idx])

# OR (if dual was solved):
reconstructed_value = dual(dual_model[:coef_cons][constraint_idx])
```
