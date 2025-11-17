struct SOSProblem{T}
    model::GenericModel{T}
    moment_matrices::Union{Nothing, Vector{Vector{Matrix{Int}}}}
end

# Decompose the matrix into the form sum_j C_αj * g_j
# j: index of the constraint
# α: the monomial (JuMP variable)
function get_Cαj(basis_dict::Dict{GenericVariableRef{T},Int}, localizing_mtx::VectorConstraint{F,S,Shape}) where {T,F,S,Shape}
    dim = get_dim(localizing_mtx)
    cis = CartesianIndices((dim, dim))

    # basis idx, row, col
    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()

    for (ci, cur_expr) in zip(cis, localizing_mtx.func)
        for (α, coeff) in cur_expr.terms
            dictionary_of_keys[(basis_dict[α], ci.I[1], ci.I[2])] = coeff
        end
    end

    return dictionary_of_keys
end

"""
    sos_dualize(moment_problem::MomentProblem{T,M}) where {T,M} -> SOSProblem

Convert a moment problem into its dual SOS (Sum of Squares) problem formulation.

This function takes a moment problem and constructs the corresponding dual optimization
problem by introducing matrix variables for each constraint and formulating the dual
constraints that ensure polynomial equality.

# Arguments
- `moment_problem::MomentProblem{T,M}`: The primal moment problem to dualize

# Returns
- `SOSProblem`: The dual SOS problem with matrix variables and constraints

# Details
The dualization process involves:
1. Creating matrix variables (G_j) for each constraint, either in symmetric matrix space
   or positive semidefinite cone depending on the constraint type
2. Introducing a scalar variable `b` to bound the minimum value of the primal problem
3. Setting up polynomial equality constraints by matching coefficients of monomials
4. Handling symmetrization of the monomial basis to ensure proper polynomial comparison

The resulting dual problem maximizes `b` subject to the constraint that the sum of
matrix variables weighted by coefficient matrices equals the objective polynomial.
"""
function sos_dualize(moment_problem::MomentProblem{T,M}) where {T,M}
    dual_model = GenericModel{T}()

    # Initialize Gj as variables
    dual_variables = map(constraint_object.(moment_problem.constraints)) do cons
        G_dim = get_dim(cons)
        @variable(dual_model, [1:G_dim, 1:G_dim] in ((cons.set isa MOI.Zeros) ? SymmetricMatrixSpace() : PSDCone()))
    end
    # b: to bound the minimum value of the primal problem
    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    primal_objective_terms = objective_function(moment_problem.model).terms

    # NOTE: objective is Symmetric, hence when comparing polynomials, we need to canonicalize them first
    unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))

    symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(moment_problem.sa)))

    # JuMP variables corresponding to symmetric_basis
    symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)


    fα_constraints = [GenericAffExpr{T,VariableRef}(get(primal_objective_terms, α, zero(T))) for α in symmetric_variables]

    # Monomial ->  index in symmetric_basis
    symmetrized_α2cons_dict = Dict(zip(unsymmetrized_basis, map(x -> searchsortedfirst(symmetric_basis, canonicalize(x, moment_problem.sa)), unsymmetrized_basis)))

    # model's jump var -> index in unsymmetrized_basis
    unsymmetrized_basis_vals_dict = Dict(zip(getindex.(Ref(moment_problem.monomap), unsymmetrized_basis), 1:length(unsymmetrized_basis)))

    add_to_expression!(fα_constraints[1], -one(T), b)

    for (i, sdp_constraint) in enumerate(moment_problem.constraints)
        Cαjs = get_Cαj(unsymmetrized_basis_vals_dict, constraint_object(sdp_constraint))
        for (ky, coef) in Cαjs
            add_to_expression!(fα_constraints[symmetrized_α2cons_dict[unsymmetrized_basis[ky[1]]]], -coef, dual_variables[i][ky[2], ky[3]])
        end
    end

    dual_model[:coef_cons] = @constraint(dual_model, fα_constraints .== 0)

    # Build moment matrix index structure
    moment_matrices_indices = build_moment_matrices_indices(
        moment_problem,
        symmetrized_α2cons_dict,
        unsymmetrized_basis_vals_dict,
        unsymmetrized_basis
    )

    return SOSProblem(dual_model, moment_matrices_indices)
end

function sos_dualize(cmp::ComplexMomentProblem{T,P}) where {T,P}
    dual_model = GenericModel{real(T)}()

    dual_variables = map(cmp.constraints) do (type,cons)
        G_dim = size(cons,1)
        @variable(dual_model, [1:2*G_dim, 1:2*G_dim] in (type == :Zero ? SymmetricMatrixSpace() : PSDCone()))
    end
    dual_variable_dims = map(dual_variables) do dv
        size(dv, 1) ÷ 2
    end

    # a little allocation may go a long way?
    # X1 + X2 and X3 - X3^T
    Xs = [[
            begin
                dim = dual_variable_dims[i]
                dv[1:dim, 1:dim] .+ dv[dim+1:2*dim, dim+1:2*dim]
            end for (i, dv) in enumerate(dual_variables)
        ], [
            begin
                dim = dual_variable_dims[i]
                dv[dim+1:2*dim, 1:dim] .- dv[1:dim, 1+dim:2*dim]
            end for (i, dv) in enumerate(dual_variables)
        ]
    ]

    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    symmetric_basis = sort(cmp.total_basis)


    # real and imag parts of fα constraints
    fα_constraints = [[zero(GenericAffExpr{T,VariableRef}) for _ in 1:length(symmetric_basis)],[zero(GenericAffExpr{T,VariableRef}) for _ in 1:length(symmetric_basis)]]

    for (coef,mono) in terms(cmp.objective)
        for (fα_constraints_part, part_func) in zip(fα_constraints, [real, imag])
            fα_constraints_part[searchsortedfirst(symmetric_basis, mono)] += part_func(coef)
        end
    end

    add_to_expression!(fα_constraints[1][1], -one(T), b)

    for (i, (_,sdp_constraint)) in enumerate(cmp.constraints)
        Cαjs = get_Cαj(cmp.total_basis, sdp_constraint)
        for (ky, coef) in Cαjs
            for (X_part, coef_part, sign, part_func) in zip([1, 2, 2, 1], [1, 1, 2, 2], [-1, -1, -1, 1], [real, imag, real, imag])
                add_to_expression!(fα_constraints[coef_part][ky[1]], sign*part_func(coef), Xs[X_part][i][ky[2], ky[3]])
            end
        end
    end
    @constraint(dual_model, fα_constraints[1] .== 0)
    @constraint(dual_model, fα_constraints[2] .== 0)
    # TODO: Implement moment_matrices_indices for ComplexMomentProblem
    return SOSProblem(dual_model, nothing)
end

function get_Cαj(unsymmetrized_basis::Vector{M}, localizing_mtx::Matrix{P}) where {T,M,P<:AbstractPolynomial{T}}
    dim = size(localizing_mtx,1)
    cis = CartesianIndices((dim, dim))

    # basis idx, row, col
    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()
    for ci in cis
        for (coeff,α) in terms(localizing_mtx[ci])
            dictionary_of_keys[(searchsortedfirst(unsymmetrized_basis,α), ci.I[1], ci.I[2])] = coeff
        end
    end
    return dictionary_of_keys
end

"""
    build_moment_matrices_indices(
        moment_problem::MomentProblem,
        symmetrized_α2cons_dict::Dict,
        unsymmetrized_basis_vals_dict::Dict,
        unsymmetrized_basis::Vector
    ) -> Vector{Vector{Matrix{Int}}}

Build index matrices that map moment matrix positions to constraint indices in the dual problem.

# Arguments
- `moment_problem::MomentProblem`: The primal moment problem
- `symmetrized_α2cons_dict`: Dictionary mapping monomials to constraint indices in symmetric_basis
- `unsymmetrized_basis_vals_dict`: Dictionary mapping JuMP variables to indices in unsymmetrized_basis
- `unsymmetrized_basis`: Vector of all monomials before symmetrization

# Returns
- `Vector{Vector{Matrix{Int}}}`: A nested structure where:
  - First level: cliques
  - Second level: blocks within each clique
  - Third level: Matrix of constraint indices, where `result[clq][blk][i,j]` gives
    the index in `dual_model[:coef_cons]` for the moment at position (i,j)

# Description
This function examines all moment matrix constraints in the primal problem and builds
an index structure that allows efficient reconstruction of moment matrix values from
dual constraint values after optimization.

Only moment matrices (not localizing matrices) are included in the output structure.
The naming pattern `mom_mtx_clique_\$(clq_idx)_block_\$(blk_idx)` is used to identify
and organize moment matrices.
"""
function build_moment_matrices_indices(
    moment_problem::MomentProblem,
    symmetrized_α2cons_dict::Dict,
    unsymmetrized_basis_vals_dict::Dict,
    unsymmetrized_basis::Vector
)
    model = moment_problem.model

    # Collect all moment matrix constraint names and parse their structure
    moment_matrix_info = Dict{Tuple{Int,Int}, Symbol}()  # (clq_idx, blk_idx) => constraint_name

    for (name, obj) in model.obj_dict
        name_str = string(name)
        # Match pattern: mom_mtx_clique_X_block_Y
        if startswith(name_str, "mom_mtx_clique_")
            m = match(r"mom_mtx_clique_(\d+)_block_(\d+)", name_str)
            if m !== nothing
                clq_idx = parse(Int, m.captures[1])
                blk_idx = parse(Int, m.captures[2])
                moment_matrix_info[(clq_idx, blk_idx)] = name
            end
        end
    end

    if isempty(moment_matrix_info)
        return Vector{Vector{Matrix{Int}}}()  # No moment matrices found
    end

    # Determine structure: max clique index and blocks per clique
    max_clq = maximum(k[1] for k in keys(moment_matrix_info))

    # Initialize result structure
    result = Vector{Vector{Matrix{Int}}}(undef, max_clq)
    for i in 1:max_clq
        result[i] = Matrix{Int}[]
    end

    # Process each moment matrix
    for ((clq_idx, blk_idx), constraint_name) in sort(collect(moment_matrix_info))
        constraint_ref = model[constraint_name]
        constraint_obj = constraint_object(constraint_ref)
        dim = get_dim(constraint_obj)

        # Build index matrix for this moment matrix
        index_matrix = zeros(Int, dim, dim)
        cis = CartesianIndices((dim, dim))

        for (ci, cur_expr) in zip(cis, constraint_obj.func)
            i, j = ci.I[1], ci.I[2]

            # Extract the JuMP variable(s) from the expression
            # For moment matrices (poly_idx == 1), each entry should be a single variable
            if length(cur_expr.terms) == 1
                a = first(keys(cur_expr.terms))  # Get the JuMP variable

                # Map variable -> monomial -> constraint index
                unsym_idx = unsymmetrized_basis_vals_dict[a]
                monomial = unsymmetrized_basis[unsym_idx]
                constraint_idx = symmetrized_α2cons_dict[monomial]

                index_matrix[i, j] = constraint_idx
            elseif length(cur_expr.terms) == 0
                # Empty expression (zero entry)
                index_matrix[i, j] = 0
            else
                # Multiple terms - shouldn't happen for moment matrices
                # Store -1 to indicate this case
                index_matrix[i, j] = -1
            end
        end

        # Ensure we have enough space in the clique's block vector
        while length(result[clq_idx]) < blk_idx
            push!(result[clq_idx], zeros(Int, 0, 0))
        end

        result[clq_idx][blk_idx] = index_matrix
    end

    return result
end

"""
    get_moment_matrices(sos_problem::SOSProblem) -> Vector{Vector{Matrix{Float64}}}

Reconstruct moment matrix values from a solved dual SOS problem.

# Arguments
- `sos_problem::SOSProblem`: A solved SOS problem (must have been optimized)

# Returns
- `Vector{Vector{Matrix{Float64}}}`: Moment matrix values with the same structure as
  `moment_matrices_indices`, where `result[clq][blk][i,j]` contains the optimal
  moment value for position (i,j)

# Description
After solving the dual SOS problem, this function extracts the moment values from the
dual values of the polynomial equality constraints. The dual values of these constraints
correspond to the optimal values of the primal moment variables.

# Throws
- Error if `sos_problem.moment_matrices` is `nothing`
- Error if the model has not been optimized

# Example
```julia
# Solve the dual problem
sos_problem = sos_dualize(moment_problem)
set_optimizer(sos_problem.model, Clarabel.Optimizer)
optimize!(sos_problem.model)

# Extract moment matrices
moment_matrices = get_moment_matrices(sos_problem)

# Access values
value_at_pos = moment_matrices[clique_idx][block_idx][i, j]
```
"""
function get_moment_matrices(sos_problem::SOSProblem{T}) where {T}
    if sos_problem.moment_matrices === nothing
        error("moment_matrices is nothing. Cannot reconstruct moment matrices.")
    end

    if termination_status(sos_problem.model) == MOI.OPTIMIZE_NOT_CALLED
        error("Model has not been optimized. Call optimize!(sos_problem.model) first.")
    end

    moment_matrices_indices = sos_problem.moment_matrices
    coef_cons = sos_problem.model[:coef_cons]

    # Reconstruct values using dual values of constraints
    result = Vector{Vector{Matrix{T}}}(undef, length(moment_matrices_indices))

    for clq_idx in 1:length(moment_matrices_indices)
        clique_matrices = Vector{Matrix{T}}(undef, length(moment_matrices_indices[clq_idx]))

        for blk_idx in 1:length(moment_matrices_indices[clq_idx])
            index_matrix = moment_matrices_indices[clq_idx][blk_idx]
            dim = size(index_matrix, 1)

            value_matrix = zeros(T, dim, dim)

            for i in 1:dim
                for j in 1:dim
                    cons_idx = index_matrix[i, j]
                    if cons_idx > 0
                        # Extract dual value of the polynomial constraint
                        # The dual of the constraint gives the primal variable value
                        value_matrix[i, j] = dual(coef_cons[cons_idx])
                    end
                    # If cons_idx <= 0, leave as zero
                end
            end

            clique_matrices[blk_idx] = value_matrix
        end

        result[clq_idx] = clique_matrices
    end

    return result
end
