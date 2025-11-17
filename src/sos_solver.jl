struct SOSProblem{T}
    model::GenericModel{T}
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

    symmetrized_α2cons_dict = Dict(zip(unsymmetrized_basis, map(x -> searchsortedfirst(symmetric_basis, canonicalize(x, moment_problem.sa)), unsymmetrized_basis)))

    unsymmetrized_basis_vals_dict = Dict(zip(getindex.(Ref(moment_problem.monomap), unsymmetrized_basis), 1:length(unsymmetrized_basis)))

    add_to_expression!(fα_constraints[1], -one(T), b)

    for (i, sdp_constraint) in enumerate(moment_problem.constraints)
        Cαjs = get_Cαj(unsymmetrized_basis_vals_dict, constraint_object(sdp_constraint))
        for (ky, coef) in Cαjs
            add_to_expression!(fα_constraints[symmetrized_α2cons_dict[unsymmetrized_basis[ky[1]]]], -coef, dual_variables[i][ky[2], ky[3]])
        end
    end
    @constraint(dual_model, fα_constraints .== 0)

    return SOSProblem(dual_model)
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
    return SOSProblem(dual_model)
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
