# =============================================================================
# SOS Problem Type
# =============================================================================

"""
    SOSProblem{T}

A Sum of Squares (SOS) optimization problem in JuMP form.

This is the dual formulation of a moment problem. The SOS problem is typically
smaller and faster to solve than the primal moment problem.

# Type Parameters
- `T`: Coefficient type (Float64 for real problems)

# Fields
- `model::GenericModel{T}`: The JuMP model ready to be optimized

# Usage
```julia
# From moment problem
mp = moment_relax(pop, corr_sparsity, term_sparsities)
sos = sos_dualize(mp)

# Solve
set_optimizer(sos.model, Clarabel.Optimizer)
optimize!(sos.model)
obj = objective_value(sos.model)
```

See also: [`sos_dualize`](@ref), [`MomentProblem`](@ref)
"""
struct SOSProblem{T}
    model::GenericModel{T}
end


# =============================================================================
# Coefficient Extraction for Constraints
# =============================================================================

"""
    get_Cαj(unsymmetrized_basis::Vector{M}, localizing_mtx::Matrix{P}) where {T,M,P<:AbstractPolynomial{T}}

Extract coefficient matrix for dualization from a symbolic polynomial constraint matrix.

# Arguments
- `unsymmetrized_basis`: Sorted vector of monomials as the basis
- `localizing_mtx`: Polynomial-valued constraint matrix

# Returns
- Dictionary mapping (basis_idx, row, col) to coefficient values

This function decomposes the constraint matrix into a sparse representation
where each non-zero coefficient is indexed by the basis monomial position
and the matrix position (row, col).
"""
function get_Cαj(unsymmetrized_basis::Vector{M}, localizing_mtx::Matrix{P}) where {T,M,P<:AbstractPolynomial{T}}
    dim = size(localizing_mtx, 1)
    cis = CartesianIndices((dim, dim))

    # basis idx, row, col -> coefficient
    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()

    for ci in cis
        for (coeff, alpha) in terms(localizing_mtx[ci])
            dictionary_of_keys[(searchsortedfirst(unsymmetrized_basis, alpha), ci.I[1], ci.I[2])] = coeff
        end
    end

    return dictionary_of_keys
end


# =============================================================================
# SOS Dualization (Unified)
# =============================================================================

"""
    sos_dualize(mp::MomentProblem{A,T,M,P}) where {A,T,M,P} -> SOSProblem

Convert a symbolic moment problem into its dual SOS (Sum of Squares) problem.

# Arguments
- `mp::MomentProblem{A,T,M,P}`: The symbolic primal moment problem to dualize

# Returns
- `SOSProblem`: The dual SOS problem with matrix variables and constraints

# Description
The dualization process involves:
1. Creating matrix variables (G_j) for each constraint:
   - `:Zero` constraints -> SymmetricMatrixSpace (equality)
   - `:PSD` constraints -> PSDCone (real positive semidefinite)
   - `:HPSD` constraints -> PSDCone of size 2n x 2n (Hermitian embedding)
2. Introducing scalar variable `b` to bound the minimum value of the primal
3. Setting up polynomial equality constraints by matching coefficients
4. Returning the maximization of `b`

For complex algebras (Pauli, Fermionic, Bosonic), Hermitian PSD constraints
are embedded as real 2n x 2n PSD constraints using the standard construction:
```
H in HPSD <=> [Re(H), -Im(H); Im(H), Re(H)] in PSD
```

# Examples
```julia
mp = moment_relax(pop, corr_sparsity, term_sparsities)
sos = sos_dualize(mp)
set_optimizer(sos.model, Clarabel.Optimizer)
optimize!(sos.model)
obj = objective_value(sos.model)
```

See also: [`MomentProblem`](@ref), [`moment_relax`](@ref), [`SOSProblem`](@ref)
"""
function sos_dualize(mp::MomentProblem{A,TI,M,P}) where {A<:AlgebraType, TI<:Integer, M<:Monomial{A,TI}, P<:Polynomial{A,TI}}
    # Determine if complex based on algebra type
    is_complex = _is_complex_problem(A)

    if is_complex
        return _sos_dualize_hermitian(mp)
    else
        return _sos_dualize_real(mp)
    end
end


# =============================================================================
# Real SOS Dualization
# =============================================================================

"""
    _sos_dualize_real(mp::MomentProblem{A,TI,M,P}) -> SOSProblem

Internal: Dualize a real-valued symbolic moment problem.

For real algebras (NonCommutative, Projector, Unipotent), constraints use
standard PSD cones without complex embedding.
"""
function _sos_dualize_real(mp::MomentProblem{A,TI,M,P}) where {A<:AlgebraType, TI<:Integer, M<:Monomial{A,TI}, P<:Polynomial{A,TI}}
    # Get coefficient type from polynomial
    C = eltype(coefficients(mp.objective))

    dual_model = GenericModel{C}()

    # Create matrix variables for each constraint
    dual_variables = map(mp.constraints) do (cone, mat)
        G_dim = size(mat, 1)
        if cone == :Zero
            @variable(dual_model, [1:G_dim, 1:G_dim] in SymmetricMatrixSpace())
        else  # :PSD
            @variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
        end
    end

    # Scalar variable b to bound minimum
    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    # Symmetrize basis: moment problem uses unsymmetrized basis, but
    # we need to match coefficients of symmetric (canonicalized) monomials
    symmetric_basis = sorted_unique(symmetric_canon.(mp.total_basis))
    n_basis = length(symmetric_basis)

    # Create mapping from unsymmetrized basis to symmetric basis position
    basis_to_sym_idx = Dict(
        m => searchsortedfirst(symmetric_basis, symmetric_canon(m))
        for m in mp.total_basis
    )

    # Initialize constraint expressions: sum of coefficients must equal objective
    # For each basis monomial alpha, we have:
    #   sum_j sum_{k,l} C_alpha_jkl * G_j[k,l] = c_alpha - delta_{alpha,1} * b
    fα_constraints = [zero(GenericAffExpr{C,VariableRef}) for _ in 1:n_basis]

    # Add objective polynomial coefficients
    for (coef, mono) in zip(coefficients(mp.objective), monomials(mp.objective))
        sym_idx = searchsortedfirst(symmetric_basis, symmetric_canon(mono))
        if sym_idx <= n_basis && symmetric_basis[sym_idx] == symmetric_canon(mono)
            add_to_expression!(fα_constraints[sym_idx], coef)
        end
    end

    # Subtract b from the constant term (identity monomial)
    add_to_expression!(fα_constraints[1], -one(C), b)

    # Process each constraint matrix
    for (i, (_, mat)) in enumerate(mp.constraints)
        Cαjs = get_Cαj(sort(mp.total_basis), mat)
        for (ky, coef) in Cαjs
            basis_idx, row, col = ky
            # Map unsymmetrized basis index to symmetric index
            if basis_idx <= length(mp.total_basis)
                unsym_mono = sort(mp.total_basis)[basis_idx]
                sym_idx = basis_to_sym_idx[unsym_mono]
                add_to_expression!(fα_constraints[sym_idx], -coef, dual_variables[i][row, col])
            end
        end
    end

    # All coefficient constraints
    @constraint(dual_model, fα_constraints .== 0)

    return SOSProblem(dual_model)
end


# =============================================================================
# Hermitian SOS Dualization
# =============================================================================

"""
    _sos_dualize_hermitian(mp::MomentProblem{A,TI,M,P}) -> SOSProblem

Internal: Dualize a complex/Hermitian symbolic moment problem.

For complex algebras (Pauli, Fermionic, Bosonic), Hermitian PSD constraints
are embedded as real 2n x 2n PSD constraints.

The Hermitian embedding:
```
H in HPSD <=> [Re(H), -Im(H); Im(H), Re(H)] in PSD
```

For the dual SOS problem, we create 2n x 2n matrix variables and extract
the real and imaginary parts of the SOS multiplier from the block structure.
"""
function _sos_dualize_hermitian(mp::MomentProblem{A,TI,M,P}) where {A<:AlgebraType, TI<:Integer, M<:Monomial{A,TI}, P<:Polynomial{A,TI}}
    # Get coefficient type (should be complex for these algebras)
    C = eltype(coefficients(mp.objective))
    RC = real(C)  # Real coefficient type for dual model

    dual_model = GenericModel{RC}()

    # Create 2n x 2n matrix variables for Hermitian embedding
    dual_variables = map(mp.constraints) do (cone, mat)
        G_dim = size(mat, 1)
        if cone == :Zero
            @variable(dual_model, [1:2*G_dim, 1:2*G_dim] in SymmetricMatrixSpace())
        else  # :HPSD
            @variable(dual_model, [1:2*G_dim, 1:2*G_dim] in PSDCone())
        end
    end

    # Store dimensions for extracting blocks
    dual_variable_dims = [size(mp.constraints[i][2], 1) for i in 1:length(mp.constraints)]

    # Extract X1 = Re(G) = top-left + bottom-right blocks
    # Extract X2 = Im(G) = bottom-left - top-right blocks
    Xs = [
        # X1 (real part): top-left + bottom-right
        [
            begin
                dim = dual_variable_dims[i]
                dv[1:dim, 1:dim] .+ dv[dim+1:2*dim, dim+1:2*dim]
            end
            for (i, dv) in enumerate(dual_variables)
        ],
        # X2 (imaginary part): bottom-left - top-right
        [
            begin
                dim = dual_variable_dims[i]
                dv[dim+1:2*dim, 1:dim] .- dv[1:dim, dim+1:2*dim]
            end
            for (i, dv) in enumerate(dual_variables)
        ]
    ]

    # Scalar variable b
    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    # Sort basis for indexing
    symmetric_basis = sort(mp.total_basis)
    n_basis = length(symmetric_basis)

    # Real and imaginary parts of constraint expressions
    fα_constraints_re = [zero(GenericAffExpr{RC,VariableRef}) for _ in 1:n_basis]
    fα_constraints_im = [zero(GenericAffExpr{RC,VariableRef}) for _ in 1:n_basis]

    # Add objective polynomial coefficients (real and imaginary parts)
    for (coef, mono) in zip(coefficients(mp.objective), monomials(mp.objective))
        idx = searchsortedfirst(symmetric_basis, mono)
        if idx <= n_basis && symmetric_basis[idx] == mono
            add_to_expression!(fα_constraints_re[idx], real(coef))
            add_to_expression!(fα_constraints_im[idx], imag(coef))
        end
    end

    # Subtract b from the constant term (real part only)
    add_to_expression!(fα_constraints_re[1], -one(RC), b)

    # Process each constraint matrix
    for (i, (_, mat)) in enumerate(mp.constraints)
        Cαjs = get_Cαj(symmetric_basis, mat)
        for (ky, coef) in Cαjs
            basis_idx, row, col = ky

            # For complex coefficient coef = c_re + i*c_im
            # The contribution to real constraint from Re(G) is -c_re
            # The contribution to real constraint from Im(G) is +c_im
            # The contribution to imag constraint from Re(G) is -c_im
            # The contribution to imag constraint from Im(G) is -c_re

            # This comes from: coef * G = (c_re + i*c_im) * (X1 + i*X2)
            # Real part: c_re*X1 - c_im*X2
            # Imag part: c_im*X1 + c_re*X2

            c_re = real(coef)
            c_im = imag(coef)

            # Real constraint gets -c_re*X1[row,col] + c_im*X2[row,col]
            add_to_expression!(fα_constraints_re[basis_idx], -c_re, Xs[1][i][row, col])
            add_to_expression!(fα_constraints_re[basis_idx], +c_im, Xs[2][i][row, col])

            # Imag constraint gets -c_im*X1[row,col] - c_re*X2[row,col]
            add_to_expression!(fα_constraints_im[basis_idx], -c_im, Xs[1][i][row, col])
            add_to_expression!(fα_constraints_im[basis_idx], -c_re, Xs[2][i][row, col])
        end
    end

    # All coefficient constraints (real and imaginary parts)
    @constraint(dual_model, fα_constraints_re .== 0)
    @constraint(dual_model, fα_constraints_im .== 0)

    return SOSProblem(dual_model)
end


# =============================================================================
# Legacy Compatibility
# =============================================================================

# The old API had:
# - sos_dualize(moment_problem::MomentProblem{T,M}) for real (with JuMP model inside)
# - sos_dualize(cmp::ComplexMomentProblem{T,P}) for complex (symbolic)
#
# The new unified API has:
# - sos_dualize(mp::MomentProblem{A,T,M,P}) for both (all symbolic)
#
# ComplexMomentProblem is now an alias, so the old complex call signature still works.

"""
    ComplexMomentProblem{T, M, P}

Type alias for `MomentProblem{A,TI,M,P}` for backward compatibility.

In the legacy API, `ComplexMomentProblem` was a separate type for Hermitian
moment problems. In the new unified design, `MomentProblem{A,T,M,P}` handles
both real and complex cases based on the algebra type `A`.

# Deprecated
Use `MomentProblem` directly. This alias will be removed in Phase 4.

See also: [`MomentProblem`](@ref)
"""
const ComplexMomentProblem{T,M,P} = MomentProblem{A,TI,M,P} where {A<:AlgebraType, TI<:Integer}
