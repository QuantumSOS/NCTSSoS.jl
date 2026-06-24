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
- `n_unique_elements::Int`: Number of unique moment variables after canonicalization

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
    n_unique_elements::Int  # Number of unique moment variables after canonicalization
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

Note: Monomials not found in the basis (e.g., odd-parity fermionic monomials
filtered by superselection) are skipped since they have expectation value 0.
"""
function get_Cαj(unsymmetrized_basis::Vector{M}, localizing_mtx::Matrix{P}) where {T,M,P<:AbstractPolynomial{T}}
    dim = size(localizing_mtx, 1)
    cis = CartesianIndices((dim, dim))

    # basis idx, row, col -> coefficient
    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()

    n_basis = length(unsymmetrized_basis)

    for ci in cis
        for (coeff, alpha) in terms(localizing_mtx[ci])
            # `terms(::Polynomial)` yields `(coefficient, Monomial)` pairs; internal SOS basis
            # for moment problems is still in `NormalMonomial` space.
            alpha_key = alpha
            idx = searchsortedfirst(unsymmetrized_basis, alpha_key)
            # Skip monomials not in basis (e.g., odd-parity fermionic monomials)
            # These have expectation value 0 and don't contribute
            if idx <= n_basis && unsymmetrized_basis[idx] == alpha_key
                dictionary_of_keys[(idx, ci.I[1], ci.I[2])] = coeff
            end
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
1. Creating dual variables for each constraint:
   - `:Zero` constraints -> equality multipliers
   - `:PSD` constraints -> PSDCone (real positive semidefinite)
   - `:HPSD` constraints -> lifted real PSDCone of size `2n × 2n`
2. Building coefficient expressions for `objective - Aᴴ(dual)`
3. Maximizing the affine identity coefficient directly
4. Constraining all remaining coefficient expressions to zero

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
function sos_dualize(mp::MomentProblem{A,TI,M,P}) where {A<:AlgebraType, TI<:Integer, M<:NormalMonomial{A,TI}, P<:Polynomial{A,TI}}
    if _is_real_moment_problem(mp)
        return _sos_dualize_real(mp)
    elseif _is_complex_problem(A)
        return _sos_dualize_hermitian(mp)
    else
        return _sos_dualize_real(mp)
    end
end


# =============================================================================
# SOS Dualization Helpers for MomentLinearData
# =============================================================================

@inline _sos_coeff_type(::MomentLinearData{K,C,M}) where {K,C,M} = C
@inline _sos_real_type(::MomentLinearData{K,C,M}) where {K,C,M} = typeof(real(one(C)))

function _sos_moment_index(L::MomentLinearData{K}, key::K) where {K}
    return _get_key_value(L.moment_index, key, "moment index")
end

function _check_real_sos_cones!(mp::MomentProblem)
    for (cone, _) in mp.constraints
        (cone == :Zero || cone == :PSD) || error("Unexpected cone type $cone for real SOS dualization")
    end
    for block in mp.linear.psd_blocks_lin
        block.meta.cone == :PSD || error("Unexpected cached cone type $(block.meta.cone) for real SOS dualization")
    end
    return nothing
end

function _check_hermitian_sos_cones!(mp::MomentProblem)
    for (cone, mat) in mp.constraints
        if cone == :Zero
            _is_hermitian_poly_matrix(mat) || throw(ArgumentError(
                "Complex SOS dualization requires Hermitian zero-cone matrices. " *
                "Pass moment problems built through `moment_relax`, which splits " *
                "non-Hermitian equalities into Hermitian real/imaginary components first."
            ))
        elseif cone != :HPSD
            error("Unexpected cone type $cone for complex SOS dualization")
        end
    end
    for block in mp.linear.psd_blocks_lin
        block.meta.cone == :HPSD || error("Unexpected cached cone type $(block.meta.cone) for complex SOS dualization")
    end
    return nothing
end

function _accumulate_dual_contribution!(
    eqs::AbstractVector,
    idx::Integer,
    coef,
    dual_block,
    row::Integer,
    col::Integer,
    cone::Symbol,
)
    cone == :PSD || error("Real SOS dualization expected :PSD block, got $cone")
    add_to_expression!(eqs[idx], -coef, dual_block[row, col])
    return nothing
end

"""
    _accumulate_dual_contribution!(eqs_re, eqs_im, idx, coef, lifted, row, col, :HPSD, dim)

Accumulate one cached linear-form term into the Hermitian SOS coefficient
matching equations. This is the single home for the Hermitian real-lift adjoint
scaling: for a lifted dual matrix `Z`, the effective complex multiplier is
`X₁ + im*X₂` with `X₁ = Z₁₁ + Z₂₂` and `X₂ = Z₂₁ - Z₁₂`. The `Z₁₁ + Z₂₂`
sum is intentional; dropping either block is the classic factor-of-2 bug.
"""
function _accumulate_dual_contribution!(
    eqs_re::AbstractVector,
    eqs_im::AbstractVector,
    idx::Integer,
    coef,
    lifted,
    row::Integer,
    col::Integer,
    cone::Symbol,
    dim::Integer,
)
    cone == :HPSD || error("Hermitian SOS dualization expected :HPSD block, got $cone")

    n = Int(dim)
    i = Int(row)
    j = Int(col)
    c_re = real(coef)
    c_im = imag(coef)

    X1 = lifted[i, j] + lifted[n + i, n + j]
    X2 = lifted[n + i, j] - lifted[i, n + j]

    # Coefficient expressions are objective - Aᴴ(dual); the real identity
    # coefficient is maximized directly and the rest are constrained to zero.
    # Real(c * (X1 + im*X2)) = c_re*X1 - c_im*X2
    # Imag(c * (X1 + im*X2)) = c_im*X1 + c_re*X2
    add_to_expression!(eqs_re[idx], -c_re, X1)
    add_to_expression!(eqs_re[idx], +c_im, X2)
    add_to_expression!(eqs_im[idx], -c_im, X1)
    add_to_expression!(eqs_im[idx], -c_re, X2)
    return nothing
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
function _sos_dualize_real(mp::MomentProblem{A,TI,M,P}) where {A<:AlgebraType, TI<:Integer, M<:NormalMonomial{A,TI}, P<:Polynomial{A,TI}}
    _check_real_sos_cones!(mp)

    L = mp.linear
    C = _sos_coeff_type(L)
    dual_model = GenericModel{C}()

    psd_duals = Any[]
    for block in L.psd_blocks_lin
        push!(psd_duals, @variable(dual_model, [1:block.size, 1:block.size] in PSDCone()))
    end
    zero_duals = [@variable(dual_model) for _ in L.zero_constraints]

    fα_constraints = [zero(GenericAffExpr{C,VariableRef}) for _ in L.moments]

    for (key, coef) in L.objective_lin
        add_to_expression!(fα_constraints[_sos_moment_index(L, key)], coef)
    end

    for (block_idx, block) in enumerate(L.psd_blocks_lin)
        dual_block = psd_duals[block_idx]
        for i in 1:block.size, j in 1:block.size
            for (key, coef) in block.entries[i, j]
                _accumulate_dual_contribution!(
                    fα_constraints,
                    _sos_moment_index(L, key),
                    coef,
                    dual_block,
                    i,
                    j,
                    block.meta.cone,
                )
            end
        end
    end

    for (zc_idx, zc) in enumerate(L.zero_constraints)
        λ = zero_duals[zc_idx]
        for (key, coef) in zc.form
            add_to_expression!(fα_constraints[_sos_moment_index(L, key)], -coef, λ)
        end
    end

    # Eliminate the usual scalar bound variable: the real identity coefficient is b.
    identity_idx = _sos_moment_index(L, L.identity)
    @objective(dual_model, Max, fα_constraints[identity_idx])

    coefficient_indices = [i for i in eachindex(fα_constraints) if i != identity_idx]
    isempty(coefficient_indices) || @constraint(dual_model, fα_constraints[coefficient_indices] .== 0)

    return SOSProblem(dual_model, length(L.moments))
end


# =============================================================================
# Hermitian SOS Dualization
# =============================================================================

"""
    sos_dualize(mp::StateMomentProblem{A,ST,TI,M,P}) -> SOSProblem

Convert a symbolic state moment problem into its dual SOS problem.

State polynomial optimization problems (with Unipotent, Projector algebras) use
real PSD constraints without complex embedding.

# Arguments
- `mp::StateMomentProblem`: The symbolic state moment problem

# Returns
- `SOSProblem`: The dual SOS problem
"""
function sos_dualize(mp::StateMomentProblem{A,ST,TI,M,P}) where {A<:AlgebraType, ST<:StateType, TI<:Integer, M<:NCStateWord{ST,A,TI}, P<:NCStatePolynomial}
    # State polynomial optimization uses real PSD constraints
    return _sos_dualize_state(mp)
end

"""
    _sos_dualize_state(mp::StateMomentProblem) -> SOSProblem

Internal: Dualize a state polynomial moment problem.

State polynomial optimization with algebras like UnipotentAlgebra and ProjectorAlgebra
uses real-valued SDP with standard PSD constraints.
"""
function _sos_dualize_state(mp::StateMomentProblem{A,ST,TI,M,P}) where {A<:AlgebraType, ST<:StateType, TI<:Integer, M<:NCStateWord{ST,A,TI}, P<:NCStatePolynomial}
    # Get coefficient type from polynomial
    C = eltype(coefficients(mp.objective))

    dual_model = GenericModel{C}()

    # Create matrix variables for each constraint (now includes block_basis in tuple)
    dual_variables = map(mp.constraints) do (cone, mat, _block_basis)
        G_dim = size(mat, 1)
        if cone == :Zero
            @variable(dual_model, [1:G_dim, 1:G_dim] in SymmetricMatrixSpace())
        else  # :PSD
            @variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
        end
    end

    # For state polynomial optimization, we work with expectation values (StateWords)
    # The key insight: expval(<I>*xy) = <xy> = expval(<xy>*I)
    # So we convert NCStateWords to StateWords via expval for comparison

    # Create canonical StateWord basis from total_basis
    # Apply expval to convert NCStateWord to StateWord, then symmetric_canon
    SW = StateWord{ST,A,TI}
    state_basis = _sorted_stateword_basis_from_ncsw(mp.total_basis)
    n_basis = length(state_basis)
    sw_to_idx = Dict(sw => i for (i, sw) in enumerate(state_basis))

    identity_sw = one(SW)
    identity_idx = searchsortedfirst(state_basis, identity_sw)
    if identity_idx > n_basis || state_basis[identity_idx] != identity_sw
        error("Identity StateWord not found in basis - this shouldn't happen")
    end

    # Initialize constraint expressions
    fα_constraints = [zero(GenericAffExpr{C,VariableRef}) for _ in 1:n_basis]

    # Add objective polynomial coefficients
    # Objective NCStateWords have form <xy>*I, we need to match them to basis StateWords
    missing_objective_words = SW[]
    for (coef, ncsw) in zip(coefficients(mp.objective), monomials(mp.objective))
        canon_sw = symmetric_canon(expval(ncsw))
        sym_idx = get(sw_to_idx, canon_sw, 0)
        if iszero(sym_idx)
            push!(missing_objective_words, canon_sw)
            continue
        end
        add_to_expression!(fα_constraints[sym_idx], coef)
    end
    _throw_missing_state_words(missing_objective_words, "the objective"; source="Relaxation basis")

    # Process each constraint matrix by iterating directly over block (row, col) pairs
    # This ensures all StateWord contributions are accumulated correctly across blocks
    # (matching NCTSSOS's on-the-fly loop approach)
    missing_constraint_words = SW[]
    for (i, (_, mat, _block_basis)) in enumerate(mp.constraints)
        dim = size(mat, 1)
        # Iterate over upper triangle to avoid double-counting
        for row in 1:dim
            for col in row:dim
                poly = mat[row, col]
                for (coeff, ncsw) in zip(coefficients(poly), monomials(poly))
                    # Convert NCStateWord to canonical StateWord
                    sw = symmetric_canon(expval(ncsw))
                    sw_idx = get(sw_to_idx, sw, 0)
                    if iszero(sw_idx)
                        push!(missing_constraint_words, sw)
                        continue
                    end
                    # Multiply off-diagonal entries by 2 for symmetric matrix contribution
                    effective_coeff = (row == col) ? coeff : 2 * coeff
                    add_to_expression!(fα_constraints[sw_idx], -effective_coeff, dual_variables[i][row, col])
                end
            end
        end
    end
    _throw_missing_state_words(missing_constraint_words, "the constraint matrices"; source="Relaxation basis")

    # Eliminate the usual scalar bound variable: the identity coefficient is b.
    @objective(dual_model, Max, fα_constraints[identity_idx])

    coefficient_indices = [i for i in eachindex(fα_constraints) if i != identity_idx]
    isempty(coefficient_indices) || @constraint(dual_model, fα_constraints[coefficient_indices] .== 0)

    return SOSProblem(dual_model, n_basis)
end

"""
    _sos_dualize_hermitian(mp::MomentProblem{A,TI,M,P}) -> SOSProblem

Internal: Dualize a complex/Hermitian symbolic moment problem using the cached
`mp.linear::MomentLinearData` forms.

For complex algebras (Pauli, Fermionic, Bosonic), each cached `:HPSD` block gets
a lifted real PSD dual matrix. Cached scalar zero constraints get real free
multipliers. The Hermitian factor-of-2 convention lives only in
`_accumulate_dual_contribution!`.
"""
function _sos_dualize_hermitian(mp::MomentProblem{A,TI,M,P}) where {A<:AlgebraType, TI<:Integer, M<:NormalMonomial{A,TI}, P<:Polynomial{A,TI}}
    _check_hermitian_sos_cones!(mp)

    L = mp.linear
    RC = _sos_real_type(L)
    dual_model = GenericModel{RC}()

    psd_duals = Any[]
    for block in L.psd_blocks_lin
        dim = block.size
        push!(psd_duals, (
            dim = dim,
            lifted = @variable(dual_model, [1:2*dim, 1:2*dim] in PSDCone()),
        ))
    end
    zero_duals = [@variable(dual_model) for _ in L.zero_constraints]

    fα_constraints_re = [zero(GenericAffExpr{RC,VariableRef}) for _ in L.moments]
    fα_constraints_im = [zero(GenericAffExpr{RC,VariableRef}) for _ in L.moments]

    for (key, coef) in L.objective_lin
        idx = _sos_moment_index(L, key)
        add_to_expression!(fα_constraints_re[idx], real(coef))
        add_to_expression!(fα_constraints_im[idx], imag(coef))
    end

    for (block_idx, block) in enumerate(L.psd_blocks_lin)
        dual = psd_duals[block_idx]
        for i in 1:block.size, j in 1:block.size
            for (key, coef) in block.entries[i, j]
                _accumulate_dual_contribution!(
                    fα_constraints_re,
                    fα_constraints_im,
                    _sos_moment_index(L, key),
                    coef,
                    dual.lifted,
                    i,
                    j,
                    block.meta.cone,
                    dual.dim,
                )
            end
        end
    end

    for (zc_idx, zc) in enumerate(L.zero_constraints)
        λ = zero_duals[zc_idx]
        for (key, coef) in zc.form
            idx = _sos_moment_index(L, key)
            add_to_expression!(fα_constraints_re[idx], -real(coef), λ)
            add_to_expression!(fα_constraints_im[idx], -imag(coef), λ)
        end
    end

    # Eliminate only the real identity equation; the imaginary identity equation
    # still enforces a real scalar bound.
    identity_idx = _sos_moment_index(L, L.identity)
    @objective(dual_model, Max, fα_constraints_re[identity_idx])

    real_coefficient_indices = [i for i in eachindex(fα_constraints_re) if i != identity_idx]
    isempty(real_coefficient_indices) || @constraint(dual_model, fα_constraints_re[real_coefficient_indices] .== 0)
    @constraint(dual_model, fα_constraints_im .== 0)

    return SOSProblem(dual_model, length(L.moments))
end
