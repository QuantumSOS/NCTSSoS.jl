# =============================================================================
# Unified Symbolic MomentProblem
# =============================================================================

"""
    MomentProblem{A<:AlgebraType, T<:Integer, M<:NormalMonomial{A,T}, P<:Polynomial{A,T}}

A symbolic representation of a moment relaxation problem.

This unified type handles both real and complex (Hermitian) moment problems.
The algebra type `A` determines which cone type to use when solving:
- Real algebras (NonCommutative, Projector, Unipotent): PSD cone
- Complex algebras (Pauli, Fermionic, Bosonic): Hermitian PSD cone

# Type Parameters
- `A`: Algebra type determining simplification rules and cone type
- `T`: Integer type for monomial word representation
- `M`: Monomial type for basis elements (may expand for PBW algebras)
- `P`: Polynomial type `Polynomial{A,T,C}` for some coefficient type `C`

# Fields
- `objective::P`: The polynomial objective function
- `constraints::Vector{Tuple{Symbol, Matrix{P}}}`: Constraint matrices with cone types
  - `:Zero` - equality constraint (zeros cone)
  - `:PSD` - real positive semidefinite cone
  - `:HPSD` - Hermitian positive semidefinite cone
- `total_basis::Vector{M}`: Union of all basis monomials across constraints

# Notes
This is a purely symbolic representation with no JuMP model. Use `solve_moment_problem`
to instantiate and solve, or `sos_dualize` to convert to the dual SOS problem.

# Examples
```julia
# After moment_relax, moment_problem is symbolic:
mp = moment_relax(pop, corr_sparsity, cliques_term_sparsities)

# Dualize to SOS form and solve
sos = sos_dualize(mp)
set_optimizer(sos.model, Clarabel.Optimizer)
optimize!(sos.model)
```

See also: [`moment_relax`](@ref), [`sos_dualize`](@ref)
"""
struct MomentProblem{A<:AlgebraType, T<:Integer, M<:NormalMonomial{A,T}, P<:Polynomial{A,T}}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}
    total_basis::Vector{M}
    n_unique_moment_matrix_elements::Int
end


# =============================================================================
# Constraint Matrix Construction (Symbolic)
# =============================================================================

"""
    _build_constraint_matrix(poly::P, local_basis::Vector{M}, cone::Symbol) where {T, P<:AbstractPolynomial{T}, M}

Build a symbolic constraint matrix for the moment relaxation.

# Arguments
- `poly`: The polynomial multiplier (1 for moment matrix, constraint poly for localizing)
- `local_basis`: Vector of Monomial elements indexing rows/columns
- `cone`: Cone type symbol (:Zero, :PSD, or :HPSD)

# Returns
- `Tuple{Symbol, Matrix{P}}`: The cone type and the polynomial-valued constraint matrix

The matrix element at (i,j) is the bilinear expansion:
  M[i,j] = Σ_{k,l} conj(c_ik) * c_jl * Σ_m coef_m * simplify(word_ik† * mono_m * word_jl)

For MonoidAlgebra/TwistedGroupAlgebra, each Monomial has a single term (single word with coefficient).
For PBWAlgebra, Monomials may have multiple terms after normalization.
"""
function _build_constraint_matrix(
    poly::Polynomial{A,T,C},
    local_basis::Vector{M},
    cone::Symbol
) where {A<:AlgebraType,T<:Integer,C<:Number,M<:NormalMonomial{A,T}}
    # Each matrix element is a polynomial: bilinear expansion over basis terms
    moment_mtx = Matrix{Polynomial{A,T,C}}(undef, length(local_basis), length(local_basis))

    for (i, row_mono) in enumerate(local_basis)
        for (j, col_mono) in enumerate(local_basis)
            # Build polynomial for this matrix element with full bilinear expansion
            # For each (c_row, row_word) in row_mono and (c_col, col_word) in col_mono,
            # compute conj(c_row) * c_col * sum over poly terms
            element_poly = sum(
                _conj_coef(A, c_row) * c_col * coef * Polynomial(simplify(A, _neat_dot3(row_word, mono, col_word)))
                for (c_row, row_word) in row_mono
                for (c_col, col_word) in col_mono
                for (coef, mono) in poly.terms
            )
            moment_mtx[i, j] = element_poly
        end
    end

    return (cone, moment_mtx)
end

"""
    _conj_coef(::Type{A}, c) where {A<:AlgebraType}

Conjugate a Monomial coefficient for the bilinear form.

For real algebras (MonoidAlgebra, PBWAlgebra with integer coefficients), returns c unchanged.
For complex algebras (TwistedGroupAlgebra with phase encoding), conjugates the phase.
"""
_conj_coef(::Type{<:AlgebraType}, c) = c  # Default: identity (real coefficients)

# For Pauli algebra, phase is encoded as UInt8: 0=1, 1=i, 2=-1, 3=-i
# conj(i^k) = i^(-k) = i^(4-k mod 4)
_conj_coef(::Type{PauliAlgebra}, c::UInt8) = (0x04 - c) & 0x03


# =============================================================================
# Moment Relaxation (Unified)
# =============================================================================

"""
    moment_relax(pop::PolyOpt{A,TI,P}, corr_sparsity::CorrelativeSparsity, cliques_term_sparsities::Vector{Vector{TermSparsity{M}}})

Construct a symbolic moment relaxation of a polynomial optimization problem.

# Arguments
- `pop::PolyOpt{A,TI,P}`: The polynomial optimization problem
- `corr_sparsity::CorrelativeSparsity`: Correlative sparsity structure with cliques
- `cliques_term_sparsities`: Term sparsity for each clique

# Returns
- `MomentProblem{A,T,M,P}`: Symbolic moment problem ready for dualization or direct solving

# Description
This function creates a symbolic representation of the moment relaxation by:
1. Computing total basis from all clique term sparsities
2. Building constraint matrices as polynomial-valued matrices
3. Selecting cone type based on algebra (PSD for real, HPSD for complex)

The cone type is automatically determined by `_is_complex_problem(A)`:
- Real algebras (NonCommutative, Projector, Unipotent): `:PSD` cone
- Complex algebras (Pauli, Fermionic, Bosonic): `:HPSD` cone

# Examples
```julia
# Build moment relaxation
mp = moment_relax(pop, corr_sparsity, cliques_term_sparsities)

# Can then dualize or solve directly
sos = sos_dualize(mp)
```

See also: [`MomentProblem`](@ref), [`sos_dualize`](@ref)
"""
function moment_relax(
    pop::PolyOpt{A,TI,P},
    corr_sparsity::CorrelativeSparsity{A,TI,P,M,Nothing},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {A<:AlgebraType,TI<:Integer,C<:Number,P<:Polynomial{A,TI,C},M<:NormalMonomial{A,TI}}

    # Unique moment variables are determined by moment matrices only (poly = 1),
    # not by the full set of localizing matrices.
    # For Monomial bases, we expand over all term pairs from the bilinear form.
    one_word = one(NormalMonomial{A,TI})
    moment_matrix_basis = sorted_union(map(cliques_term_sparsities) do term_sparsities
        reduce(vcat, [
            [
                mono
                for row_mono in basis
                for col_mono in basis
                for (_, row_word) in row_mono  # Iterate over terms in row Monomial
                for (_, col_word) in col_mono  # Iterate over terms in col Monomial
                for (_, mono) in Polynomial(simplify(A, _neat_dot3(row_word, one_word, col_word))).terms
            ]
            for basis in term_sparsities[1].block_bases
        ])
    end...)
    n_unique_moment_matrix_elements = length(_sorted_symmetric_basis(moment_matrix_basis))

    # Compute total basis: union of all moment matrix entry monomials
    # Expands over all term pairs from Monomial bases
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        reduce(vcat, [
            [mono
             for m in last.(poly.terms)
             for (_, row_word) in row_mono  # Iterate over terms in row Monomial
             for (_, col_word) in col_mono  # Iterate over terms in col Monomial
             for (_, mono) in Polynomial(simplify(A, _neat_dot3(row_word, m, col_word))).terms]
            for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities)
            for basis in term_sparsity.block_bases
            for row_mono in basis
            for col_mono in basis
        ])
    end...)

    # NOTE: For fermionic algebras, we do NOT filter odd-parity monomials from the basis.
    # The parity superselection rule is enforced via constraints (see _add_parity_constraints!)
    # rather than basis filtering. This is because moment matrix entries M[i,j] = <basis[i]^dag * op * basis[j]>
    # can have even total parity even when basis[i] and basis[j] are individually odd-parity.

    # Determine cone type based on algebra
    is_complex = _is_complex_problem(A)
    psd_cone = is_complex ? :HPSD : :PSD

    # Build constraint matrices symbolically
    constraints = Vector{Tuple{Symbol, Matrix{P}}}()

    # Process clique constraints
    for (term_sparsities, cons_idx) in zip(cliques_term_sparsities, corr_sparsity.clq_cons)
        polys = [one(pop.objective); corr_sparsity.cons[cons_idx]...]

        for (term_sparsity, poly) in zip(term_sparsities, polys)
            for ts_sub_basis in term_sparsity.block_bases
                # Determine cone: Zero for equality constraints, PSD/HPSD otherwise
                cone = poly in pop.eq_constraints ? :Zero : psd_cone
                constraint = _build_constraint_matrix(poly, ts_sub_basis, cone)
                push!(constraints, constraint)
            end
        end
    end

    # Process global constraints
    for global_con in corr_sparsity.global_cons
        poly = corr_sparsity.cons[global_con]
        cone = poly in pop.eq_constraints ? :Zero : psd_cone
        # Global constraints use identity basis (scalar moment)
        constraint = _build_constraint_matrix(poly, [one(M)], cone)
        push!(constraints, constraint)
    end

    mp = MomentProblem{A, TI, M, P}(pop.objective, constraints, total_basis, n_unique_moment_matrix_elements)

    # Add parity superselection constraints for fermionic algebras
    # This enforces that odd-parity moment entries are zero
    _add_parity_constraints!(mp)

    return mp
end


# =============================================================================
# Direct Solving Interface
# =============================================================================

"""
    solve_moment_problem(mp::MomentProblem{A,T,M,P}, optimizer; silent::Bool=true) where {A,T,M,P}

Directly solve a symbolic moment problem by instantiating a JuMP model.

# Arguments
- `mp`: Symbolic moment problem from `moment_relax`
- `optimizer`: JuMP-compatible optimizer (e.g., Clarabel.Optimizer)

# Keyword Arguments
- `silent`: Suppress optimizer output (default: true)

# Returns
- `NamedTuple` with:
  - `objective`: Optimal objective value
  - `model`: The JuMP model (for extracting dual values, etc.)
  - `monomap`: Dictionary mapping monomials to JuMP variables

# Description
This function instantiates the symbolic moment problem as a JuMP model:
1. Creates variables for each monomial in total_basis
2. Sets y[1] = 1 (normalization)
3. Adds constraint matrices in appropriate cones
4. Minimizes the objective
5. Solves and returns results

For most use cases, `sos_dualize` is preferred (smaller SDP, faster solving).
Use this function when you need moment problem solutions directly.

# Examples
```julia
mp = moment_relax(pop, corr_sparsity, term_sparsities)
result = solve_moment_problem(mp, Clarabel.Optimizer)
println("Optimal value: ", result.objective)
```

See also: [`moment_relax`](@ref), [`MomentProblem`](@ref), [`sos_dualize`](@ref)
"""
function solve_moment_problem(
    mp::MomentProblem{A,T,M,P},
    optimizer;
    silent::Bool=true
) where {A<:AlgebraType, T<:Integer, M<:NormalMonomial{A,T}, P<:Polynomial{A,T}}

    # Determine coefficient type from polynomial
    C = eltype(coefficients(mp.objective))

    # For complex algebras, we need real-valued JuMP model
    # The constraint matrices are Hermitian, dualized to real 2x2 blocks
    is_complex = _is_complex_problem(A)

    if is_complex
        return _solve_complex_moment_problem(mp, optimizer, silent)
    else
        return _solve_real_moment_problem(mp, optimizer, silent)
    end
end

"""
    _solve_real_moment_problem(mp, optimizer, silent)

Internal: Solve a real-valued moment problem directly.
"""
function _solve_real_moment_problem(
    mp::MomentProblem{A,T,M,P},
    optimizer,
    silent::Bool
) where {A<:AlgebraType, T<:Integer, M<:NormalMonomial{A,T}, P}

    # Get coefficient type
    C = eltype(coefficients(mp.objective))
    model = GenericModel{C}()

    basis = [symmetric_canon(expval(m)) for m in mp.total_basis]
    sorted_unique!(basis)

    @variable(model, y[1:length(basis)], set_string_name=false)

    one_sym = expval(one(M))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("Expected identity moment to be present in basis")
    @constraint(model, y[idx_one] == 1)  # Normalization

    monomap = Dict(zip(basis, y))

    # Add constraints
    for (cone, mat) in mp.constraints
        dim = size(mat, 1)
        # Convert polynomial matrix to JuMP expression matrix
        jump_mat = [
            _substitute_poly(mat[i,j], monomap)
            for i in 1:dim, j in 1:dim
        ]

        if cone == :Zero
            @constraint(model, jump_mat in Zeros())
        elseif cone == :PSD
            @constraint(model, jump_mat in PSDCone())
        else
            error("Unexpected cone type $cone for real problem")
        end
    end

    # Set objective
    obj_expr = _substitute_poly(mp.objective, monomap)
    @objective(model, Min, obj_expr)

    # Solve
    set_optimizer(model, optimizer)
    silent && set_silent(model)
    optimize!(model)

    n_unique = mp.n_unique_moment_matrix_elements

    return (objective=objective_value(model), model=model, monomap=monomap, n_unique_elements=n_unique)
end

"""
    _solve_complex_moment_problem(mp, optimizer, silent)

Internal: Solve a complex/Hermitian moment problem directly.

For Hermitian PSD constraints, we use the standard embedding:
H in HPSD <==> [Re(H), -Im(H); Im(H), Re(H)] in PSD (real 2n x 2n)
"""
function _solve_complex_moment_problem(
    mp::MomentProblem{A,T,M,P},
    optimizer,
    silent::Bool
) where {A<:AlgebraType, T<:Integer, M<:NormalMonomial{A,T}, P}

    # Get real coefficient type
    C = real(eltype(coefficients(mp.objective)))
    model = GenericModel{C}()

    basis = [symmetric_canon(expval(m)) for m in mp.total_basis]
    sorted_unique!(basis)

    n_basis = length(basis)
    @variable(model, y_re[1:n_basis], set_string_name=false)
    @variable(model, y_im[1:n_basis], set_string_name=false)

    one_sym = expval(one(M))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("Expected identity moment to be present in basis")

    # Normalization: y[1] = 1 (real), y_im[1] = 0
    @constraint(model, y_re[idx_one] == 1)
    @constraint(model, y_im[idx_one] == 0)

    basis_to_idx = Dict(m => i for (i, m) in enumerate(basis))

    # Add constraints with Hermitian embedding
    for (cone, mat) in mp.constraints
        dim = size(mat, 1)

        # Build real and imaginary parts of constraint matrix
        mat_re = Matrix{Any}(undef, dim, dim)
        mat_im = Matrix{Any}(undef, dim, dim)

        for i in 1:dim, j in 1:dim
            re_expr, im_expr = _substitute_complex_poly(mat[i,j], basis_to_idx, y_re, y_im)
            mat_re[i,j] = re_expr
            mat_im[i,j] = im_expr
        end

        if cone == :Zero
            # Zero cone: both real and imaginary parts must be zero
            @constraint(model, [mat_re[i,j] for i in 1:dim, j in 1:dim] .== 0)
            @constraint(model, [mat_im[i,j] for i in 1:dim, j in 1:dim] .== 0)
        elseif cone == :HPSD
            # Hermitian PSD: embed as real 2n x 2n PSD
            # [Re(H), -Im(H); Im(H), Re(H)] in PSD
            embedded = [
                [mat_re[i,j] for i in 1:dim, j in 1:dim]   [-mat_im[i,j] for i in 1:dim, j in 1:dim]
                [mat_im[i,j] for i in 1:dim, j in 1:dim]   [mat_re[i,j] for i in 1:dim, j in 1:dim]
            ]
            @constraint(model, embedded in PSDCone())
        else
            error("Unexpected cone type $cone for complex problem")
        end
    end

    # Set objective (minimize real part, imaginary should be zero for valid problem)
    obj_re, _ = _substitute_complex_poly(mp.objective, basis_to_idx, y_re, y_im)
    @objective(model, Min, obj_re)

    # Solve
    set_optimizer(model, optimizer)
    silent && set_silent(model)
    optimize!(model)

    # Build monomap returning complex values
    monomap = Dict(m => Complex(value(y_re[i]), value(y_im[i])) for (m, i) in basis_to_idx)

    n_unique = mp.n_unique_moment_matrix_elements

    return (objective=objective_value(model), model=model, monomap=monomap, n_unique_elements=n_unique)
end


# =============================================================================
# Helper Functions
# =============================================================================

"""
    _substitute_poly(poly::P, monomap::Dict{M,V}) where {T, P<:AbstractPolynomial{T}, M, V}

Substitute monomials in a polynomial with JuMP variables.
Returns a JuMP affine expression.

Monomials not in monomap are treated as having expectation value 0.
"""
function _substitute_poly(
    poly::P,
    monomap::Dict{M,V}
) where {T, P<:AbstractPolynomial{T}, M, V}
    if iszero(poly)
        return zero(T) * first(values(monomap))
    end
    # Skip monomials not in monomap (expectation value 0)
    return sum(
        (
            canon_mono = symmetric_canon(expval(mono));
            haskey(monomap, canon_mono) ? coef * monomap[canon_mono] : zero(coef) * first(values(monomap))
        )
        for (coef, mono) in zip(coefficients(poly), monomials(poly))
    )
end

"""
    _substitute_complex_poly(poly, basis_to_idx, y_re, y_im)

Substitute monomials in a polynomial with separate real/imaginary JuMP variables.
Returns (real_expr, imag_expr) tuple.

Monomials not in basis_to_idx are treated as having expectation value 0.

# Type Suggestions for mat_re/mat_im
For optimal type stability in complex moment solving, prefer:
- `Matrix{JuMP.GenericAffExpr{Cr,JuMP.VariableRef}}` over `Matrix{Any}`
where `Cr = real(eltype(coefficients(poly)))` (typically Float64).
This ensures type-stable matrix operations during constraint construction.
Currently Matrix{Any} is used for simplicity, but typed matrices would
reduce allocation overhead and enable better JIT optimization.
"""
function _substitute_complex_poly(
    poly::P,
    basis_to_idx::Dict{M,Int},
    y_re::Vector{V},
    y_im::Vector{V}
) where {T, P<:AbstractPolynomial{T}, M, V}

    # For zero polynomial, return correctly typed zero expressions
    # rather than literal (0.0, 0.0) which causes type instability
    if iszero(poly)
        # Use the JuMP variable type to construct proper zero expressions
        zero_re = zero(eltype(y_re)) * y_re[1]
        zero_im = zero(eltype(y_im)) * y_im[1]
        return (zero_re, zero_im)
    end

    # Initialize with properly typed zero expressions
    re_expr = zero(eltype(y_re)) * y_re[1]
    im_expr = zero(eltype(y_im)) * y_im[1]

    for (coef, mono) in zip(coefficients(poly), monomials(poly))
        canon_mono = symmetric_canon(expval(mono))

        # Skip monomials not in basis (expectation value 0)
        if !haskey(basis_to_idx, canon_mono)
            continue
        end

        idx = basis_to_idx[canon_mono]

        # (a + bi)(x + yi) = (ax - by) + (ay + bx)i
        c_re = real(coef)
        c_im = imag(coef)

        re_expr += c_re * y_re[idx] - c_im * y_im[idx]
        im_expr += c_im * y_re[idx] + c_re * y_im[idx]
    end

    return (re_expr, im_expr)
end


# =============================================================================
# Fermionic Parity Superselection Constraints
# =============================================================================

"""
    _has_odd_parity_only(poly::Polynomial{FermionicAlgebra,T,C}) where {T,C}

Check if a polynomial's expectation value must be zero due to parity superselection.
Returns `true` if all non-zero terms have odd parity after canonicalization.

For fermionic systems, only operators with even total fermion parity can have
non-zero expectation values (parity superselection rule). This function checks
whether ALL terms in a polynomial have odd parity, meaning the entire polynomial
must have zero expectation value.

# Arguments
- `poly`: A fermionic polynomial to check

# Returns
- `true` if all non-zero terms have odd parity (expectation value must be 0)
- `false` if at least one term has even parity (expectation value may be non-zero)
"""
function _has_odd_parity_only(
    poly::Polynomial{FermionicAlgebra,T,C}
) where {T<:Integer,C<:Number}
    has_nonzero_term = false

    for t in terms(poly)
        if !iszero(t.coefficient)
            has_nonzero_term = true
            canon_sym = symmetric_canon(expval(t.monomial))
            if has_even_parity(canon_sym.mono)
                return false  # Found even-parity term
            end
        end
    end

    return has_nonzero_term  # true only if all terms are odd parity
end

# Fallback for non-fermionic algebras (always returns false)
_has_odd_parity_only(poly::Polynomial) = false

"""
    _add_parity_constraints!(mp::MomentProblem{A,T,M,P})

Add zero constraints for moment matrix entries with odd fermion parity.

For fermionic algebras, the parity superselection rule requires that expectation
values of odd-parity operators be zero. This function scans all constraint matrices
and adds explicit Zero cone constraints for entries where the polynomial has only
odd-parity terms.

This approach is correct because:
- We keep ALL monomials in the basis (including odd-parity ones)
- Moment matrix entry M[i,j] = <basis[i]^dag * op * basis[j]>
- Even if basis[i] and basis[j] are individually odd-parity, the product
  basis[i]^dag * basis[j] may have even total parity
- We only constrain entries where the TOTAL expression has odd parity

For non-fermionic algebras, this function is a no-op.

# Arguments
- `mp`: A MomentProblem to process (modified in place for fermionic algebras)
"""
function _add_parity_constraints!(
    mp::MomentProblem{A,T,M,P}
) where {A<:AlgebraType,T<:Integer,M,P}
    # Only FermionicAlgebra needs parity constraints
    A === FermionicAlgebra || return nothing

    parity_constraints = Tuple{Symbol, Matrix{P}}[]

    for (cone, mat) in mp.constraints
        dim = size(mat, 1)
        for i in 1:dim, j in 1:dim
            poly = mat[i,j]
            if _has_odd_parity_only(poly)
                # This entry must be zero - add as 1x1 Zero constraint
                push!(parity_constraints, (:Zero, reshape([poly], 1, 1)))
            end
        end
    end

    # Add all parity constraints to the problem
    append!(mp.constraints, parity_constraints)
end

# =============================================================================
# State Moment Problem
# =============================================================================

"""
    StateMomentProblem{A<:AlgebraType, ST<:StateType, T<:Integer, M<:NCStateWord{ST,A,T}, P<:NCStatePolynomial}

A symbolic representation of a state polynomial moment relaxation problem.

Similar to `MomentProblem` but for state polynomial optimization.

# Type Parameters
- `A`: Algebra type
- `ST`: State type (Arbitrary or MaxEntangled)
- `T`: Integer type for word representation
- `M`: NCStateWord type
- `P`: NCStatePolynomial type

# Fields
- `objective::P`: The state polynomial objective function
- `constraints::Vector{Tuple{Symbol, Matrix{P}, Vector{M}}}`: Constraint matrices with cone types and block bases
- `total_basis::Vector{M}`: Union of all basis NCStateWords
"""
struct StateMomentProblem{A<:AlgebraType, ST<:StateType, T<:Integer, M<:NCStateWord{ST,A,T}, P<:NCStatePolynomial}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}, Vector{M}}}  # (cone, matrix, block_basis)
    total_basis::Vector{M}
    n_unique_moment_matrix_elements::Int
end

"""
    _build_state_constraint_matrix(poly, local_basis, cone) -> Tuple{Symbol, Matrix}

Build a symbolic constraint matrix for state polynomial moment relaxation.

# Arguments
- `poly`: The NCStatePolynomial multiplier (1 for moment matrix, constraint poly for localizing)
- `local_basis`: Vector of NCStateWords indexing rows/columns
- `cone`: Cone type symbol (:Zero or :PSD)

# Returns
- `Tuple{Symbol, Matrix{NCStatePolynomial}}`: The cone type and state polynomial-valued matrix
"""
function _build_state_constraint_matrix(
    poly::P,
    local_basis::Vector{M},
    cone::Symbol
) where {ST<:StateType, A<:AlgebraType, T<:Integer, C<:Number, P<:NCStatePolynomial{C,ST,A,T}, M<:NCStateWord{ST,A,T}}
    # Each matrix element is an NCStatePolynomial
    moment_mtx = Matrix{P}(undef, length(local_basis), length(local_basis))

    for (i, row_idx) in enumerate(local_basis)
        for (j, col_idx) in enumerate(local_basis)
            # Build NCStatePolynomial for this matrix element
            # _neat_dot3 returns NCStateWord, simplify to get NCStatePolynomial
            element_poly = zero(P)
            for (coef, ncsw) in zip(coefficients(poly), monomials(poly))
                prod_poly = simplify(_neat_dot3(row_idx, ncsw, col_idx))
                element_poly = element_poly + coef * prod_poly
            end
            moment_mtx[i, j] = element_poly
        end
    end

    return (cone, moment_mtx)
end

"""
    moment_relax(pop::PolyOpt, corr_sparsity, cliques_term_sparsities) -> StateMomentProblem

Construct a symbolic moment relaxation of a state polynomial optimization problem.

# Arguments
- `pop::PolyOpt{A,TI,P}`: The polynomial optimization problem with NCStatePolynomial objective
- `corr_sparsity::CorrelativeSparsity`: Correlative sparsity structure
- `cliques_term_sparsities`: Term sparsity for each clique

# Returns
- `StateMomentProblem{A,ST,T,M,P}`: Symbolic state moment problem
"""
function moment_relax(
    pop::PolyOpt{A,TI,P},
    corr_sparsity::CorrelativeSparsity{A,TI,P,M,ST},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {A<:AlgebraType,TI<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,A,TI},M<:NCStateWord{ST,A,TI}}

    # Unique moment variables are determined by moment matrices only (poly = 1),
    # not by the full set of localizing matrices.
    moment_matrix_basis = sorted_union(map(cliques_term_sparsities) do term_sparsities
        reduce(vcat, [
            [
                ncsw_result
                for row_idx in basis
                for col_idx in basis
                for ncsw_result in monomials(simplify(_neat_dot3(row_idx, one(M), col_idx)))
            ]
            for basis in term_sparsities[1].block_bases
        ])
    end...)
    n_unique_moment_matrix_elements = length(_sorted_stateword_basis_from_ncsw(moment_matrix_basis))

    # Compute total basis: union of all moment matrix entry NCStateWords
    # _neat_dot3 returns NCStateWord, simplify to get NCStatePolynomial
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        reduce(vcat, [
            [ncsw_result
             for ncsw in monomials(poly)
             for basis in term_sparsity.block_bases
             for row_idx in basis
             for col_idx in basis
             for ncsw_result in monomials(simplify(_neat_dot3(row_idx, ncsw, col_idx)))]
            for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities)
        ])
    end...)

    # State polynomial optimization uses real PSD cone (unipotent, projector algebras)
    # These are "real" algebras that don't produce complex phases
    psd_cone = :PSD

    # Build constraint matrices symbolically, storing block basis with each constraint
    constraints = Vector{Tuple{Symbol, Matrix{P}, Vector{M}}}()

    # Process clique constraints
    for (term_sparsities, cons_idx) in zip(cliques_term_sparsities, corr_sparsity.clq_cons)
        polys = [one(pop.objective); corr_sparsity.cons[cons_idx]...]

        for (term_sparsity, poly) in zip(term_sparsities, polys)
            for ts_sub_basis in term_sparsity.block_bases
                # Determine cone: Zero for equality constraints, PSD otherwise
                cone = poly in pop.eq_constraints ? :Zero : psd_cone
                (cone_type, mat) = _build_state_constraint_matrix(poly, ts_sub_basis, cone)
                # Store block basis with constraint for correct coefficient extraction in SOS dualization
                push!(constraints, (cone_type, mat, ts_sub_basis))
            end
        end
    end

    # Process global constraints
    for global_con in corr_sparsity.global_cons
        poly = corr_sparsity.cons[global_con]
        cone = poly in pop.eq_constraints ? :Zero : psd_cone
        # Global constraints use identity basis (scalar moment)
        global_basis = [one(M)]
        (cone_type, mat) = _build_state_constraint_matrix(poly, global_basis, cone)
        push!(constraints, (cone_type, mat, global_basis))
    end

    return StateMomentProblem{A, ST, TI, M, P}(pop.objective, constraints, total_basis, n_unique_moment_matrix_elements)
end


# =============================================================================
# Direct Solving for State Moment Problems
# =============================================================================

"""
    solve_moment_problem(mp::StateMomentProblem{A,ST,T,M,P}, optimizer; silent::Bool=true)

Directly solve a symbolic state moment problem by instantiating a JuMP model.

# Arguments
- `mp`: Symbolic state moment problem from `moment_relax`
- `optimizer`: JuMP-compatible optimizer (e.g., Clarabel.Optimizer)

# Keyword Arguments
- `silent`: Suppress optimizer output (default: true)

# Returns
- `NamedTuple` with:
  - `objective`: Optimal objective value
  - `model`: The JuMP model (for extracting dual values, etc.)
  - `monomap`: Dictionary mapping StateWords to JuMP variable values

# Description
This function instantiates the symbolic state moment problem as a JuMP model:
1. Creates variables for each unique StateWord (via expval of NCStateWord)
2. Sets y[identity] = 1 (normalization)
3. Adds constraint matrices in appropriate cones
4. Minimizes the objective
5. Solves and returns results

State polynomial optimization uses real-valued SDP (no complex embedding needed).

# Examples
```julia
mp = moment_relax(spop, corr_sparsity, term_sparsities)
result = solve_moment_problem(mp, Clarabel.Optimizer)
println("Optimal value: ", result.objective)
```

See also: [`moment_relax`](@ref), [`StateMomentProblem`](@ref), [`sos_dualize`](@ref)
"""
function solve_moment_problem(
    mp::StateMomentProblem{A,ST,T,M,P},
    optimizer;
    silent::Bool=true
) where {A<:AlgebraType, ST<:StateType, T<:Integer, M<:NCStateWord{ST,A,T}, P<:NCStatePolynomial}

    # Get coefficient type from polynomial
    C = eltype(coefficients(mp.objective))

    # State polynomial optimization uses real-valued model
    model = GenericModel{C}()

    # Build the StateWord basis from total_basis NCStateWords
    # We need to convert NCStateWord -> StateWord via expval for moment variables
    SW = StateWord{ST,A,T}
    state_basis = _sorted_stateword_basis_from_ncsw(mp.total_basis)
    n_basis = length(state_basis)

    # Create variables for basis StateWords
    @variable(model, y[1:n_basis], set_string_name=false)

    # Normalization: y[identity] = 1
    identity_sw = one(SW)
    identity_idx = searchsortedfirst(state_basis, identity_sw)
    if identity_idx <= n_basis && state_basis[identity_idx] == identity_sw
        @constraint(model, y[identity_idx] == 1)
    else
        error("Identity StateWord not found in basis - this shouldn't happen")
    end

    # Map StateWord to variable index
    sw_to_idx = Dict(sw => i for (i, sw) in enumerate(state_basis))

    # Add constraints (block_basis not needed here - used only during construction)
    for (cone, mat, _) in mp.constraints
        dim = size(mat, 1)
        # Convert NCStatePolynomial matrix to JuMP expression matrix
        jump_mat = [
            _substitute_state_poly(mat[i,j], sw_to_idx, y)
            for i in 1:dim, j in 1:dim
        ]

        if cone == :Zero
            @constraint(model, jump_mat in Zeros())
        elseif cone == :PSD
            @constraint(model, jump_mat in PSDCone())
        else
            error("Unexpected cone type $cone for state polynomial problem")
        end
    end

    # Set objective
    obj_expr = _substitute_state_poly(mp.objective, sw_to_idx, y)
    @objective(model, Min, obj_expr)

    # Solve
    set_optimizer(model, optimizer)
    silent && set_silent(model)
    optimize!(model)

    # Build monomap returning StateWord -> value
    monomap = Dict(sw => value(y[i]) for (sw, i) in sw_to_idx)

    n_unique = mp.n_unique_moment_matrix_elements

    return (objective=objective_value(model), model=model, monomap=monomap, n_unique_elements=n_unique)
end

"""
    _substitute_state_poly(poly::NCStatePolynomial, sw_to_idx::Dict, y::Vector) -> AffExpr

Substitute StateWords in an NCStatePolynomial with JuMP variables.
Returns a JuMP affine expression.

NCStateWords are converted to StateWords via expval, then canonicalized.
StateWords not in sw_to_idx are treated as having expectation value 0.
"""
function _substitute_state_poly(
    poly::P,
    sw_to_idx::Dict{SW,Int},
    y::Vector{V}
) where {C<:Number, ST<:StateType, A<:AlgebraType, T<:Integer, P<:NCStatePolynomial{C,ST,A,T}, SW<:StateWord{ST,A,T}, V}

    if iszero(poly)
        return zero(eltype(y))
    end

    expr = zero(eltype(y))
    for (coef, ncsw) in zip(coefficients(poly), monomials(poly))
        # Convert NCStateWord to StateWord via expval, then canonicalize
        canon_sw = symmetric_canon(expval(ncsw))

        # Skip StateWords not in basis (expectation value 0)
        if haskey(sw_to_idx, canon_sw)
            idx = sw_to_idx[canon_sw]
            expr += coef * y[idx]
        end
    end

    return expr
end
