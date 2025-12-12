# =============================================================================
# Unified Symbolic MomentProblem
# =============================================================================

"""
    MomentProblem{A<:AlgebraType, T<:Integer, M<:Monomial{A,T}, P<:Polynomial{A,T}}

A symbolic representation of a moment relaxation problem.

This unified type handles both real and complex (Hermitian) moment problems.
The algebra type `A` determines which cone type to use when solving:
- Real algebras (NonCommutative, Projector, Unipotent): PSD cone
- Complex algebras (Pauli, Fermionic, Bosonic): Hermitian PSD cone

# Type Parameters
- `A`: Algebra type determining simplification rules and cone type
- `T`: Integer type for monomial word representation
- `M`: Monomial type `Monomial{A,T}`
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
struct MomentProblem{A<:AlgebraType, T<:Integer, M<:Monomial{A,T}, P<:Polynomial{A,T}}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}
    total_basis::Vector{M}
end


# =============================================================================
# Constraint Matrix Construction (Symbolic)
# =============================================================================

"""
    _build_constraint_matrix(poly::P, local_basis::Vector{M}, cone::Symbol) where {T, P<:AbstractPolynomial{T}, M}

Build a symbolic constraint matrix for the moment relaxation.

# Arguments
- `poly`: The polynomial multiplier (1 for moment matrix, constraint poly for localizing)
- `local_basis`: Vector of monomials indexing rows/columns
- `cone`: Cone type symbol (:Zero, :PSD, or :HPSD)

# Returns
- `Tuple{Symbol, Matrix{P}}`: The cone type and the polynomial-valued constraint matrix

The matrix element at (i,j) is sum over terms: coef * basis[i]' * mono * basis[j]
"""
function _build_constraint_matrix(
    poly::P,
    local_basis::Vector{M},
    cone::Symbol
) where {T, P<:AbstractPolynomial{T}, M}
    # Each matrix element is a polynomial: sum of coef * neat_dot(row, mono, col)
    moment_mtx = Matrix{P}(undef, length(local_basis), length(local_basis))

    for (i, row_idx) in enumerate(local_basis)
        for (j, col_idx) in enumerate(local_basis)
            # Build polynomial for this matrix element
            element_poly = sum(
                coef * _neat_dot3(row_idx, mono, col_idx)
                for (coef, mono) in zip(coefficients(poly), monomials(poly))
            )
            moment_mtx[i, j] = element_poly
        end
    end

    return (cone, moment_mtx)
end


# =============================================================================
# Moment Relaxation (Unified)
# =============================================================================

"""
    moment_relax(pop::PolyOpt{A,P}, corr_sparsity::CorrelativeSparsity, cliques_term_sparsities::Vector{Vector{TermSparsity{M}}})

Construct a symbolic moment relaxation of a polynomial optimization problem.

# Arguments
- `pop::PolyOpt{A,P}`: The polynomial optimization problem
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
    pop::PolyOpt{A,P},
    corr_sparsity::CorrelativeSparsity{A,TI,P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {A<:AlgebraType, TI<:Integer, C<:Number, P<:Polynomial{A,TI,C}, M<:Monomial{A,TI}}

    # Compute total basis: union of all moment matrix entry monomials
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        reduce(vcat, [
            map(monomials(poly)) do m
                expval(_neat_dot3(rol_idx, m, col_idx))
            end
            for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities)
            for basis in term_sparsity.block_bases
            for rol_idx in basis
            for col_idx in basis
        ])
    end...)

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

    return MomentProblem{A, TI, M, P}(pop.objective, constraints, total_basis)
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
) where {A<:AlgebraType, T<:Integer, M<:Monomial{A,T}, P<:Polynomial{A,T}}

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
) where {A<:AlgebraType, T<:Integer, M<:Monomial{A,T}, P}

    # Get coefficient type
    C = eltype(coefficients(mp.objective))
    model = GenericModel{C}()

    # Create variables for basis monomials
    @variable(model, y[1:length(mp.total_basis)], set_string_name=false)
    @constraint(model, first(y) == 1)  # Normalization

    monomap = Dict(zip(mp.total_basis, y))

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

    return (objective=objective_value(model), model=model, monomap=monomap)
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
) where {A<:AlgebraType, T<:Integer, M<:Monomial{A,T}, P}

    # Get real coefficient type
    C = real(eltype(coefficients(mp.objective)))
    model = GenericModel{C}()

    # For complex problems, we need separate real/imag variables
    n_basis = length(mp.total_basis)
    @variable(model, y_re[1:n_basis], set_string_name=false)
    @variable(model, y_im[1:n_basis], set_string_name=false)

    # Normalization: y[1] = 1 (real), y_im[1] = 0
    @constraint(model, y_re[1] == 1)
    @constraint(model, y_im[1] == 0)

    # Map monomials to complex JuMP expressions
    # We'll create a map that returns (re_expr, im_expr) tuples
    basis_to_idx = Dict(m => i for (i, m) in enumerate(mp.total_basis))

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
    monomap = Dict(
        m => Complex(value(y_re[i]), value(y_im[i]))
        for (m, i) in basis_to_idx
    )

    return (objective=objective_value(model), model=model, monomap=monomap)
end


# =============================================================================
# Helper Functions
# =============================================================================

"""
    _substitute_poly(poly::P, monomap::Dict{M,V}) where {T, P<:AbstractPolynomial{T}, M, V}

Substitute monomials in a polynomial with JuMP variables.
Returns a JuMP affine expression.
"""
function _substitute_poly(
    poly::P,
    monomap::Dict{M,V}
) where {T, P<:AbstractPolynomial{T}, M, V}
    if iszero(poly)
        return zero(T) * first(values(monomap))
    end
    return sum(
        coef * monomap[symmetric_canon(expval(mono))]
        for (coef, mono) in zip(coefficients(poly), monomials(poly))
    )
end

"""
    _substitute_complex_poly(poly, basis_to_idx, y_re, y_im)

Substitute monomials in a polynomial with separate real/imaginary JuMP variables.
Returns (real_expr, imag_expr) tuple.
"""
function _substitute_complex_poly(
    poly::P,
    basis_to_idx::Dict{M,Int},
    y_re::Vector{V},
    y_im::Vector{V}
) where {T, P<:AbstractPolynomial{T}, M, V}

    if iszero(poly)
        return (0.0, 0.0)
    end

    re_expr = zero(eltype(y_re))
    im_expr = zero(eltype(y_im))

    for (coef, mono) in zip(coefficients(poly), monomials(poly))
        canon_mono = symmetric_canon(expval(mono))
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
# Backward Compatibility
# =============================================================================

# Legacy type alias for gradual migration
# NOTE: The old MomentProblem had different type parameters:
# MomentProblem{T,M,CR<:ConstraintRef,JS<:AbstractJuMPScalar}
# This alias helps during transition but callers may need updates.

"""
    complex_moment_relax(pop, corr_sparsity, cliques_term_sparsities)

Backward compatibility alias for `moment_relax`.

In the new unified design, `moment_relax` handles both real and complex algebras
based on `_is_complex_problem(A)`. This function is provided for compatibility
with existing code that explicitly called `complex_moment_relax`.

# Deprecated
Use `moment_relax` directly. This alias will be removed in Phase 4.

See also: [`moment_relax`](@ref)
"""
function complex_moment_relax(
    pop::PolyOpt{A,P},
    corr_sparsity::CorrelativeSparsity{A,TI,P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {A<:AlgebraType, TI<:Integer, C<:Number, P<:Polynomial{A,TI,C}, M<:Monomial{A,TI}}
    # Delegate to unified moment_relax
    return moment_relax(pop, corr_sparsity, cliques_term_sparsities)
end
