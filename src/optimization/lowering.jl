# =============================================================================
# MomentProblem -> JuMP lowering
# =============================================================================

struct Pivot
    key::Any
    constraint_idx::Int
    block_idx::Int
    i::Int
    j::Int
    phase::ComplexF64
    cone::Symbol
end

struct AffineResolver{D,Z}
    values::D
    zero_value::Z
end

struct PivotResolver{P,B,Z}
    pivots::P
    blocks::B
    zero_value::Z
end

function (resolver::AffineResolver)(key)
    return get(resolver.values, key, resolver.zero_value)
end

function (resolver::PivotResolver)(key)
    pivot = resolver.pivots[key]
    return resolver.blocks[pivot.block_idx][pivot.i, pivot.j] / pivot.phase
end

@inline _resolver_zero(resolver::AffineResolver) = copy(resolver.zero_value)
@inline _resolver_zero(resolver::PivotResolver) = copy(resolver.zero_value)

function substitute(poly::Polynomial, resolver)
    expr = _resolver_zero(resolver)
    iszero(poly) && return expr

    for (coef, mono) in poly
        expr += coef * resolver(symmetric_canon(expval(mono)))
    end
    return expr
end

@inline _is_pivot_cone(cone::Symbol) = cone == :PSD || cone == :HPSD

function _strict_unit_phase(coef)
    one_coef = one(coef)
    coef == one_coef && return ComplexF64(1)
    coef == -one_coef && return ComplexF64(-1)
    coef == im * one_coef && return ComplexF64(im)
    coef == -im * one_coef && return ComplexF64(-im)
    return nothing
end

function _pivot_candidate(poly::Polynomial)
    simplified = simplify(poly)
    length(simplified.terms) == 1 || return nothing

    coef, mono = only(simplified.terms)
    phase = _strict_unit_phase(coef)
    phase === nothing && return nothing

    return (symmetric_canon(expval(mono)), phase)
end

function _collect_moment_key!(ordered_keys::Vector{Any}, seen::Set{Any}, key)
    key in seen && return nothing
    push!(seen, key)
    push!(ordered_keys, key)
    return nothing
end

function _collect_polynomial_keys!(ordered_keys::Vector{Any}, seen::Set{Any}, poly::Polynomial)
    simplified = simplify(poly)
    for (coef, mono) in simplified
        iszero(coef) && continue
        _collect_moment_key!(ordered_keys, seen, symmetric_canon(expval(mono)))
    end
    return nothing
end

function _all_moment_keys(mp::MomentProblem)
    ordered_keys = Any[]
    seen = Set{Any}()

    _collect_polynomial_keys!(ordered_keys, seen, mp.objective)
    for (_, mat) in mp.constraints
        for entry in mat
            _collect_polynomial_keys!(ordered_keys, seen, entry)
        end
    end

    return ordered_keys
end

function _discover_pivots_unchecked(mp::MomentProblem)
    pivots = Dict{Any,Pivot}()
    block_idx = 0

    for (constraint_idx, (cone, mat)) in pairs(mp.constraints)
        _is_pivot_cone(cone) || continue
        size(mat, 1) == size(mat, 2) ||
            throw(DimensionMismatch("$cone constraint $constraint_idx must be square, got size $(size(mat))"))

        block_idx += 1
        for i in axes(mat, 1), j in axes(mat, 2)
            candidate = _pivot_candidate(mat[i, j])
            candidate === nothing && continue

            key, phase = candidate
            haskey(pivots, key) && continue  # C1: first lexicographic match wins.
            pivots[key] = Pivot(key, constraint_idx, block_idx, i, j, phase, cone)
        end
    end

    return pivots
end

function orphan_keys(mp::MomentProblem, pivots::Dict{Any,Pivot})
    return [key for key in _all_moment_keys(mp) if !haskey(pivots, key)]
end

function orphan_keys(mp::MomentProblem)
    return orphan_keys(mp, _discover_pivots_unchecked(mp))
end

function _summarize_canonical_keys(keys; limit::Int=5)
    shown = join((sprint(show, key) for key in Iterators.take(keys, limit)), ", ")
    length(keys) > limit && (shown *= ", ...")
    return "[" * shown * "]"
end

"""
    discover_pivots(mp::MomentProblem; orphan_policy=:error) -> Dict{Any,Pivot}

Discover deterministic PSD/HPSD entry pivots for canonical moment keys.
A pivot candidate is exactly a single-term polynomial with coefficient in
`±1, ±im`. The first candidate in `(constraint_idx, i, j)` order wins.
"""
function discover_pivots(mp::MomentProblem; orphan_policy::Symbol=:error)
    orphan_policy in (:error, :aux_psd_free) ||
        throw(ArgumentError("Unsupported orphan_policy $(repr(orphan_policy)); expected :error or :aux_psd_free"))

    pivots = _discover_pivots_unchecked(mp)
    orphans = orphan_keys(mp, pivots)
    if !isempty(orphans)
        if orphan_policy == :error
            throw(ArgumentError(
                "$(length(orphans)) canonical moment(s) have no qualifying PSD/HPSD pivot: " *
                _summarize_canonical_keys(orphans)
            ))
        else
            throw(ArgumentError("orphan_policy=:aux_psd_free is not implemented yet"))
        end
    end

    return pivots
end

"""
    build_jump_model(mp::MomentProblem;
        formulation=:moment_variables,
        representation=:real,
        orphan_policy=:error,
    ) -> (model, extract_monomap)

Lower a symbolic `MomentProblem` into a JuMP model without attaching an
optimizer. The returned `extract_monomap` closure reads the solved model values
and returns the canonical-moment map used by the existing direct moment path.

The default `formulation=:moment_variables, representation=:real` preserves the
historical lowering exactly: real algebras use one free real moment variable per
canonical key, and complex algebras use split real/imaginary moment variables
with Hermitian blocks embedded into real PSD cones.
"""
function build_jump_model(
    mp::MomentProblem{A,T,M,P};
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P<:Polynomial{A,T}}
    formulation in (:moment_variables, :psd_blocks) ||
        throw(ArgumentError("Unsupported formulation $(repr(formulation)); expected :moment_variables or :psd_blocks"))
    representation in (:real, :complex) ||
        throw(ArgumentError("Unsupported representation $(repr(representation)); expected :real or :complex"))
    orphan_policy in (:error, :aux_psd_free) ||
        throw(ArgumentError("Unsupported orphan_policy $(repr(orphan_policy)); expected :error or :aux_psd_free"))

    if formulation == :psd_blocks
        _is_complex_problem(A) ||
            throw(ArgumentError("formulation=:psd_blocks is currently implemented only for complex algebras"))
        representation == :complex ||
            throw(ArgumentError("formulation=:psd_blocks currently supports only representation=:complex"))
        return _build_complex_psd_block_model(mp; orphan_policy=orphan_policy)
    end

    if _is_complex_problem(A)
        representation == :real ||
            throw(ArgumentError("formulation=:moment_variables currently supports representation=:real for complex algebras"))
        return _build_complex_moment_variable_model(mp)
    else
        # Real algebras ignore representation; their moment coordinates are real.
        return _build_real_moment_variable_model(mp)
    end
end

function _moment_variable_basis(mp::MomentProblem{A,T,M,P}) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    basis = [symmetric_canon(expval(m)) for m in mp.total_basis]
    sorted_unique!(basis)
    return basis
end

"""
    _build_real_moment_variable_model(mp) -> (model, extract_monomap)

Historical real-algebra direct moment lowering, factored out of
`solve_moment_problem` without semantic changes.
"""
function _build_real_moment_variable_model(
    mp::MomentProblem{A,T,M,P},
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    C = eltype(coefficients(mp.objective))
    model = GenericModel{C}()

    basis = _moment_variable_basis(mp)

    @variable(model, y[1:length(basis)], set_string_name=false)

    one_sym = symmetric_canon(expval(one(M)))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("Expected identity moment to be present in basis")
    @constraint(model, y[idx_one] == 1)  # Normalization

    monovars = Dict(zip(basis, y))
    resolver_values = Dict(m => one(C) * monovars[m] for m in basis)
    resolver = AffineResolver(resolver_values, zero(C) * y[1])

    # Add constraints.
    for (cone, mat) in mp.constraints
        dim = size(mat, 1)
        jump_mat = [
            substitute(mat[i, j], resolver)
            for i in 1:dim, j in 1:dim
        ]

        if cone == :Zero
            @constraint(model, jump_mat in Zeros())
        elseif cone == :PSD
            @constraint(model, _checked_symmetric(jump_mat; context="real PSD constraint") in PSDCone())
        else
            error("Unexpected cone type $cone for real problem")
        end
    end

    obj_expr = substitute(mp.objective, resolver)
    @objective(model, Min, obj_expr)

    extract_monomap = function ()
        return Dict(m => value(monovars[m]) for m in basis)
    end

    return model, extract_monomap
end

function _solve_real_moment_problem(
    mp::MomentProblem{A,T,M,P},
    optimizer,
    silent::Bool,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    model, extract_monomap = _build_real_moment_variable_model(mp)
    set_optimizer(model, optimizer)
    silent && set_silent(model)
    optimize!(model)
    return (
        objective=objective_value(model),
        model=model,
        monomap=extract_monomap(),
        n_unique_elements=mp.n_unique_moment_matrix_elements,
    )
end

"""
    _build_complex_moment_variable_model(mp) -> (model, extract_monomap)

Historical complex/Hermitian direct moment lowering, factored out of
`solve_moment_problem` without semantic changes.
"""
function _build_complex_moment_variable_model(
    mp::MomentProblem{A,T,M,P},
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    C = real(eltype(coefficients(mp.objective)))
    model = GenericModel{C}()

    basis = _moment_variable_basis(mp)
    n_basis = length(basis)

    @variable(model, y_re[1:n_basis], set_string_name=false)
    @variable(model, y_im[1:n_basis], set_string_name=false)

    one_sym = symmetric_canon(expval(one(M)))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("Expected identity moment to be present in basis")

    # Normalization: y[1] = 1 (real), y_im[1] = 0.
    @constraint(model, y_re[idx_one] == 1)
    @constraint(model, y_im[idx_one] == 0)

    basis_to_idx = Dict(m => i for (i, m) in enumerate(basis))
    zero_complex_moment = zero(C) * y_re[1] + im * (zero(C) * y_im[1])
    resolver_values = Dict(m => y_re[i] + im * y_im[i] for (m, i) in basis_to_idx)
    resolver = AffineResolver(resolver_values, zero_complex_moment)

    # Add constraints with Hermitian embedding.
    for (cone, mat) in mp.constraints
        dim = size(mat, 1)

        zero_expr = substitute(zero(eltype(mat)), resolver)
        Aff = typeof(real(zero_expr))
        mat_re = Matrix{Aff}(undef, dim, dim)
        mat_im = Matrix{Aff}(undef, dim, dim)

        for i in 1:dim, j in 1:dim
            expr = substitute(mat[i, j], resolver)
            mat_re[i, j] = real(expr)
            mat_im[i, j] = imag(expr)
        end

        if cone == :Zero
            # Zero cone: both real and imaginary parts must be zero.
            @constraint(model, [mat_re[i, j] for i in 1:dim, j in 1:dim] .== 0)
            @constraint(model, [mat_im[i, j] for i in 1:dim, j in 1:dim] .== 0)
        elseif cone == :HPSD
            # Hermitian PSD: embed as real 2n x 2n PSD.
            # [Re(H), -Im(H); Im(H), Re(H)] in PSD.
            embedded = [
                [mat_re[i, j] for i in 1:dim, j in 1:dim]   [-mat_im[i, j] for i in 1:dim, j in 1:dim]
                [mat_im[i, j] for i in 1:dim, j in 1:dim]   [mat_re[i, j] for i in 1:dim, j in 1:dim]
            ]
            @constraint(model, embedded in PSDCone())
        else
            error("Unexpected cone type $cone for complex problem")
        end
    end

    # `polyopt` validates that complex-algebra objectives are Hermitian, so their
    # expectation values are real. Optimize the real part only.
    obj_expr = substitute(mp.objective, resolver)
    @objective(model, Min, real(obj_expr))

    extract_monomap = function ()
        return Dict(m => Complex(value(y_re[i]), value(y_im[i])) for (m, i) in basis_to_idx)
    end

    return model, extract_monomap
end

function _solve_complex_moment_problem(
    mp::MomentProblem{A,T,M,P},
    optimizer,
    silent::Bool,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    model, extract_monomap = _build_complex_moment_variable_model(mp)
    set_optimizer(model, optimizer)
    silent && set_silent(model)
    optimize!(model)
    return (
        objective=objective_value(model),
        model=model,
        monomap=extract_monomap(),
        n_unique_elements=mp.n_unique_moment_matrix_elements,
    )
end

function _declare_psd_block_variable!(model, cone::Symbol, n::Int)
    if cone == :HPSD
        return @variable(model, [1:n, 1:n] in HermitianPSDCone(), set_string_name=false)
    elseif cone == :PSD
        return @variable(model, [1:n, 1:n] in PSDCone(), set_string_name=false)
    else
        error("Unexpected PSD-block cone $cone")
    end
end

function _declare_psd_block_variables!(model, mp::MomentProblem)
    blocks = Dict{Int,Any}()
    constraint_to_block = Dict{Int,Int}()
    block_idx = 0

    for (constraint_idx, (cone, mat)) in pairs(mp.constraints)
        _is_pivot_cone(cone) || continue
        size(mat, 1) == size(mat, 2) ||
            throw(DimensionMismatch("$cone constraint $constraint_idx must be square, got size $(size(mat))"))

        block_idx += 1
        constraint_to_block[constraint_idx] = block_idx
        blocks[block_idx] = _declare_psd_block_variable!(model, cone, size(mat, 1))
    end

    isempty(blocks) && throw(ArgumentError("formulation=:psd_blocks requires at least one PSD/HPSD constraint"))
    return blocks, constraint_to_block
end

function _pivot_entry_lookup(pivots::Dict{Any,Pivot})
    lookup = Dict{Tuple{Int,Int,Int},Pivot}()
    for pivot in values(pivots)
        lookup[(pivot.constraint_idx, pivot.i, pivot.j)] = pivot
    end
    return lookup
end

@inline _upper_triangle_indices(n::Int) = ((i, j) for j in 1:n for i in 1:j)

function _add_psd_block_bindings!(model, X, mat, resolver, pivot_entries, constraint_idx::Int)
    n = size(mat, 1)
    for (i, j) in _upper_triangle_indices(n)
        haskey(pivot_entries, (constraint_idx, i, j)) && continue
        @constraint(model, X[i, j] == substitute(mat[i, j], resolver))
    end
    return nothing
end

function _add_zero_bindings!(model, mat, resolver, constraint_idx::Int)
    size(mat, 1) == size(mat, 2) ||
        throw(DimensionMismatch("Zero constraint $constraint_idx must be square, got size $(size(mat))"))
    _is_hermitian_poly_matrix(mat) ||
        throw(ArgumentError("Zero constraint $constraint_idx is not Hermitian; MomentProblem zero matrices must be split before lowering"))

    n = size(mat, 1)
    for (i, j) in _upper_triangle_indices(n)
        iszero(mat[i, j]) && continue
        @constraint(model, substitute(mat[i, j], resolver) == 0)
    end
    return nothing
end

function _build_complex_psd_block_model(
    mp::MomentProblem{A,T,M,P};
    orphan_policy::Symbol=:error,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    C = real(eltype(coefficients(mp.objective)))
    model = GenericModel{C}()

    pivots = discover_pivots(mp; orphan_policy=orphan_policy)
    blocks, constraint_to_block = _declare_psd_block_variables!(model, mp)

    anchor_block = blocks[first(sort!(collect(keys(blocks))))]
    anchor = anchor_block[1, 1]
    zero_complex_moment = zero(C) * real(anchor) + im * (zero(C) * real(anchor))
    resolver = PivotResolver(pivots, blocks, zero_complex_moment)

    one_key = symmetric_canon(expval(one(M)))
    haskey(pivots, one_key) ||
        throw(ArgumentError("Identity moment has no qualifying PSD/HPSD pivot"))
    @constraint(model, resolver(one_key) == 1)

    pivot_entries = _pivot_entry_lookup(pivots)
    for (constraint_idx, (cone, mat)) in pairs(mp.constraints)
        if _is_pivot_cone(cone)
            block_idx = constraint_to_block[constraint_idx]
            _add_psd_block_bindings!(model, blocks[block_idx], mat, resolver, pivot_entries, constraint_idx)
        elseif cone == :Zero
            _add_zero_bindings!(model, mat, resolver, constraint_idx)
        else
            error("Unexpected cone type $cone for complex PSD-block problem")
        end
    end

    obj_expr = substitute(mp.objective, resolver)
    @objective(model, Min, real(obj_expr))

    extract_monomap = function ()
        return Dict(key => value(resolver(key)) for key in keys(pivots))
    end

    return model, extract_monomap
end
