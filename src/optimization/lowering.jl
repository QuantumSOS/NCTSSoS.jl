# =============================================================================
# MomentProblem -> JuMP lowering
# =============================================================================

const AUX_ORPHANS_PER_BLOCK = 32

struct AffineResolver{D,Z}
    values::D
    zero_value::Z
end

struct BlockMomentResolver{L,B,F,A,Z}
    linear::L
    blocks::B
    free_values::F
    aux_pivots::A
    zero_value::Z
end

function _dict_get_value_or(dict::Dict, key, default)
    haskey(dict, key) && return dict[key]

    # Current canonical keys may be vectors. Keep lookup value-based so lowering
    # is not hostage to object identity if a form carries an equal-but-distinct key.
    for (candidate, value) in dict
        key_isequal(candidate, key) && return value
    end
    return default
end

function (resolver::AffineResolver)(key)
    return _dict_get_value_or(resolver.values, key, resolver.zero_value)
end

function _pivot_value(pivot::Pivot, blocks)
    X = blocks[pivot.block]
    if pivot.adjoint
        # The pivot is stored at the Hermitian upper-triangle position for α*, so
        # ⟨α*⟩ = phase * X[j,i], not conj(phase) * X[i,j].
        return pivot.phase * X[pivot.col, pivot.row]
    else
        # X[i,j] = phase * ⟨α⟩ for unit phases.
        return conj(pivot.phase) * X[pivot.row, pivot.col]
    end
end

function (resolver::BlockMomentResolver)(key)
    pivot = _dict_get_value_or(resolver.linear.pivots, key, nothing)
    pivot !== nothing && return _pivot_value(pivot, resolver.blocks)

    aux_pivot = _dict_get_value_or(resolver.aux_pivots, key, nothing)
    aux_pivot !== nothing && return _pivot_value(aux_pivot, resolver.blocks)

    free_value = _dict_get_value_or(resolver.free_values, key, nothing)
    free_value !== nothing && return free_value

    return resolver.zero_value
end

@inline _resolver_zero(resolver::AffineResolver) = copy(resolver.zero_value)
@inline _resolver_zero(resolver::BlockMomentResolver) = copy(resolver.zero_value)

function _eval_form(form::LinearMomentForm, resolver)
    expr = _resolver_zero(resolver)
    for (key, coef) in form
        add_to_expression!(expr, coef, resolver(key))
    end
    return expr
end

# Kept for state-polynomial helpers in moment.jl. New MomentProblem lowering uses
# the cached linear forms in `mp.linear` instead of substituting polynomials.
function substitute(poly::Polynomial, resolver)
    expr = _resolver_zero(resolver)
    iszero(poly) && return expr

    for (coef, mono) in simplify(poly)
        iszero(coef) && continue
        add_to_expression!(expr, coef, resolver(symmetric_canon(expval(mono))))
    end
    return expr
end

@inline _is_pivot_cone(cone::Symbol) = cone == :PSD || cone == :HPSD

function _summarize_canonical_keys(keys; limit::Int=5)
    shown = join((sprint(show, key) for key in Iterators.take(keys, limit)), ", ")
    length(keys) > limit && (shown *= ", ...")
    return "[" * shown * "]"
end

function _n_aux_blocks(n_orphans::Int; orphans_per_block::Int=AUX_ORPHANS_PER_BLOCK)
    n_orphans == 0 && return 0
    return cld(n_orphans, orphans_per_block)
end

function orphan_keys(mp::MomentProblem)
    return copy(mp.linear.free_keys)
end

"""
    build_jump_model(mp::MomentProblem;
        formulation=:moment_variables,
        representation=:real,
        orphan_policy=:error,
    ) -> (model, extract_monomap)

Lower a symbolic `MomentProblem` into a JuMP model without attaching an
optimizer. The public API is unchanged; the implementation consumes the cached
`mp.linear` view instead of rediscovering moment pivots and key universes.
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
    orphan_policy in (:error, :free_variables, :aux_psd_free) ||
        throw(ArgumentError("Unsupported orphan_policy $(repr(orphan_policy)); expected :error, :free_variables, or :aux_psd_free"))

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
        return _build_real_moment_variable_model(mp)
    end
end

function _lowering_coeff_type(::MomentLinearData{K,C,M}) where {K,C,M}
    return C
end

function _lowering_real_type(::MomentLinearData{K,C,M}) where {K,C,M}
    return typeof(real(one(C)))
end

function _check_real_lowering_cones!(mp::MomentProblem)
    for (cone, _) in mp.constraints
        (cone == :Zero || cone == :PSD) || error("Unexpected cone type $cone for real problem")
    end
    return nothing
end

function _check_complex_moment_variable_cones!(mp::MomentProblem)
    for (cone, _) in mp.constraints
        (cone == :Zero || cone == :HPSD) || error("Unexpected cone type $cone for complex problem")
    end
    return nothing
end

function _check_complex_psd_block_cones!(mp::MomentProblem)
    for (cone, _) in mp.constraints
        (cone == :Zero || _is_pivot_cone(cone)) ||
            error("Unexpected cone type $cone for complex PSD-block problem")
    end
    return nothing
end

function _moment_values_resolver(model, L::MomentLinearData{K}, ::Type{C}) where {K,C}
    @variable(model, y[1:length(L.moments)], set_string_name=false)
    values = Dict{K,Any}(key => one(C) * y[idx] for (idx, key) in enumerate(L.moments))
    return y, AffineResolver(values, zero(C) * y[1])
end

function _complex_moment_values_resolver(model, L::MomentLinearData{K}, ::Type{C}) where {K,C}
    @variable(model, y_re[1:length(L.moments)], set_string_name=false)
    @variable(model, y_im[1:length(L.moments)], set_string_name=false)
    values = Dict{K,Any}(key => y_re[idx] + im * y_im[idx] for (idx, key) in enumerate(L.moments))
    zero_complex = zero(C) * y_re[1] + im * (zero(C) * y_im[1])
    return y_re, y_im, AffineResolver(values, zero_complex)
end

function _identity_index(L::MomentLinearData)
    return _get_key_value(L.moment_index, L.identity, "identity moment index")
end

"""
    _build_real_moment_variable_model(mp) -> (model, extract_monomap)

Historical real-algebra direct moment lowering, now driven by `mp.linear`.
"""
function _build_real_moment_variable_model(
    mp::MomentProblem{A,T,M,P},
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    _check_real_lowering_cones!(mp)

    L = mp.linear
    C = _lowering_real_type(L)
    model = GenericModel{C}()
    y, resolver = _moment_values_resolver(model, L, C)

    @constraint(model, y[_identity_index(L)] == one(C))

    for block in L.psd_blocks_lin
        block.meta.cone == :PSD || error("Unexpected cone type $(block.meta.cone) for real problem")
        n = block.size
        jump_mat = [_eval_form(block.entries[i, j], resolver) for i in 1:n, j in 1:n]
        @constraint(model, _checked_symmetric(jump_mat; context="real PSD constraint") in PSDCone())
    end

    for zc in L.zero_constraints
        @constraint(model, _eval_form(zc.form, resolver) == zero(C))
    end

    @objective(model, Min, _eval_form(L.objective_lin, resolver))

    extract_monomap = function ()
        return Dict(key => value(y[idx]) for (idx, key) in enumerate(L.moments))
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

Historical complex/Hermitian direct moment lowering, now driven by `mp.linear`.
"""
function _build_complex_moment_variable_model(
    mp::MomentProblem{A,T,M,P},
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    _check_complex_moment_variable_cones!(mp)

    L = mp.linear
    C = _lowering_real_type(L)
    model = GenericModel{C}()
    y_re, y_im, resolver = _complex_moment_values_resolver(model, L, C)

    idx_one = _identity_index(L)
    @constraint(model, y_re[idx_one] == one(C))
    @constraint(model, y_im[idx_one] == zero(C))

    for block in L.psd_blocks_lin
        block.meta.cone == :HPSD || error("Unexpected cone type $(block.meta.cone) for complex problem")
        n = block.size

        zero_expr = _resolver_zero(resolver)
        Aff = typeof(real(zero_expr))
        mat_re = Matrix{Aff}(undef, n, n)
        mat_im = Matrix{Aff}(undef, n, n)

        for i in 1:n, j in 1:n
            expr = _eval_form(block.entries[i, j], resolver)
            mat_re[i, j] = real(expr)
            mat_im[i, j] = imag(expr)
        end

        embedded = [
            [mat_re[i, j] for i in 1:n, j in 1:n]   [-mat_im[i, j] for i in 1:n, j in 1:n]
            [mat_im[i, j] for i in 1:n, j in 1:n]   [mat_re[i, j] for i in 1:n, j in 1:n]
        ]
        @constraint(model, embedded in PSDCone())
    end

    for zc in L.zero_constraints
        @constraint(model, real(_eval_form(zc.form, resolver)) == zero(C))
    end

    obj_expr = _eval_form(L.objective_lin, resolver)
    @objective(model, Min, real(obj_expr))

    extract_monomap = function ()
        return Dict(key => Complex(value(y_re[idx]), value(y_im[idx])) for (idx, key) in enumerate(L.moments))
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
    if cone == :HPSD || cone == :AuxHPSD
        return @variable(model, [1:n, 1:n] in HermitianPSDCone(), set_string_name=false)
    elseif cone == :PSD
        return @variable(model, [1:n, 1:n] in PSDCone(), set_string_name=false)
    else
        error("Unexpected PSD-block cone $cone")
    end
end

function _declare_psd_block_variables!(model, L::MomentLinearData)
    isempty(L.psd_blocks_lin) && throw(ArgumentError("formulation=:psd_blocks requires at least one PSD/HPSD constraint"))
    return Any[_declare_psd_block_variable!(model, block.meta.cone, block.size) for block in L.psd_blocks_lin]
end

function _declare_aux_orphan_blocks!(
    model,
    blocks::Vector{Any},
    n_orphans::Int;
    orphans_per_block::Int=AUX_ORPHANS_PER_BLOCK,
)
    n_orphans == 0 && return 0

    n_aux = _n_aux_blocks(n_orphans; orphans_per_block=orphans_per_block)
    for aux_offset in 0:(n_aux - 1)
        n_in_block = min(orphans_per_block, n_orphans - aux_offset * orphans_per_block)
        push!(blocks, _declare_psd_block_variable!(model, :AuxHPSD, n_in_block + 1))
    end

    return n_aux
end

function _declare_free_orphan_moments!(model, orphan_keys::AbstractVector)
    free_values = Dict{eltype(orphan_keys),Any}()
    isempty(orphan_keys) && return free_values

    n_orphans = length(orphan_keys)
    @variable(model, orphan_re[1:n_orphans], set_string_name=false)
    @variable(model, orphan_im[1:n_orphans], set_string_name=false)
    for (idx, key) in enumerate(orphan_keys)
        free_values[key] = orphan_re[idx] + im * orphan_im[idx]
    end
    return free_values
end

function _aux_orphan_pivots(
    ::Type{C},
    orphan_keys::AbstractVector{K},
    first_aux_block_idx::Int;
    orphans_per_block::Int=AUX_ORPHANS_PER_BLOCK,
) where {C,K}
    aux_pivots = Dict{K,Pivot{C}}()
    for (offset, key) in enumerate(orphan_keys)
        local_idx = mod1(offset, orphans_per_block)
        block_offset = (offset - 1) ÷ orphans_per_block
        block_idx = first_aux_block_idx + block_offset
        aux_pivots[key] = Pivot{C}(block_idx, 1, 1 + local_idx, one(C), false)
    end
    return aux_pivots
end

@inline _upper_triangle_indices(n::Int) = ((i, j) for j in 1:n for i in 1:j)

function _add_psd_block_bindings!(model, X, block::PSDBlockLin, resolver, block_idx::Int)
    for (i, j) in _upper_triangle_indices(block.size)
        haskey(resolver.linear.pivot_at, (block_idx, i, j)) && continue
        @constraint(model, X[i, j] == _eval_form(block.entries[i, j], resolver))
    end
    return nothing
end

function _build_complex_psd_block_model(
    mp::MomentProblem{A,T,M,P};
    orphan_policy::Symbol=:error,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    _check_complex_psd_block_cones!(mp)

    L = mp.linear
    C = _lowering_real_type(L)
    LC = _lowering_coeff_type(L)
    model = GenericModel{C}()

    if !isempty(L.free_keys) && orphan_policy == :error
        throw(ArgumentError(
            "$(length(L.free_keys)) canonical moment(s) have no qualifying PSD/HPSD pivot: " *
            _summarize_canonical_keys(L.free_keys)
        ))
    end

    blocks = _declare_psd_block_variables!(model, L)
    free_values = Dict{eltype(L.free_keys),Any}()
    aux_pivots = Dict{eltype(L.free_keys),Pivot{LC}}()

    if !isempty(L.free_keys)
        if orphan_policy == :aux_psd_free
            first_aux_block_idx = length(blocks) + 1
            _declare_aux_orphan_blocks!(model, blocks, length(L.free_keys))
            aux_pivots = _aux_orphan_pivots(LC, L.free_keys, first_aux_block_idx)
        elseif orphan_policy == :free_variables
            free_values = _declare_free_orphan_moments!(model, L.free_keys)
        end
    end

    anchor = blocks[1][1, 1]
    zero_complex_moment = zero(C) * real(anchor) + im * (zero(C) * real(anchor))
    resolver = BlockMomentResolver(L, blocks, free_values, aux_pivots, zero_complex_moment)

    @constraint(model, resolver(L.identity) == one(C))

    for (block_idx, block) in enumerate(L.psd_blocks_lin)
        _add_psd_block_bindings!(model, blocks[block_idx], block, resolver, block_idx)
    end

    for zc in L.zero_constraints
        @constraint(model, real(_eval_form(zc.form, resolver)) == zero(C))
    end

    obj_expr = _eval_form(L.objective_lin, resolver)
    @objective(model, Min, real(obj_expr))

    extract_monomap = function ()
        return Dict(key => value(resolver(key)) for key in L.moments)
    end

    return model, extract_monomap
end
