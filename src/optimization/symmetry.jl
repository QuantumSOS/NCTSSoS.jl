# =============================================================================
# Symmetry Reduction for Dense Monoid Moment Relaxations
# =============================================================================

import SymbolicWedderburn
import GroupsCore

const _SYMMETRY_ATOL = 1e-10

"""
    SignedPermutation

Finite signed permutation acting on registry indices.

Stored images are `(sign, target_index)` pairs with `sign ∈ {-1, +1}`.
Indices not listed in `images` are fixed.

This is the only symmetry action supported by the current MVP reduction path.
"""
struct SignedPermutation{T<:Integer}
    images::Dict{T,Tuple{Int,T}}

    function SignedPermutation{T}(images::Dict{T,Tuple{Int,T}}) where {T<:Integer}
        normalized = Dict{T,Tuple{Int,T}}()
        for (src, (sign, dst)) in images
            sign ∈ (-1, 1) || throw(ArgumentError(
                "SignedPermutation signs must be ±1; got $sign for source index $src."
            ))
            if sign != 1 || dst != src
                normalized[src] = (sign, dst)
            end
        end
        return new{T}(normalized)
    end
end

function SignedPermutation(images::AbstractDict{T,<:Tuple{<:Integer,T}}) where {T<:Integer}
    normalized = Dict{T,Tuple{Int,T}}()
    for (src, (sign, dst)) in images
        normalized[src] = (Int(sign), dst)
    end
    return SignedPermutation{T}(normalized)
end

function _signed_permutation_require_integer(value, label::AbstractString; src=nothing)
    value isa Integer && return value
    context = isnothing(src) ? "" : " for source index $src"
    throw(ArgumentError(
        "SignedPermutation $label must be an integer$context; got $(repr(value)) of type $(typeof(value))."
    ))
end

function _parse_signed_permutation_pair(pair::Pair)
    src = _signed_permutation_require_integer(pair.first, "source index")
    value = pair.second

    if value isa Tuple
        length(value) == 2 || throw(ArgumentError(
            "SignedPermutation tuple images must be `(sign, target_index)`; got $(repr(value)) for source index $src."
        ))

        sign_raw, dst_raw = value
        sign_raw isa Integer || throw(ArgumentError(
            "SignedPermutation signs must be integer ±1; got $(repr(sign_raw)) for source index $src."
        ))

        sign = Int(sign_raw)
        sign ∈ (-1, 1) || throw(ArgumentError(
            "SignedPermutation signs must be ±1; got $sign for source index $src."
        ))

        dst = _signed_permutation_require_integer(dst_raw, "target index"; src)
        return src, sign, dst
    end

    dst = _signed_permutation_require_integer(value, "target index"; src)
    return src, 1, dst
end

function SignedPermutation(pairs::Pair...)
    isempty(pairs) && return SignedPermutation(Dict{Int,Tuple{Int,Int}}())

    parsed = map(_parse_signed_permutation_pair, pairs)
    types = Type[]
    for (src, _, dst) in parsed
        push!(types, typeof(src), typeof(dst))
    end
    T = promote_type(types...)

    images = Dict{T,Tuple{Int,T}}()
    for (src, sign, dst) in parsed
        images[convert(T, src)] = (sign, convert(T, dst))
    end
    return SignedPermutation(images)
end

"""
    SymmetrySpec

Symmetry configuration for the dense signed-permutation MVP.

# Fields
- `generators`: signed-permutation generators on registry indices
- `check_invariance`: whether to fail fast on objective/constraint non-invariance
"""
@kwdef struct SymmetrySpec{T<:Integer}
    generators::Vector{SignedPermutation{T}}
    check_invariance::Bool = true
end

function SymmetrySpec(generators::AbstractVector{<:SignedPermutation}; check_invariance::Bool=true)
    isempty(generators) && throw(ArgumentError("`SymmetrySpec` needs at least one generator."))

    T = Int
    for generator in generators
        for (src, (_, dst)) in generator.images
            T = promote_type(T, typeof(src), typeof(dst))
        end
    end

    converted = SignedPermutation{T}[
        _convert_signed_permutation(T, generator) for generator in generators
    ]
    return SymmetrySpec{T}(converted, check_invariance)
end

function SymmetrySpec(generators::SignedPermutation...; check_invariance::Bool=true)
    return SymmetrySpec(collect(generators); check_invariance)
end

"""
    SymmetryReport

Summary of the symmetry reduction applied to a relaxation.
"""
struct SymmetryReport
    group_order::Int
    invariant_moment_count::Int
    psd_block_sizes::Vector{Int}
    basis_half_size::Int
    basis_full_size::Int
end

function Base.show(io::IO, report::SymmetryReport)
    print(
        io,
        "SymmetryReport(group_order=$(report.group_order), " *
        "invariant_moment_count=$(report.invariant_moment_count), " *
        "psd_block_sizes=$(report.psd_block_sizes), " *
        "basis_half_size=$(report.basis_half_size), " *
        "basis_full_size=$(report.basis_full_size))"
    )
end

struct _OrbitReducer{A<:AlgebraType,T<:Integer}
    representative::Dict{NormalMonomial{A,T},NormalMonomial{A,T}}
    phase::Dict{NormalMonomial{A,T},Int}
end

struct NCWordSignedPermutationAction{A<:MonoidAlgebra,T<:Integer} <: SymbolicWedderburn.BySignedPermutations
    algebra::Type{A}
end

struct _SignedPermutationGroup{T<:Integer} <: GroupsCore.Group
    domain::Vector{T}
    elements::Vector{SignedPermutation{T}}
    generator_indices::Vector{Int}
    mult_table::Matrix{Int}
    inv_table::Vector{Int}
    element_orders::Vector{Int}
end

struct _SignedPermutationGroupElement{T<:Integer} <: GroupsCore.GroupElement
    group::_SignedPermutationGroup{T}
    index::Int
end

@inline _signed_image(g::SignedPermutation{T}, idx::T) where {T<:Integer} = get(g.images, idx, (1, idx))
@inline _signed_image_signature(g::SignedPermutation{T}, idx::T) where {T<:Integer} = begin
    sign, dst = _signed_image(g, idx)
    return sign * Int(dst)
end

function _convert_signed_permutation(::Type{T}, g::SignedPermutation) where {T<:Integer}
    images = Dict{T,Tuple{Int,T}}()
    for (src, (sign, dst)) in g.images
        images[convert(T, src)] = (sign, convert(T, dst))
    end
    return SignedPermutation(images)
end

function _convert_symmetry_spec(::Type{T}, spec::SymmetrySpec) where {T<:Integer}
    converted = SignedPermutation{T}[_convert_signed_permutation(T, g) for g in spec.generators]
    return SymmetrySpec{T}(converted, spec.check_invariance)
end

function _validate_signed_permutation(g::SignedPermutation{T}, domain::AbstractVector{T}) where {T<:Integer}
    domain_set = Set(domain)
    targets = T[]
    sizehint!(targets, length(domain))

    for idx in domain
        sign, dst = _signed_image(g, idx)
        sign ∈ (-1, 1) || throw(ArgumentError(
            "SignedPermutation signs must be ±1; got $sign for index $idx."
        ))
        dst in domain_set || throw(ArgumentError(
            "SignedPermutation sends index $idx outside the active symmetry domain: $dst."
        ))
        push!(targets, dst)
    end

    length(unique(targets)) == length(domain) || throw(ArgumentError(
        "SignedPermutation must be bijective on the active symmetry domain."
    ))

    return nothing
end

function _compose_signed_permutations(
    left::SignedPermutation{T},
    right::SignedPermutation{T},
    domain::AbstractVector{T},
) where {T<:Integer}
    images = Dict{T,Tuple{Int,T}}()
    for idx in domain
        sign_right, mid = _signed_image(right, idx)
        sign_left, dst = _signed_image(left, mid)
        sign = sign_left * sign_right
        if sign != 1 || dst != idx
            images[idx] = (sign, dst)
        end
    end
    return SignedPermutation(images)
end

@inline _identity_signed_permutation(::Type{T}) where {T<:Integer} =
    SignedPermutation(Dict{T,Tuple{Int,T}}())

function _symmetry_key(g::SignedPermutation{T}, domain::AbstractVector{T}) where {T<:Integer}
    return Tuple(_signed_image_signature(g, idx) for idx in domain)
end

function _enumerate_symmetry_group(spec::SymmetrySpec{T}, domain::Vector{T}) where {T<:Integer}
    isempty(spec.generators) && throw(ArgumentError("`SymmetrySpec` needs at least one generator."))
    isempty(domain) && throw(ArgumentError("Symmetry reduction needs a non-empty active variable domain."))

    for generator in spec.generators
        _validate_signed_permutation(generator, domain)
    end

    identity = _identity_signed_permutation(T)
    seen = Dict{Any,SignedPermutation{T}}(_symmetry_key(identity, domain) => identity)
    queue = SignedPermutation{T}[identity]

    while !isempty(queue)
        current = popfirst!(queue)
        for generator in spec.generators
            candidate = _compose_signed_permutations(generator, current, domain)
            key = _symmetry_key(candidate, domain)
            if !haskey(seen, key)
                seen[key] = candidate
                push!(queue, candidate)
            end
        end
    end

    ordered_keys = sort!(collect(keys(seen)))
    return [seen[key] for key in ordered_keys]
end

function _inverse_signed_permutation(g::SignedPermutation{T}, domain::AbstractVector{T}) where {T<:Integer}
    images = Dict{T,Tuple{Int,T}}()
    for idx in domain
        sign, dst = _signed_image(g, idx)
        if sign != 1 || idx != dst
            images[dst] = (sign, idx)
        end
    end
    return SignedPermutation(images)
end

@inline _sw_group_value(g::_SignedPermutationGroupElement{T}) where {T<:Integer} = g.group.elements[g.index]

Base.eltype(::Type{_SignedPermutationGroup{T}}) where {T<:Integer} = _SignedPermutationGroupElement{T}
Base.IteratorSize(::Type{<:_SignedPermutationGroup}) = Base.HasLength()
Base.length(G::_SignedPermutationGroup) = length(G.elements)

function Base.iterate(G::_SignedPermutationGroup{T}, state::Int=1) where {T<:Integer}
    state > length(G.elements) && return nothing
    return _SignedPermutationGroupElement{T}(G, state), state + 1
end

Base.one(G::_SignedPermutationGroup{T}) where {T<:Integer} = _SignedPermutationGroupElement{T}(G, 1)
Base.parent(g::_SignedPermutationGroupElement) = g.group
Base.copy(g::_SignedPermutationGroupElement) = g
Base.isone(g::_SignedPermutationGroupElement) = g.index == 1

function Base.:(==)(a::_SignedPermutationGroupElement{T}, b::_SignedPermutationGroupElement{T}) where {T<:Integer}
    return a.group.domain == b.group.domain && _sw_group_value(a) == _sw_group_value(b)
end

function Base.hash(g::_SignedPermutationGroupElement, h::UInt)
    return hash((g.group.domain, _sw_group_value(g)), h)
end

function GroupsCore.order(::Type{I}, G::_SignedPermutationGroup) where {I<:Integer}
    return convert(I, length(G.elements))
end

function GroupsCore.order(::Type{I}, g::_SignedPermutationGroupElement) where {I<:Integer}
    return convert(I, g.group.element_orders[g.index])
end

function GroupsCore.gens(G::_SignedPermutationGroup{T}) where {T<:Integer}
    return [_SignedPermutationGroupElement{T}(G, idx) for idx in G.generator_indices]
end

function Base.inv(g::_SignedPermutationGroupElement{T}) where {T<:Integer}
    return _SignedPermutationGroupElement{T}(g.group, g.group.inv_table[g.index])
end

function Base.:(*)(a::_SignedPermutationGroupElement{T}, b::_SignedPermutationGroupElement{T}) where {T<:Integer}
    a.group.domain == b.group.domain && a.group.elements == b.group.elements || throw(ArgumentError(
        "Cannot multiply elements from different signed-permutation symmetry groups."
    ))
    return _SignedPermutationGroupElement{T}(a.group, a.group.mult_table[a.index, b.index])
end

function _sw_signed_permutation_group(
    spec::SymmetrySpec{T},
    group::Vector{SignedPermutation{T}},
    domain::Vector{T},
) where {T<:Integer}
    identity = _identity_signed_permutation(T)
    identity_key = _symmetry_key(identity, domain)
    elements = SignedPermutation{T}[identity]
    for g in group
        _symmetry_key(g, domain) == identity_key && continue
        push!(elements, g)
    end

    lookup = Dict{Any,Int}(_symmetry_key(g, domain) => i for (i, g) in enumerate(elements))
    generator_indices = Int[]
    for generator in spec.generators
        idx = get(lookup, _symmetry_key(generator, domain), 0)
        idx == 0 && throw(ArgumentError(
            "Internal symmetry error: failed to locate a generator in the enumerated symmetry group."
        ))
        push!(generator_indices, idx)
    end

    n = length(elements)
    mult_table = Matrix{Int}(undef, n, n)
    for j in 1:n, i in 1:n
        product = _compose_signed_permutations(elements[j], elements[i], domain)
        idx = get(lookup, _symmetry_key(product, domain), 0)
        idx == 0 && throw(ArgumentError(
            "Internal symmetry error: enumerated symmetry group is not closed under multiplication."
        ))
        mult_table[i, j] = idx
    end

    inv_table = Vector{Int}(undef, n)
    for i in 1:n
        inverse = _inverse_signed_permutation(elements[i], domain)
        idx = get(lookup, _symmetry_key(inverse, domain), 0)
        idx == 0 && throw(ArgumentError(
            "Internal symmetry error: enumerated symmetry group is missing an inverse element."
        ))
        inv_table[i] = idx
    end

    element_orders = ones(Int, n)
    for i in 2:n
        order_i = 1
        power = i
        while power != 1
            power = mult_table[power, i]
            order_i += 1
            order_i <= n || throw(ArgumentError(
                "Internal symmetry error: failed to compute an element order for the signed-permutation group."
            ))
        end
        element_orders[i] = order_i
    end

    return _SignedPermutationGroup(domain, elements, generator_indices, mult_table, inv_table, element_orders)
end

function SymbolicWedderburn.action(
    ::NCWordSignedPermutationAction{A,T},
    g::_SignedPermutationGroupElement{T},
    mono::NormalMonomial{A,T},
) where {A<:MonoidAlgebra,T<:Integer}
    sign, image = _act_monomial(_sw_group_value(g), mono)
    return image, sign
end

function _sw_decompose_half_basis(
    basis::Vector{M},
    group::_SignedPermutationGroup{T},
) where {A<:MonoidAlgebra,T<:Integer,M<:NormalMonomial{A,T}}
    action = NCWordSignedPermutationAction{A,T}(A)
    blocks = SymbolicWedderburn.symmetry_adapted_basis(Float64, group, action, basis)

    for block in blocks
        SymbolicWedderburn.multiplicity(block) == 1 || throw(ArgumentError(
            "Symmetry reduction MVP currently supports multiplicity-free reductions only. " *
            "SymbolicWedderburn found a direct summand with multiplicity $(SymbolicWedderburn.multiplicity(block))."
        ))
        SymbolicWedderburn.issimple(block) || throw(ArgumentError(
            "Symmetry reduction MVP currently supports only scalar reduced PSD blocks. " *
            "SymbolicWedderburn found a non-simple direct summand of size $(size(block, 1))."
        ))
        size(block, 1) == 1 || throw(ArgumentError(
            "Symmetry reduction MVP expected a 1×1 simple direct summand, got size $(size(block, 1))."
        ))
    end

    return blocks
end

function _symmetry_domain(
    pop::PolyOpt{A,T,P},
    corr_sparsity::CorrelativeSparsity{A,T,P,M,Nothing},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
) where {A<:MonoidAlgebra,T<:Integer,C<:Number,P<:Polynomial{A,T,C},M<:NormalMonomial{A,T}}
    used = Set{T}()

    for poly in (pop.objective,)
        union!(used, variable_indices(poly))
    end
    for poly in pop.eq_constraints
        union!(used, variable_indices(poly))
    end
    for poly in pop.ineq_constraints
        union!(used, variable_indices(poly))
    end
    for poly in pop.moment_eq_constraints
        union!(used, variable_indices(poly))
    end

    for basis in corr_sparsity.clq_mom_mtx_bases
        for mono in basis
            union!(used, variable_indices(mono))
        end
    end
    for term_sparsities in cliques_term_sparsities, term_sparsity in term_sparsities, block_basis in term_sparsity.block_bases, mono in block_basis
        union!(used, variable_indices(mono))
    end

    return sort!(collect(used))
end

function _act_monomial(
    g::SignedPermutation{T},
    mono::NormalMonomial{A,T},
) where {A<:MonoidAlgebra,T<:Integer}
    sign = 1
    word = similar(mono.word)
    for (i, idx) in pairs(mono.word)
        letter_sign, dst = _signed_image(g, idx)
        sign *= letter_sign
        word[i] = dst
    end
    image = NormalMonomial{A,T}(simplify(A, word))
    return sign, image
end

function _act_polynomial(
    g::SignedPermutation{T},
    poly::Polynomial{A,T,C},
) where {A<:MonoidAlgebra,T<:Integer,C<:Number}
    acted_terms = Tuple{C,NormalMonomial{A,T}}[]
    sizehint!(acted_terms, length(poly.terms))
    for (coef, mono) in poly.terms
        sign, image = _act_monomial(g, mono)
        push!(acted_terms, (coef * sign, image))
    end
    return Polynomial(acted_terms)
end

function _check_polynomial_invariance(
    label::AbstractString,
    poly::Polynomial{A,T,C},
    group::AbstractVector{<:SignedPermutation{T}},
) where {A<:MonoidAlgebra,T<:Integer,C<:Number}
    for g in group
        _act_polynomial(g, poly) == poly || throw(ArgumentError(
            "Symmetry reduction requires $label to be invariant under the supplied symmetry."
        ))
    end
    return nothing
end

function _check_symmetry_invariance(
    pop::PolyOpt{A,T,P},
    group::AbstractVector{<:SignedPermutation{T}},
) where {A<:MonoidAlgebra,T<:Integer,C<:Number,P<:Polynomial{A,T,C}}
    _check_polynomial_invariance("the objective", pop.objective, group)

    for (i, poly) in pairs(pop.eq_constraints)
        _check_polynomial_invariance("equality constraint $i", poly, group)
    end
    for (i, poly) in pairs(pop.ineq_constraints)
        _check_polynomial_invariance("inequality constraint $i", poly, group)
    end
    for (i, poly) in pairs(pop.moment_eq_constraints)
        _check_polynomial_invariance("moment equality constraint $i", poly, group)
    end

    return nothing
end

function _check_basis_closure(
    label::AbstractString,
    basis::Vector{M},
    group::AbstractVector{<:SignedPermutation{T}},
) where {A<:MonoidAlgebra,T<:Integer,M<:NormalMonomial{A,T}}
    lookup = Set(basis)
    for g in group, mono in basis
        _, image = _act_monomial(g, mono)
        image in lookup || throw(ArgumentError(
            "Symmetry reduction requires closure of the action on $label. " *
            "Basis element $mono maps outside the basis to $image."
        ))
    end
    return nothing
end

@inline function _symmetric_monomial(mono::NormalMonomial{A,T}) where {A<:MonoidAlgebra,T<:Integer}
    return NormalMonomial{A,T}(symmetric_canon(mono))
end

function _build_orbit_reducer(
    basis::AbstractVector{<:NormalMonomial{A,T}},
    group::AbstractVector{<:SignedPermutation{T}},
) where {A<:MonoidAlgebra,T<:Integer}
    sym_basis = sorted_unique!(_symmetric_monomial.(collect(basis)))
    sym_basis_set = Set(sym_basis)

    representative = Dict{NormalMonomial{A,T},NormalMonomial{A,T}}()
    phase = Dict{NormalMonomial{A,T},Int}()
    unvisited = Set(sym_basis)

    while !isempty(unvisited)
        start = first(unvisited)
        orbit_phase = Dict(start => 1)
        orbit = Set([start])
        queue = NormalMonomial{A,T}[start]
        zero_orbit = false

        while !isempty(queue)
            current = popfirst!(queue)
            for g in group
                sign, image = _act_monomial(g, current)
                image_sym = _symmetric_monomial(image)
                image_sym in sym_basis_set || throw(ArgumentError(
                    "Symmetry action does not preserve the constructed moment basis. " *
                    "Missing orbit image $image_sym."
                ))

                candidate_phase = orbit_phase[current] * sign
                if haskey(orbit_phase, image_sym)
                    orbit_phase[image_sym] == candidate_phase || (zero_orbit = true)
                else
                    orbit_phase[image_sym] = candidate_phase
                    push!(orbit, image_sym)
                    push!(queue, image_sym)
                end
            end
        end

        orbit_rep = minimum(collect(orbit))
        rep_phase = orbit_phase[orbit_rep]
        for mono in orbit
            representative[mono] = orbit_rep
            phase[mono] = zero_orbit ? 0 : orbit_phase[mono] * rep_phase
            delete!(unvisited, mono)
        end
    end

    return _OrbitReducer{A,T}(representative, phase)
end

function _chop_polynomial(
    poly::Polynomial{A,T,C};
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:MonoidAlgebra,T<:Integer,C<:Number}
    chopped_terms = Tuple{C,NormalMonomial{A,T}}[]
    sizehint!(chopped_terms, length(poly.terms))
    for (coef, mono) in poly.terms
        abs(coef) > atol && push!(chopped_terms, (coef, mono))
    end
    return Polynomial(chopped_terms)
end

function _orbit_reduce_polynomial(
    poly::Polynomial{A,T,C},
    reducer::_OrbitReducer{A,T};
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:MonoidAlgebra,T<:Integer,C<:Number}
    reduced_terms = Tuple{C,NormalMonomial{A,T}}[]
    sizehint!(reduced_terms, length(poly.terms))

    for (coef, mono) in poly.terms
        sym_mono = _symmetric_monomial(mono)
        haskey(reducer.representative, sym_mono) || throw(ArgumentError(
            "Symmetry reducer is missing the monomial $sym_mono."
        ))
        orbit_phase = reducer.phase[sym_mono]
        orbit_phase == 0 && continue
        push!(reduced_terms, (coef * orbit_phase, reducer.representative[sym_mono]))
    end

    return _chop_polynomial(Polynomial(reduced_terms); atol)
end

function _orbit_reduce_matrix(
    mat::Matrix{P},
    reducer::_OrbitReducer{A,T};
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:MonoidAlgebra,T<:Integer,C<:Number,P<:Polynomial{A,T,C}}
    reduced = Matrix{P}(undef, size(mat, 1), size(mat, 2))
    for j in axes(mat, 2), i in axes(mat, 1)
        reduced[i, j] = _orbit_reduce_polynomial(_chop_polynomial(mat[i, j]; atol), reducer; atol)
    end
    return reduced
end

function _poly_is_small(
    poly::Polynomial{A,T,C};
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:MonoidAlgebra,T<:Integer,C<:Number}
    return iszero(_chop_polynomial(poly; atol))
end

function _assert_poly_small(
    poly::Polynomial{A,T,C},
    context::AbstractString;
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:MonoidAlgebra,T<:Integer,C<:Number}
    _poly_is_small(poly; atol) || throw(ArgumentError(
        "$context should vanish after symmetry reduction, but got $poly."
    ))
    return nothing
end

function _assert_poly_equal(
    lhs::Polynomial{A,T,C},
    rhs::Polynomial{A,T,C},
    context::AbstractString;
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:MonoidAlgebra,T<:Integer,C<:Number}
    diff = _chop_polynomial(lhs - rhs; atol)
    iszero(diff) || throw(ArgumentError(
        "$context should match after symmetry reduction, but differs by $diff."
    ))
    return nothing
end

function _transform_polynomial_block(
    mat::Matrix{P},
    U_left::AbstractMatrix{<:Real},
    U_right::AbstractMatrix{<:Real};
    scale::Real=1,
    atol::Float64=_SYMMETRY_ATOL,
) where {P<:Polynomial}
    transformed = Matrix{P}(undef, size(U_left, 1), size(U_right, 1))

    for j in axes(transformed, 2), i in axes(transformed, 1)
        acc = zero(P)
        for b in axes(mat, 2), a in axes(mat, 1)
            coeff = U_left[i, a] * U_right[j, b]
            abs(coeff) <= atol && continue
            acc += coeff * mat[a, b]
        end
        transformed[i, j] = _chop_polynomial(scale == 1 ? acc : scale * acc; atol)
    end

    return transformed
end

function _reduce_to_scalar_blocks_sw(
    mat::Matrix{P},
    blocks,
    reducer::_OrbitReducer{A,T};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="constraint matrix",
) where {A<:MonoidAlgebra,T<:Integer,C<:Number,P<:Polynomial{A,T,C}}
    row_bases = [Matrix(block) for block in blocks]

    for i in 1:length(row_bases), j in (i+1):length(row_bases)
        off_block = _orbit_reduce_matrix(
            _transform_polynomial_block(mat, row_bases[i], row_bases[j]; atol),
            reducer;
            atol,
        )
        for col in axes(off_block, 2), row in axes(off_block, 1)
            _assert_poly_small(
                off_block[row, col],
                "$context off-block entry ($i, $j)[$row, $col]";
                atol,
            )
        end

        off_block_t = _orbit_reduce_matrix(
            _transform_polynomial_block(mat, row_bases[j], row_bases[i]; atol),
            reducer;
            atol,
        )
        for col in axes(off_block_t, 2), row in axes(off_block_t, 1)
            _assert_poly_small(
                off_block_t[row, col],
                "$context off-block entry ($j, $i)[$row, $col]";
                atol,
            )
        end
    end

    scalar_blocks = Matrix{P}[]
    for (block_idx, block) in enumerate(blocks)
        diagonal = _orbit_reduce_matrix(
            _transform_polynomial_block(
                mat,
                row_bases[block_idx],
                row_bases[block_idx];
                scale=SymbolicWedderburn.degree(block),
                atol,
            ),
            reducer;
            atol,
        )
        size(diagonal) == (1, 1) || throw(ArgumentError(
            "$context should reduce to a scalar PSD block, got a block of size $(size(diagonal))."
        ))
        push!(scalar_blocks, reshape([_chop_polynomial(only(diagonal); atol)], 1, 1))
    end

    return scalar_blocks
end

function _reduce_constraint_matrix_symmetric(
    mat::Matrix{P},
    basis::Vector{M},
    sw_group::_SignedPermutationGroup{T},
    reducer::_OrbitReducer{A,T};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="constraint matrix",
) where {A<:MonoidAlgebra,T<:Integer,C<:Number,P<:Polynomial{A,T,C},M<:NormalMonomial{A,T}}
    if length(basis) == 1
        reduced = _orbit_reduce_matrix(mat, reducer; atol)
        return [reshape([_chop_polynomial(only(reduced); atol)], 1, 1)]
    end

    blocks = _sw_decompose_half_basis(basis, sw_group)
    return _reduce_to_scalar_blocks_sw(mat, blocks, reducer; atol, context)
end

function _collect_reduced_basis(
    objective::Polynomial{A,T,C},
    constraints::Vector{Tuple{Symbol, Matrix{Polynomial{A,T,C}}}},
) where {A<:MonoidAlgebra,T<:Integer,C<:Number}
    basis = NormalMonomial{A,T}[]
    append!(basis, monomials(objective))
    for (_, mat) in constraints, poly in mat
        append!(basis, monomials(poly))
    end
    return sorted_unique!(basis)
end

function _add_moment_eq_constraints_symmetric!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    pop::PolyOpt{A,T,PP},
    moment_eq_row_bases::Vector{M},
    moment_eq_row_basis_degrees::Vector{Int},
    reducer::_OrbitReducer{A,T},
    ::Type{MP};
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:MonoidAlgebra,T<:Integer,CMP<:Number,CPP<:Number,M<:NormalMonomial{A,T},MP<:Polynomial{A,T,CMP},PP<:Polynomial{A,T,CPP}}
    isempty(pop.moment_eq_constraints) && return nothing

    one_mono = one(NormalMonomial{A,T})
    for g in pop.moment_eq_constraints
        row_bases = _truncate_moment_eq_row_bases(moment_eq_row_bases, moment_eq_row_basis_degrees, g)
        isempty(row_bases) && continue

        for row_mono in row_bases
            poly = sum(
                _conj_coef(A, c_row) * coef * Polynomial(_simplified_to_terms(A, simplify(A, _neat_dot3(row_word, mono, one_mono)), T))
                for (c_row, row_word) in row_mono
                for (coef, mono) in g.terms;
                init=zero(MP)
            )
            poly = _orbit_reduce_polynomial(_chop_polynomial(poly; atol), reducer; atol)
            iszero(poly) && continue
            push!(constraints, (:Zero, reshape([poly], 1, 1)))
        end
    end

    return nothing
end


"""
    moment_relax_symmetric(pop, corr_sparsity, cliques_term_sparsities, symmetry)

Construct a symmetry-reduced symbolic moment problem for the dense monoid MVP.
"""
function moment_relax_symmetric(
    pop::PolyOpt{A,T,P},
    corr_sparsity::CorrelativeSparsity{A,T,P,M,Nothing},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
    symmetry::SymmetrySpec;
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:MonoidAlgebra,T<:Integer,C<:Number,P<:Polynomial{A,T,C},M<:NormalMonomial{A,T}}
    spec = _convert_symmetry_spec(T, symmetry)
    domain = _symmetry_domain(pop, corr_sparsity, cliques_term_sparsities)
    group = _enumerate_symmetry_group(spec, domain)
    sw_group = _sw_signed_permutation_group(spec, group, domain)

    spec.check_invariance && _check_symmetry_invariance(pop, group)

    for (clique_idx, basis) in enumerate(corr_sparsity.clq_mom_mtx_bases)
        _check_basis_closure("moment basis of clique $clique_idx", basis, group)
    end
    for (clique_idx, term_sparsities) in enumerate(cliques_term_sparsities)
        for (poly_idx, term_sparsity) in enumerate(term_sparsities)
            _check_basis_closure(
                "constraint block basis $(poly_idx) of clique $clique_idx",
                only(term_sparsity.block_bases),
                group,
            )
        end
    end

    moment_matrix_basis = _moment_matrix_basis(cliques_term_sparsities)
    total_basis, moment_eq_row_bases, moment_eq_row_basis_degrees =
        _polynomial_total_basis(pop, corr_sparsity, cliques_term_sparsities)
    _validate_polynomial_relaxation_support(pop, total_basis; source="Constructed relaxation basis")

    reducer = _build_orbit_reducer(total_basis, group)

    psd_cone = _is_complex_problem(A) ? :HPSD : :PSD
    MP_C = _moment_problem_coeff_type(A, C)
    MP_P = Polynomial{A,T,MP_C}
    objective_mp = _orbit_reduce_polynomial(convert(MP_P, pop.objective), reducer; atol)

    constraints = Tuple{Symbol, Matrix{MP_P}}[]
    reduced_moment_terms = NormalMonomial{A,T}[]
    reduced_moment_block_sizes = Int[]

    for (clique_idx, (term_sparsities, cons_idx)) in enumerate(zip(cliques_term_sparsities, corr_sparsity.clq_cons))
        polys = [one(pop.objective); corr_sparsity.cons[cons_idx]...]

        for (poly_idx, (term_sparsity, poly)) in enumerate(zip(term_sparsities, polys))
            for (block_idx, basis) in enumerate(term_sparsity.block_bases)
                cone = poly in pop.eq_constraints ? :Zero : psd_cone
                _, mat = _build_constraint_matrix(poly, basis, cone)
                scalar_blocks = _reduce_constraint_matrix_symmetric(
                    mat,
                    basis,
                    sw_group,
                    reducer;
                    atol,
                    context="clique $clique_idx block $block_idx for polynomial $poly_idx",
                )

                for scalar_block in scalar_blocks
                    _append_constraint!(constraints, cone, scalar_block, MP_P)
                    poly_idx == 1 && append!(reduced_moment_terms, monomials(only(scalar_block)))
                end

                poly_idx == 1 && append!(reduced_moment_block_sizes, fill(1, length(scalar_blocks)))
            end
        end
    end

    for global_con in corr_sparsity.global_cons
        poly = corr_sparsity.cons[global_con]
        cone = poly in pop.eq_constraints ? :Zero : psd_cone
        reduced_poly = _orbit_reduce_polynomial(convert(MP_P, poly), reducer; atol)
        _append_constraint!(constraints, cone, reshape([reduced_poly], 1, 1), MP_P)
    end

    _add_moment_eq_constraints_symmetric!(
        constraints,
        pop,
        moment_eq_row_bases,
        moment_eq_row_basis_degrees,
        reducer,
        MP_P;
        atol,
    )

    reduced_total_basis = _collect_reduced_basis(objective_mp, constraints)
    reduced_moment_basis = sorted_unique!(_symmetric_monomial.(reduced_moment_terms))
    n_unique_moment_matrix_elements = length(reduced_moment_basis)

    mp = MomentProblem{A,T,M,MP_P}(
        objective_mp,
        constraints,
        reduced_total_basis,
        n_unique_moment_matrix_elements,
    )

    report = SymmetryReport(
        length(group),
        count(mono -> !isone(mono), reduced_moment_basis),
        reduced_moment_block_sizes,
        length(only(corr_sparsity.clq_mom_mtx_bases)),
        length(moment_matrix_basis),
    )

    return mp, report
end
