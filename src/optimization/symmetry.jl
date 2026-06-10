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

Use this for monoid-algebra registry-index symmetries. Pauli problems should
use [`CliffordSymmetry`](@ref), and fermionic problems should use
[`FermionicModePermutation`](@ref) plus optional sector/spin specs.
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
    FermionicModePermutation

Finite permutation acting on physical fermionic modes.

Unlike `SignedPermutation`, this acts on unsigned physical mode ids. Creation vs
annihilation is preserved by the operator action itself, and any reordering sign
comes from re-normal-ordering the fermionic word.
"""
struct FermionicModePermutation{T<:Integer}
    images::Dict{T,T}

    function FermionicModePermutation{T}(images::Dict{T,T}) where {T<:Integer}
        normalized = Dict{T,T}()
        for (src, dst) in images
            src == dst || (normalized[src] = dst)
        end
        return new{T}(normalized)
    end
end

function FermionicModePermutation(images::AbstractDict{T,<:Integer}) where {T<:Integer}
    normalized = Dict{T,T}()
    for (src, dst) in images
        normalized[src] = convert(T, dst)
    end
    return FermionicModePermutation{T}(normalized)
end

function FermionicModePermutation(pairs::Pair...)
    isempty(pairs) && return FermionicModePermutation(Dict{Int,Int}())
    parsed = Tuple{Integer,Integer}[]
    for pair in pairs
        src = _signed_permutation_require_integer(pair.first, "source mode")
        dst = _signed_permutation_require_integer(pair.second, "target mode"; src)
        push!(parsed, (src, dst))
    end

    types = Type[]
    for (src, dst) in parsed
        push!(types, typeof(src), typeof(dst))
    end
    T = promote_type(types...)

    images = Dict{T,T}()
    for (src, dst) in parsed
        images[convert(T, src)] = convert(T, dst)
    end
    return FermionicModePermutation(images)
end

const _PAULI_X_TYPE = 0
const _PAULI_Y_TYPE = 1
const _PAULI_Z_TYPE = 2

"""
    CliffordSymmetry

Finite Clifford conjugation action on Pauli words.

The action is stored as signed images of Pauli-basis words:
`source => (sign, target)`, with `sign ∈ {-1, +1}`. Omitted basis
letters are fixed. Named constructors provide the standard gate actions:
`CliffordSymmetry(:H, q)`, `CliffordSymmetry(:S, q)`,
`CliffordSymmetry(:CNOT, control, target)`, and
`CliffordSymmetry(:SWAP, q1, q2)`.
"""
struct CliffordSymmetry{T<:Integer}
    nqubits::Int
    images::Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}

    function CliffordSymmetry{T}(
        nqubits::Integer,
        images::Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}},
    ) where {T<:Integer}
        g = new{T}(_clifford_validate_nqubits(nqubits), images)
        _validate_clifford_symmetry(g)
        return g
    end
end

function _clifford_validate_nqubits(nqubits::Integer)
    n = Int(nqubits)
    n > 0 || throw(ArgumentError("CliffordSymmetry needs a positive number of qubits; got $n."))
    return n
end

function _clifford_validate_site(site::Integer)
    q = Int(site)
    q > 0 || throw(ArgumentError("CliffordSymmetry qubit sites must be positive; got $q."))
    return q
end

function _clifford_max_site(mono::NormalMonomial{PauliAlgebra})
    max_site = 0
    for idx in mono.word
        idx > 0 || throw(ArgumentError(
            "CliffordSymmetry image dictionary contains non-positive Pauli index $idx."
        ))
        max_site = max(max_site, _pauli_site(idx))
    end
    return max_site
end

function _clifford_copy_monomial(mono::NormalMonomial{PauliAlgebra,T}) where {T<:Integer}
    return NormalMonomial{PauliAlgebra,T}(copy(mono.word))
end

function _clifford_validate_pauli_indices(
    mono::NormalMonomial{PauliAlgebra},
    nqubits::Integer,
    label::AbstractString,
)
    for idx in mono.word
        idx > 0 || throw(ArgumentError(
            "CliffordSymmetry $label contains non-positive Pauli index $idx."
        ))
        site = _pauli_site(idx)
        site <= nqubits || throw(ArgumentError(
            "CliffordSymmetry $label references site $site outside nqubits=$nqubits."
        ))
    end
    return nothing
end

function _clifford_validate_domain_indices(domain::AbstractVector{<:Integer}, nqubits::Integer)
    for idx in domain
        idx > 0 || throw(ArgumentError(
            "CliffordSymmetry active domain contains non-positive Pauli index $idx."
        ))
        site = _pauli_site(idx)
        site <= nqubits || throw(ArgumentError(
            "CliffordSymmetry active domain references site $site outside nqubits=$nqubits."
        ))
    end
    return nothing
end

@inline _clifford_sign_phase(sign::Integer) = sign == 1 ? UInt8(0) : UInt8(2)

function _pauli_anticommutes(
    a::NormalMonomial{PauliAlgebra,T},
    b::NormalMonomial{PauliAlgebra,T},
) where {T<:Integer}
    parity = false
    i = firstindex(a.word)
    j = firstindex(b.word)

    while i <= lastindex(a.word) && j <= lastindex(b.word)
        ai = a.word[i]
        bj = b.word[j]
        site_a = _pauli_site(ai)
        site_b = _pauli_site(bj)

        if site_a == site_b
            parity = xor(parity, _pauli_type(ai) != _pauli_type(bj))
            i += 1
            j += 1
        elseif site_a < site_b
            i += 1
        else
            j += 1
        end
    end

    return parity
end

function _clifford_letter_action(
    g::CliffordSymmetry{T},
    site::Integer,
    pauli_type::Integer,
) where {T<:Integer}
    idx = convert(T, _pauli_index(site, pauli_type))
    return _act_monomial(g, NormalMonomial{PauliAlgebra,T}(T[idx]))
end

function _assert_clifford_preserves_local_product(
    g::CliffordSymmetry{T},
    site::Integer,
    left_type::Int,
    right_type::Int,
) where {T<:Integer}
    product_phase, result_type = _pauli_product(left_type, right_type)
    result_type == -1 && return nothing

    left_sign, left_image = _clifford_letter_action(g, site, left_type)
    right_sign, right_image = _clifford_letter_action(g, site, right_type)
    result_sign, result_image = _clifford_letter_action(g, site, result_type)

    product_word, image_phase = simplify(PauliAlgebra, vcat(left_image.word, right_image.word))
    product_image = NormalMonomial{PauliAlgebra,T}(product_word)
    product_image == result_image || throw(ArgumentError(
        "CliffordSymmetry images must preserve Pauli multiplication on site $site."
    ))

    lhs_phase = (image_phase + _clifford_sign_phase(left_sign * right_sign)) % UInt8(4)
    rhs_phase = (product_phase + _clifford_sign_phase(result_sign)) % UInt8(4)
    lhs_phase == rhs_phase || throw(ArgumentError(
        "CliffordSymmetry images must preserve Pauli multiplication phases on site $site."
    ))

    return nothing
end

function _clifford_infer_nqubits(
    images::AbstractDict{NormalMonomial{PauliAlgebra,T},<:Tuple{<:Integer,NormalMonomial{PauliAlgebra,T}}},
    nqubits,
) where {T<:Integer}
    inferred = 0
    for (src, (_, dst)) in images
        inferred = max(inferred, _clifford_max_site(src), _clifford_max_site(dst))
    end

    if isnothing(nqubits)
        inferred > 0 || throw(ArgumentError(
            "CliffordSymmetry cannot infer `nqubits` from an empty image dictionary."
        ))
        return inferred
    end

    n = _clifford_validate_nqubits(nqubits)
    inferred <= n || throw(ArgumentError(
        "CliffordSymmetry image dictionary references site $inferred outside nqubits=$n."
    ))
    return n
end

function CliffordSymmetry(
    images::AbstractDict{NormalMonomial{PauliAlgebra,T},<:Tuple{<:Integer,NormalMonomial{PauliAlgebra,T}}};
    nqubits=nothing,
) where {T<:Integer}
    n = _clifford_infer_nqubits(images, nqubits)
    normalized = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()

    for (src_raw, (sign_raw, dst_raw)) in images
        sign = Int(sign_raw)
        sign ∈ (-1, 1) || throw(ArgumentError(
            "CliffordSymmetry signs must be ±1; got $sign for source word $src_raw."
        ))

        src = _clifford_copy_monomial(src_raw)
        dst = _clifford_copy_monomial(dst_raw)

        _clifford_validate_pauli_indices(src, n, "source word")
        _clifford_validate_pauli_indices(dst, n, "target word")

        if isone(src)
            sign == 1 && isone(dst) || throw(ArgumentError("CliffordSymmetry must fix the identity word."))
            continue
        end
        length(src.word) == 1 || throw(ArgumentError(
            "CliffordSymmetry image dictionaries currently accept only single-Pauli-letter sources; got $src."
        ))
        !isone(dst) || throw(ArgumentError(
            "CliffordSymmetry must map non-identity Pauli letters to non-identity Pauli words."
        ))

        if sign != 1 || dst != src
            normalized[src] = (sign, dst)
        end
    end

    return CliffordSymmetry{T}(n, normalized)
end

function _clifford_gate_nqubits(nqubits, sites::Integer...)
    max_site = maximum(_clifford_validate_site(site) for site in sites)
    isnothing(nqubits) && return max_site
    n = _clifford_validate_nqubits(nqubits)
    max_site <= n || throw(ArgumentError(
        "CliffordSymmetry gate references site $max_site outside nqubits=$n."
    ))
    return n
end

function _pauli_letter(::Type{T}, site::Integer, pauli_type::Integer) where {T<:Integer}
    return NormalMonomial{PauliAlgebra,T}(T[convert(T, _pauli_index(site, pauli_type))])
end

function _pauli_word(::Type{T}, letters::Pair{<:Integer,<:Integer}...) where {T<:Integer}
    raw = T[convert(T, _pauli_index(site, pauli_type)) for (site, pauli_type) in letters]
    word, phase = simplify(PauliAlgebra, raw)
    sign = _phase_to_clifford_sign(phase)
    return sign, NormalMonomial{PauliAlgebra,T}(word)
end

function _pauli_unsigned_word(::Type{T}, letters::Pair{<:Integer,<:Integer}...) where {T<:Integer}
    sign, mono = _pauli_word(T, letters...)
    sign == 1 || throw(ArgumentError(
        "Internal Clifford gate image construction produced a negative Pauli phase."
    ))
    return mono
end

function CliffordSymmetry(
    gate::Symbol,
    sites::Integer...;
    nqubits=nothing,
    integer_type::Type{T}=Int,
) where {T<:Integer}
    gate_key = Symbol(uppercase(String(gate)))

    if gate_key === :H
        length(sites) == 1 || throw(ArgumentError("CliffordSymmetry(:H, q) needs exactly one qubit site."))
        q = _clifford_validate_site(only(sites))
        n = _clifford_gate_nqubits(nqubits, q)
        x = _pauli_letter(T, q, _PAULI_X_TYPE)
        y = _pauli_letter(T, q, _PAULI_Y_TYPE)
        z = _pauli_letter(T, q, _PAULI_Z_TYPE)
        return CliffordSymmetry(Dict(x => (1, z), y => (-1, y), z => (1, x)); nqubits=n)
    elseif gate_key === :S
        length(sites) == 1 || throw(ArgumentError("CliffordSymmetry(:S, q) needs exactly one qubit site."))
        q = _clifford_validate_site(only(sites))
        n = _clifford_gate_nqubits(nqubits, q)
        x = _pauli_letter(T, q, _PAULI_X_TYPE)
        y = _pauli_letter(T, q, _PAULI_Y_TYPE)
        return CliffordSymmetry(Dict(x => (1, y), y => (-1, x)); nqubits=n)
    elseif gate_key === :CNOT
        length(sites) == 2 || throw(ArgumentError(
            "CliffordSymmetry(:CNOT, control, target) needs exactly two qubit sites."
        ))
        control = _clifford_validate_site(sites[1])
        target = _clifford_validate_site(sites[2])
        control != target || throw(ArgumentError("CliffordSymmetry(:CNOT, control, target) needs distinct sites."))
        n = _clifford_gate_nqubits(nqubits, control, target)

        x_control = _pauli_letter(T, control, _PAULI_X_TYPE)
        y_control = _pauli_letter(T, control, _PAULI_Y_TYPE)
        z_target = _pauli_letter(T, target, _PAULI_Z_TYPE)
        y_target = _pauli_letter(T, target, _PAULI_Y_TYPE)

        images = Dict(
            x_control => (1, _pauli_unsigned_word(T, control => _PAULI_X_TYPE, target => _PAULI_X_TYPE)),
            y_control => (1, _pauli_unsigned_word(T, control => _PAULI_Y_TYPE, target => _PAULI_X_TYPE)),
            z_target => (1, _pauli_unsigned_word(T, control => _PAULI_Z_TYPE, target => _PAULI_Z_TYPE)),
            y_target => (1, _pauli_unsigned_word(T, control => _PAULI_Z_TYPE, target => _PAULI_Y_TYPE)),
        )
        return CliffordSymmetry(images; nqubits=n)
    elseif gate_key === :SWAP
        length(sites) == 2 || throw(ArgumentError("CliffordSymmetry(:SWAP, q1, q2) needs exactly two qubit sites."))
        q1 = _clifford_validate_site(sites[1])
        q2 = _clifford_validate_site(sites[2])
        q1 != q2 || throw(ArgumentError("CliffordSymmetry(:SWAP, q1, q2) needs distinct sites."))
        n = _clifford_gate_nqubits(nqubits, q1, q2)

        images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
        for pauli_type in (_PAULI_X_TYPE, _PAULI_Y_TYPE, _PAULI_Z_TYPE)
            images[_pauli_letter(T, q1, pauli_type)] = (1, _pauli_letter(T, q2, pauli_type))
            images[_pauli_letter(T, q2, pauli_type)] = (1, _pauli_letter(T, q1, pauli_type))
        end
        return CliffordSymmetry(images; nqubits=n)
    end

    throw(ArgumentError(
        "Unknown Clifford gate `:$gate`. Supported gates are :H, :S, :CNOT, and :SWAP."
    ))
end

function Base.show(io::IO, g::CliffordSymmetry)
    print(io, "CliffordSymmetry(nqubits=$(g.nqubits), images=$(length(g.images)))")
end

Base.:(==)(a::CliffordSymmetry, b::CliffordSymmetry) =
    a.nqubits == b.nqubits && a.images == b.images
Base.hash(g::CliffordSymmetry, h::UInt) = hash((g.nqubits, g.images), h)

function _phase_to_clifford_sign(phase::UInt8)
    phase_mod = phase % UInt8(4)
    phase_mod == UInt8(0) && return 1
    phase_mod == UInt8(2) && return -1
    throw(ArgumentError(
        "Clifford conjugation of a Hermitian Pauli word produced a non-real phase. " *
        "Check that the supplied images define a valid Pauli Clifford action."
    ))
end

function _convert_clifford_symmetry(
    ::Type{T},
    g::CliffordSymmetry;
    nqubits::Integer=g.nqubits,
) where {T<:Integer}
    n = _clifford_validate_nqubits(nqubits)
    g.nqubits <= n || throw(ArgumentError(
        "Cannot shrink CliffordSymmetry from nqubits=$(g.nqubits) to nqubits=$n."
    ))

    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    for (src, (sign, dst)) in g.images
        src_t = NormalMonomial{PauliAlgebra,T}(T.(src.word))
        dst_t = NormalMonomial{PauliAlgebra,T}(T.(dst.word))
        images[src_t] = (sign, dst_t)
    end
    return CliffordSymmetry(images; nqubits=n)
end

function _clifford_letter_image(
    g::CliffordSymmetry{T},
    idx::T,
) where {T<:Integer}
    letter = NormalMonomial{PauliAlgebra,T}(T[idx])
    return get(g.images, letter, (1, letter))
end

function _act_monomial(
    g::CliffordSymmetry{T},
    mono::NormalMonomial{PauliAlgebra,T},
) where {T<:Integer}
    direct = get(g.images, mono, nothing)
    isnothing(direct) || return direct

    sign = 1
    word = T[]
    for idx in mono.word
        letter_sign, image = _clifford_letter_image(g, idx)
        sign *= letter_sign
        append!(word, image.word)
    end

    image_word, phase = simplify(PauliAlgebra, word)
    sign *= _phase_to_clifford_sign(phase)
    return sign, NormalMonomial{PauliAlgebra,T}(image_word)
end

function _act_monomial(
    g::CliffordSymmetry,
    mono::NormalMonomial{PauliAlgebra,T},
) where {T<:Integer}
    return _act_monomial(_convert_clifford_symmetry(T, g), mono)
end

function _act_polynomial(
    g::CliffordSymmetry,
    poly::Polynomial{PauliAlgebra,T,C},
) where {T<:Integer,C<:Number}
    acted_terms = Tuple{C,NormalMonomial{PauliAlgebra,T}}[]
    sizehint!(acted_terms, length(poly.terms))
    for (coef, mono) in poly.terms
        sign, image = _act_monomial(g, mono)
        push!(acted_terms, (coef * sign, image))
    end
    return Polynomial(acted_terms)
end

function _clifford_letter_domain(::Type{T}, nqubits::Integer) where {T<:Integer}
    n = _clifford_validate_nqubits(nqubits)
    return T[convert(T, _pauli_index(site, pauli_type)) for site in 1:n for pauli_type in 0:2]
end

function _clifford_nqubits_from_domain(domain::AbstractVector{<:Integer})
    max_site = 0
    for idx in domain
        max_site = max(max_site, _pauli_site(idx))
    end
    return max_site
end

function _clifford_max_generator_nqubits(generators::AbstractVector{<:CliffordSymmetry})
    isempty(generators) && throw(ArgumentError("CliffordSymmetryGroup needs at least one generator."))
    return maximum(g.nqubits for g in generators)
end

function _compose_clifford_symmetries(
    left::CliffordSymmetry{T},
    right::CliffordSymmetry{T},
    domain::AbstractVector{T},
) where {T<:Integer}
    n = max(left.nqubits, right.nqubits, _clifford_nqubits_from_domain(domain))
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()

    for idx in domain
        src = NormalMonomial{PauliAlgebra,T}(T[idx])
        sign_right, mid = _act_monomial(right, src)
        sign_left, dst = _act_monomial(left, mid)
        sign = sign_left * sign_right
        if sign != 1 || dst != src
            images[src] = (sign, dst)
        end
    end

    return CliffordSymmetry(images; nqubits=n)
end

function _clifford_action_signature(g::CliffordSymmetry{T}, idx::T) where {T<:Integer}
    sign, image = _act_monomial(g, NormalMonomial{PauliAlgebra,T}(T[idx]))
    return (sign, Tuple(image.word))
end

function _clifford_push_site_letters!(domain::Set{T}, site::Integer) where {T<:Integer}
    for pauli_type in 0:2
        push!(domain, convert(T, _pauli_index(site, pauli_type)))
    end
    return domain
end

function _clifford_add_site_to_domain!(domain::Set{T}, queue::Vector{T}, site::Integer) where {T<:Integer}
    for pauli_type in 0:2
        idx = convert(T, _pauli_index(site, pauli_type))
        if !(idx in domain)
            push!(domain, idx)
            push!(queue, idx)
        end
    end
    return domain
end

# Full Clifford validation only needs sites touched by a nontrivial source or
# target image. Untouched sites are fixed pointwise, so checking all `nqubits`
# would turn sparse local gates into an accidental O(n²) constructor.
function _clifford_active_validation_domain(g::CliffordSymmetry{T}) where {T<:Integer}
    touched = Set{T}()
    sizehint!(touched, 3 * length(g.images))

    for (source, (_, image)) in g.images
        for idx in source.word
            _clifford_push_site_letters!(touched, _pauli_site(idx))
        end
        for idx in image.word
            _clifford_push_site_letters!(touched, _pauli_site(idx))
        end
    end

    return sort!(collect(touched))
end

function _clifford_action_closed_domain(
    generators::AbstractVector{CliffordSymmetry{T}},
    seeds::AbstractVector{T},
) where {T<:Integer}
    domain = Set{T}()
    queue = T[]

    for idx in seeds
        idx > 0 || throw(ArgumentError(
            "Clifford symmetry domain contains non-positive Pauli index $idx."
        ))
        _clifford_add_site_to_domain!(domain, queue, _pauli_site(idx))
    end

    for generator in generators
        for idx in _clifford_active_validation_domain(generator)
            _clifford_add_site_to_domain!(domain, queue, _pauli_site(idx))
        end
    end

    if isempty(domain)
        _clifford_add_site_to_domain!(domain, queue, 1)
    end

    head = 1
    while head <= length(queue)
        idx = queue[head]
        head += 1
        source = NormalMonomial{PauliAlgebra,T}(T[idx])
        for generator in generators
            _, image = _act_monomial(generator, source)
            for image_idx in image.word
                _clifford_add_site_to_domain!(domain, queue, _pauli_site(image_idx))
            end
        end
    end

    return sort!(collect(domain))
end

function _validate_clifford_symmetry(g::CliffordSymmetry{T}) where {T<:Integer}
    return _validate_clifford_symmetry(g, T[])
end

function _validate_clifford_symmetry(g::CliffordSymmetry{T}, domain::AbstractVector{T}) where {T<:Integer}
    required_nqubits = _clifford_nqubits_from_domain(domain)
    g.nqubits >= required_nqubits || throw(ArgumentError(
        "CliffordSymmetry acts on $(g.nqubits) qubits but the active domain reaches site $required_nqubits."
    ))
    _clifford_validate_domain_indices(domain, g.nqubits)

    if !isempty(domain)
        domain_set = Set(domain)
        for idx in domain
            _, image = _act_monomial(g, NormalMonomial{PauliAlgebra,T}(T[idx]))
            for image_idx in image.word
                image_idx in domain_set || throw(ArgumentError(
                    "CliffordSymmetry sends active Pauli index $idx outside the action-closed domain: $image_idx."
                ))
            end
        end
    end

    for (source, (sign, image)) in g.images
        sign ∈ (-1, 1) || throw(ArgumentError(
            "CliffordSymmetry signs must be ±1; got $sign for source word $source."
        ))
        _clifford_validate_pauli_indices(source, g.nqubits, "source word")
        _clifford_validate_pauli_indices(image, g.nqubits, "target word")
        isone(source) && throw(ArgumentError("CliffordSymmetry image dictionaries must not store the identity word."))
        length(source.word) == 1 || throw(ArgumentError(
            "CliffordSymmetry image dictionaries currently accept only single-Pauli-letter sources; got $source."
        ))
        !isone(image) || throw(ArgumentError(
            "CliffordSymmetry must map non-identity Pauli letters to non-identity Pauli words."
        ))
    end

    validation_domain = _clifford_active_validation_domain(g)
    isempty(validation_domain) && return nothing

    letter_images = Tuple{NormalMonomial{PauliAlgebra,T},NormalMonomial{PauliAlgebra,T}}[]
    image_words = Set{Tuple{Vararg{T}}}()
    sizehint!(letter_images, length(validation_domain))
    source_sites = Set{Int}()

    for idx in validation_domain
        source = NormalMonomial{PauliAlgebra,T}(T[idx])
        _clifford_validate_pauli_indices(source, g.nqubits, "source word")

        _, image = _act_monomial(g, source)
        !isone(image) || throw(ArgumentError(
            "CliffordSymmetry must map non-identity Pauli letters to non-identity Pauli words."
        ))
        _clifford_validate_pauli_indices(image, g.nqubits, "target word")

        word_key = Tuple(image.word)
        word_key in image_words && throw(ArgumentError(
            "CliffordSymmetry must map Pauli letters injectively up to sign on the active domain."
        ))
        push!(image_words, word_key)
        push!(letter_images, (source, image))
        haskey(g.images, source) && push!(source_sites, _pauli_site(idx))
    end

    for i in eachindex(letter_images), j in (i + 1):lastindex(letter_images)
        source_i, image_i = letter_images[i]
        source_j, image_j = letter_images[j]
        _pauli_anticommutes(source_i, source_j) == _pauli_anticommutes(image_i, image_j) || throw(ArgumentError(
            "CliffordSymmetry images must preserve Pauli commutation/anticommutation relations."
        ))
    end

    for site in source_sites
        for (left_type, right_type) in ((0, 1), (1, 2), (2, 0), (1, 0), (2, 1), (0, 2))
            _assert_clifford_preserves_local_product(g, site, left_type, right_type)
        end
    end

    return nothing
end

"""
    AbelianIrrepTable

Explicit Abelian irrep composition rules used by fermionic sector splitting.
"""
struct AbelianIrrepTable{Γ,M,D}
    identity::Γ
    multiply::M
    dual::D
end

function AbelianIrrepTable(identity::Γ; multiply, dual=(x -> x)) where {Γ}
    return AbelianIrrepTable{Γ,typeof(multiply),typeof(dual)}(identity, multiply, dual)
end

"""
    FermionicModeLayout

Explicit metadata for fermionic spin-orbitals used by the symmetry layer.
"""
struct FermionicModeLayout
    orbital_of::Dict{Int,Int}
    spin2_of::Dict{Int,Int}
    irrep_of::Dict{Int,Any}
    irrep_table::Union{Nothing,AbelianIrrepTable}
end

function FermionicModeLayout(
    orbital_of::AbstractDict{<:Integer,<:Integer}=Dict{Int,Int}();
    spin2_of::AbstractDict{<:Integer,<:Integer}=Dict{Int,Int}(),
    irrep_of::AbstractDict{<:Integer,<:Any}=Dict{Int,Any}(),
    irrep_table::Union{Nothing,AbelianIrrepTable}=nothing,
)
    orbitals = Dict{Int,Int}(Int(mode) => Int(orbital) for (mode, orbital) in orbital_of)
    spins = Dict{Int,Int}(Int(mode) => Int(spin2) for (mode, spin2) in spin2_of)
    irreps = Dict{Int,Any}(Int(orbital) => irrep for (orbital, irrep) in irrep_of)
    return FermionicModeLayout(orbitals, spins, irreps, irrep_table)
end

"""
    FermionicSectorSpec

Sector-splitting specification for fermionic symmetry reduction.
"""
@kwdef struct FermionicSectorSpec
    mode_layout::FermionicModeLayout = FermionicModeLayout()
    split_parity::Bool = true
    split_number::Bool = false
    split_spin::Bool = false
    split_irrep::Bool = false
end

"""
    FermionicSpinAdaptationSpec

Casimir-based SU(2) basis adaptation layered on top of fermionic sector splits.
"""
@kwdef struct FermionicSpinAdaptationSpec
    mode_layout::FermionicModeLayout = FermionicModeLayout()
    casimir_tol::Float64 = 1e-8
    compress_multiplets::Bool = false
end

"""
    FermionicSectorLabel

Quantum-number label for a fermionic basis element under the enabled sector split.
Disabled quantum numbers are stored as `nothing`.
"""
struct FermionicSectorLabel
    parity::Union{Nothing,Int}
    number::Union{Nothing,Int}
    spin2::Union{Nothing,Int}
    irrep
end

"""
    SymmetrySpec

Symmetry configuration for dense ordinary polynomial relaxations.

Supported pieces today:
- signed permutations on monoid-algebra registry indices;
- Clifford conjugation actions on Pauli words;
- fermionic mode permutations on physical modes;
- fermionic sector splitting by parity / particle number / `S_z` / Abelian irreps;
- fermionic SU(2) spin adaptation layered on top of fixed-sector blocks.

Unsupported combinations fail loudly instead of quietly doing nonsense.

# Soundness checks

- `check_invariance::Bool = true`: verify the objective and constraints are
  fixed by every group element before reducing.
- `offblock_check::Symbol = :randomized`: how to certify that the
  symmetry-adapted off-diagonal blocks of each constraint matrix vanish (the
  property that makes replacing one big PSD block by the small diagonal
  blocks exact):
  - `:randomized` (default) — evaluate every off-diagonal block at
    deterministic pseudo-random values of the orbit-reduced moment variables
    and check the result is numerically zero. Since each block entry is
    *linear* in those variables, vanishing at independent generic points
    certifies the symbolic property; two independent draws are used. Runs in
    BLAS time instead of polynomial arithmetic time.
  - `:full` — expand every off-diagonal block symbolically and assert each
    polynomial entry vanishes (the most expensive but entry-level diagnosis;
    use this to debug a failing `:randomized` certificate).
  - `:off` — skip the certificate entirely. Only safe when the generators are
    known-correct (e.g. re-running a spec that already passed).
"""
struct SymmetrySpec
    generators::Vector{SignedPermutation}
    fermionic_generators::Vector{FermionicModePermutation}
    clifford_generators::Vector{CliffordSymmetry}
    sector::Union{Nothing,FermionicSectorSpec}
    spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec}
    check_invariance::Bool
    offblock_check::Symbol

    function SymmetrySpec(
        generators::Vector{SignedPermutation},
        fermionic_generators::Vector{FermionicModePermutation},
        clifford_generators::Vector{CliffordSymmetry},
        sector::Union{Nothing,FermionicSectorSpec},
        spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec},
        check_invariance::Bool,
        offblock_check::Symbol=:randomized,
    )
        offblock_check in (:full, :randomized, :off) || throw(ArgumentError(
            "`offblock_check` must be `:full`, `:randomized`, or `:off`; got `$(repr(offblock_check))`."
        ))
        isempty(generators) && isempty(fermionic_generators) && isempty(clifford_generators) &&
            isnothing(sector) && isnothing(spin_adaptation) && throw(ArgumentError(
                "`SymmetrySpec` needs at least one finite action, fermionic sector split, or spin adaptation."
            ))

        n_finite_actions = count(!isempty, (generators, fermionic_generators, clifford_generators))
        n_finite_actions <= 1 || throw(ArgumentError(
            "`SymmetrySpec` does not yet support mixing `SignedPermutation`, `FermionicModePermutation`, and `CliffordSymmetry` actions in the same spec."
        ))

        !isempty(generators) && (!isnothing(sector) || !isnothing(spin_adaptation)) && throw(ArgumentError(
            "Fermionic sector splitting / spin adaptation cannot be combined with raw `SignedPermutation` symmetry. Use `FermionicModePermutation` for fermionic problems."
        ))
        !isempty(clifford_generators) && (!isnothing(sector) || !isnothing(spin_adaptation)) && throw(ArgumentError(
            "Fermionic sector splitting / spin adaptation cannot be combined with `CliffordSymmetry` symmetry."
        ))

        return new(generators, fermionic_generators, clifford_generators, sector, spin_adaptation, check_invariance, offblock_check)
    end
end

function SymmetrySpec(
    generators::AbstractVector{<:SignedPermutation};
    check_invariance::Bool=true,
    offblock_check::Symbol=:randomized,
)
    isempty(generators) && throw(ArgumentError("`SymmetrySpec` needs at least one generator."))
    converted = SignedPermutation[generator for generator in generators]
    return SymmetrySpec(converted, FermionicModePermutation[], CliffordSymmetry[], nothing, nothing, check_invariance, offblock_check)
end

function SymmetrySpec(generators::SignedPermutation...; check_invariance::Bool=true, offblock_check::Symbol=:randomized)
    return SymmetrySpec(collect(generators); check_invariance, offblock_check)
end

function SymmetrySpec(
    fermionic_generators::AbstractVector{<:FermionicModePermutation};
    sector::Union{Nothing,FermionicSectorSpec}=nothing,
    spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec}=nothing,
    check_invariance::Bool=true,
    offblock_check::Symbol=:randomized,
)
    isempty(fermionic_generators) && isnothing(sector) && isnothing(spin_adaptation) && throw(ArgumentError(
        "`SymmetrySpec` needs at least one fermionic mode permutation, sector split, or spin adaptation."
    ))
    converted = FermionicModePermutation[generator for generator in fermionic_generators]
    return SymmetrySpec(SignedPermutation[], converted, CliffordSymmetry[], sector, spin_adaptation, check_invariance, offblock_check)
end

function SymmetrySpec(
    fermionic_generators::FermionicModePermutation...;
    sector::Union{Nothing,FermionicSectorSpec}=nothing,
    spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec}=nothing,
    check_invariance::Bool=true,
    offblock_check::Symbol=:randomized,
)
    return SymmetrySpec(
        collect(fermionic_generators);
        sector=sector,
        spin_adaptation=spin_adaptation,
        check_invariance=check_invariance,
        offblock_check=offblock_check,
    )
end

function SymmetrySpec(
    clifford_generators::AbstractVector{<:CliffordSymmetry};
    check_invariance::Bool=true,
    offblock_check::Symbol=:randomized,
)
    isempty(clifford_generators) && throw(ArgumentError("`SymmetrySpec` needs at least one Clifford generator."))
    converted = CliffordSymmetry[generator for generator in clifford_generators]
    return SymmetrySpec(SignedPermutation[], FermionicModePermutation[], converted, nothing, nothing, check_invariance, offblock_check)
end

function SymmetrySpec(clifford_generators::CliffordSymmetry...; check_invariance::Bool=true, offblock_check::Symbol=:randomized)
    return SymmetrySpec(collect(clifford_generators); check_invariance, offblock_check)
end

function SymmetrySpec(
    ;
    generators::AbstractVector{<:SignedPermutation}=SignedPermutation[],
    fermionic_generators::AbstractVector{<:FermionicModePermutation}=FermionicModePermutation[],
    clifford_generators::AbstractVector{<:CliffordSymmetry}=CliffordSymmetry[],
    sector::Union{Nothing,FermionicSectorSpec}=nothing,
    spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec}=nothing,
    check_invariance::Bool=true,
    offblock_check::Symbol=:randomized,
)
    return SymmetrySpec(
        SignedPermutation[generator for generator in generators],
        FermionicModePermutation[generator for generator in fermionic_generators],
        CliffordSymmetry[generator for generator in clifford_generators],
        sector,
        spin_adaptation,
        check_invariance,
        offblock_check,
    )
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
    block_labels::Vector{Any}
    block_provenance::Vector{Symbol}
end

function SymmetryReport(
    group_order::Int,
    invariant_moment_count::Int,
    psd_block_sizes::Vector{Int},
    basis_half_size::Int,
    basis_full_size::Int,
)
    return SymmetryReport(
        group_order,
        invariant_moment_count,
        psd_block_sizes,
        basis_half_size,
        basis_full_size,
        Any[],
        Symbol[],
    )
end

function Base.show(io::IO, report::SymmetryReport)
    print(
        io,
        "SymmetryReport(group_order=$(report.group_order), " *
        "invariant_moment_count=$(report.invariant_moment_count), " *
        "psd_block_sizes=$(report.psd_block_sizes), " *
        "basis_half_size=$(report.basis_half_size), " *
        "basis_full_size=$(report.basis_full_size), " *
        "block_provenance=$(report.block_provenance))"
    )
end

struct _OrbitReducer{A<:AlgebraType,T<:Integer}
    representative::Dict{NormalMonomial{A,T},NormalMonomial{A,T}}
    phase::Dict{NormalMonomial{A,T},Int}
end

struct NCWordSignedPermutationAction{A<:MonoidAlgebra,T<:Integer} <: SymbolicWedderburn.BySignedPermutations
    algebra::Type{A}
end

struct NCFermionicModePermutationAction{T<:Integer} <: SymbolicWedderburn.BySignedPermutations end

struct NCPauliCliffordAction{T<:Integer} <: SymbolicWedderburn.BySignedPermutations end

struct _FiniteSymmetryGroup{G,T<:Integer} <: GroupsCore.Group
    domain::Vector{T}
    elements::Vector{G}
    generator_indices::Vector{Int}
    mult_table::Matrix{Int}
    inv_table::Vector{Int}
    element_orders::Vector{Int}
end

struct _FiniteSymmetryGroupElement{G,T<:Integer} <: GroupsCore.GroupElement
    group::_FiniteSymmetryGroup{G,T}
    index::Int
end

const _SignedPermutationGroup{T} = _FiniteSymmetryGroup{SignedPermutation{T},T}
const _SignedPermutationGroupElement{T} = _FiniteSymmetryGroupElement{SignedPermutation{T},T}
const _FermionicModePermutationGroup{T} = _FiniteSymmetryGroup{FermionicModePermutation{T},T}
const _FermionicModePermutationGroupElement{T} = _FiniteSymmetryGroupElement{FermionicModePermutation{T},T}

"""
    CliffordSymmetryGroup(generators...; nqubits=nothing, integer_type=Int, domain=nothing)

Finite `GroupsCore` group generated by [`CliffordSymmetry`](@ref) actions.
Most users should pass `SymmetrySpec(CliffordSymmetry(...))` through
`SolverConfig`; this constructor is for direct action tests and advanced use.
Pass `domain` to restrict enumeration to the action-closed Pauli-letter domain
needed by a local high-index action.
"""
struct CliffordSymmetryGroup{T<:Integer} <: GroupsCore.Group
    inner::_FiniteSymmetryGroup{CliffordSymmetry{T},T}
end

struct CliffordSymmetryGroupElement{T<:Integer} <: GroupsCore.GroupElement
    group::CliffordSymmetryGroup{T}
    index::Int
end

function CliffordSymmetryGroup(
    generators::AbstractVector{<:CliffordSymmetry};
    nqubits=nothing,
    integer_type::Type{T}=Int,
    domain=nothing,
) where {T<:Integer}
    generator_nqubits = _clifford_max_generator_nqubits(generators)
    seed_nqubits = isnothing(domain) ? 0 : _clifford_nqubits_from_domain(collect(domain))
    inferred_nqubits = max(generator_nqubits, seed_nqubits)
    n = isnothing(nqubits) ? inferred_nqubits : _clifford_validate_nqubits(nqubits)
    inferred_nqubits <= n || throw(ArgumentError(
        "CliffordSymmetryGroup cannot shrink generators/domain acting on $inferred_nqubits qubits to nqubits=$n."
    ))

    converted = CliffordSymmetry{T}[
        _convert_clifford_symmetry(T, generator; nqubits=n) for generator in generators
    ]
    action_domain = if isnothing(domain)
        _clifford_letter_domain(T, n)
    else
        seeds = T[convert(T, idx) for idx in domain]
        _clifford_action_closed_domain(converted, seeds)
    end

    group = _enumerate_symmetry_group(converted, action_domain)
    return CliffordSymmetryGroup{T}(_sw_finite_action_group(converted, group, action_domain))
end

function CliffordSymmetryGroup(
    generators::CliffordSymmetry...;
    nqubits=nothing,
    integer_type::Type{T}=Int,
    domain=nothing,
) where {T<:Integer}
    return CliffordSymmetryGroup(collect(generators); nqubits, integer_type=T, domain)
end

@inline _clifford_group_value(g::CliffordSymmetryGroupElement{T}) where {T<:Integer} =
    g.group.inner.elements[g.index]

Base.eltype(::Type{CliffordSymmetryGroup{T}}) where {T<:Integer} = CliffordSymmetryGroupElement{T}
Base.IteratorSize(::Type{<:CliffordSymmetryGroup}) = Base.HasLength()
Base.length(G::CliffordSymmetryGroup) = length(G.inner)

function Base.iterate(G::CliffordSymmetryGroup{T}, state::Int=1) where {T<:Integer}
    state > length(G) && return nothing
    return CliffordSymmetryGroupElement{T}(G, state), state + 1
end

Base.one(G::CliffordSymmetryGroup{T}) where {T<:Integer} = CliffordSymmetryGroupElement{T}(G, 1)
Base.parent(g::CliffordSymmetryGroupElement) = g.group
Base.copy(g::CliffordSymmetryGroupElement) = g
Base.isone(g::CliffordSymmetryGroupElement) = g.index == 1

function Base.:(==)(a::CliffordSymmetryGroupElement{T}, b::CliffordSymmetryGroupElement{T}) where {T<:Integer}
    return a.group.inner.domain == b.group.inner.domain && _clifford_group_value(a) == _clifford_group_value(b)
end

function Base.hash(g::CliffordSymmetryGroupElement, h::UInt)
    return hash((g.group.inner.domain, _clifford_group_value(g)), h)
end

function GroupsCore.order(::Type{I}, G::CliffordSymmetryGroup) where {I<:Integer}
    return convert(I, length(G))
end

function GroupsCore.order(::Type{I}, g::CliffordSymmetryGroupElement) where {I<:Integer}
    return convert(I, g.group.inner.element_orders[g.index])
end

function GroupsCore.gens(G::CliffordSymmetryGroup{T}) where {T<:Integer}
    return [CliffordSymmetryGroupElement{T}(G, idx) for idx in G.inner.generator_indices]
end

function Base.inv(g::CliffordSymmetryGroupElement{T}) where {T<:Integer}
    return CliffordSymmetryGroupElement{T}(g.group, g.group.inner.inv_table[g.index])
end

function Base.:(*)(a::CliffordSymmetryGroupElement{T}, b::CliffordSymmetryGroupElement{T}) where {T<:Integer}
    a.group.inner.domain == b.group.inner.domain && a.group.inner.elements == b.group.inner.elements || throw(ArgumentError(
        "Cannot multiply elements from different CliffordSymmetryGroups."
    ))
    return CliffordSymmetryGroupElement{T}(a.group, a.group.inner.mult_table[a.index, b.index])
end

@inline _signed_image(g::SignedPermutation{T}, idx::T) where {T<:Integer} = get(g.images, idx, (1, idx))
@inline _signed_image_signature(g::SignedPermutation{T}, idx::T) where {T<:Integer} = _signed_image(g, idx)
@inline _mode_image(g::FermionicModePermutation{T}, idx::T) where {T<:Integer} = get(g.images, idx, idx)

function _convert_signed_permutation(::Type{T}, g::SignedPermutation) where {T<:Integer}
    images = Dict{T,Tuple{Int,T}}()
    for (src, (sign, dst)) in g.images
        images[convert(T, src)] = (sign, convert(T, dst))
    end
    return SignedPermutation(images)
end

function _convert_fermionic_mode_permutation(::Type{T}, g::FermionicModePermutation) where {T<:Integer}
    images = Dict{T,T}()
    for (src, dst) in g.images
        images[convert(T, src)] = convert(T, dst)
    end
    return FermionicModePermutation(images)
end

function _convert_symmetry_spec(::Type{T}, spec::SymmetrySpec) where {T<:Integer}
    converted = SignedPermutation[_convert_signed_permutation(T, g) for g in spec.generators]
    fermionic = FermionicModePermutation[_convert_fermionic_mode_permutation(Int, g) for g in spec.fermionic_generators]
    clifford = CliffordSymmetry[_convert_clifford_symmetry(T, g) for g in spec.clifford_generators]
    return SymmetrySpec(
        generators=converted,
        fermionic_generators=fermionic,
        clifford_generators=clifford,
        sector=spec.sector,
        spin_adaptation=spec.spin_adaptation,
        check_invariance=spec.check_invariance,
    )
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

function _validate_fermionic_mode_permutation(
    g::FermionicModePermutation{T},
    domain::AbstractVector{T},
) where {T<:Integer}
    domain_set = Set(domain)
    targets = T[]
    sizehint!(targets, length(domain))

    for idx in domain
        dst = _mode_image(g, idx)
        dst in domain_set || throw(ArgumentError(
            "FermionicModePermutation sends mode $idx outside the active symmetry domain: $dst."
        ))
        push!(targets, dst)
    end

    length(unique(targets)) == length(domain) || throw(ArgumentError(
        "FermionicModePermutation must be bijective on the active symmetry domain."
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

function _compose_fermionic_mode_permutations(
    left::FermionicModePermutation{T},
    right::FermionicModePermutation{T},
    domain::AbstractVector{T},
) where {T<:Integer}
    images = Dict{T,T}()
    for idx in domain
        dst = _mode_image(left, _mode_image(right, idx))
        dst == idx || (images[idx] = dst)
    end
    return FermionicModePermutation(images)
end

@inline _identity_signed_permutation(::Type{T}) where {T<:Integer} =
    SignedPermutation(Dict{T,Tuple{Int,T}}())
@inline _identity_fermionic_mode_permutation(::Type{T}) where {T<:Integer} =
    FermionicModePermutation(Dict{T,T}())

function _symmetry_key(g::SignedPermutation{T}, domain::AbstractVector{T}) where {T<:Integer}
    return Tuple(_signed_image_signature(g, idx) for idx in domain)
end

function _symmetry_key(g::FermionicModePermutation{T}, domain::AbstractVector{T}) where {T<:Integer}
    return Tuple(_mode_image(g, idx) for idx in domain)
end

function _symmetry_key(g::CliffordSymmetry{T}, domain::AbstractVector{T}) where {T<:Integer}
    return Tuple(_clifford_action_signature(g, idx) for idx in domain)
end

_validate_finite_action(g::SignedPermutation{T}, domain::AbstractVector{T}) where {T<:Integer} =
    _validate_signed_permutation(g, domain)
_validate_finite_action(g::FermionicModePermutation{T}, domain::AbstractVector{T}) where {T<:Integer} =
    _validate_fermionic_mode_permutation(g, domain)
_validate_finite_action(g::CliffordSymmetry{T}, domain::AbstractVector{T}) where {T<:Integer} =
    _validate_clifford_symmetry(g, domain)

_compose_finite_actions(left::SignedPermutation{T}, right::SignedPermutation{T}, domain::AbstractVector{T}) where {T<:Integer} =
    _compose_signed_permutations(left, right, domain)
_compose_finite_actions(left::FermionicModePermutation{T}, right::FermionicModePermutation{T}, domain::AbstractVector{T}) where {T<:Integer} =
    _compose_fermionic_mode_permutations(left, right, domain)
_compose_finite_actions(left::CliffordSymmetry{T}, right::CliffordSymmetry{T}, domain::AbstractVector{T}) where {T<:Integer} =
    _compose_clifford_symmetries(left, right, domain)

_inverse_finite_action(g::SignedPermutation{T}, domain::AbstractVector{T}) where {T<:Integer} = begin
    images = Dict{T,Tuple{Int,T}}()
    for idx in domain
        sign, dst = _signed_image(g, idx)
        if sign != 1 || idx != dst
            images[dst] = (sign, idx)
        end
    end
    return SignedPermutation(images)
end

_inverse_finite_action(g::FermionicModePermutation{T}, domain::AbstractVector{T}) where {T<:Integer} = begin
    images = Dict{T,T}()
    for idx in domain
        dst = _mode_image(g, idx)
        idx == dst || (images[dst] = idx)
    end
    return FermionicModePermutation(images)
end

_identity_finite_action(::Type{SignedPermutation{T}}) where {T<:Integer} = _identity_signed_permutation(T)
_identity_finite_action(::Type{FermionicModePermutation{T}}) where {T<:Integer} = _identity_fermionic_mode_permutation(T)
function _identity_finite_action(::Type{CliffordSymmetry{T}}, domain::AbstractVector{T}) where {T<:Integer}
    n = max(1, _clifford_nqubits_from_domain(domain))
    return CliffordSymmetry(Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}(); nqubits=n)
end
_identity_finite_action(::Type{G}, _domain::AbstractVector) where {G} = _identity_finite_action(G)

function _enumerate_symmetry_group(spec::SymmetrySpec, domain::Vector{T}) where {T<:Integer}
    isempty(domain) && throw(ArgumentError("Symmetry reduction needs a non-empty active variable domain."))

    if !isempty(spec.generators)
        generators = SignedPermutation{T}[_convert_signed_permutation(T, g) for g in spec.generators]
        return _enumerate_symmetry_group(generators, domain)
    elseif !isempty(spec.fermionic_generators)
        generators = FermionicModePermutation{T}[
            _convert_fermionic_mode_permutation(T, g) for g in spec.fermionic_generators
        ]
        return _enumerate_symmetry_group(generators, domain)
    elseif !isempty(spec.clifford_generators)
        nqubits = max(_clifford_nqubits_from_domain(domain), _clifford_max_generator_nqubits(spec.clifford_generators))
        generators = CliffordSymmetry{T}[
            _convert_clifford_symmetry(T, g; nqubits) for g in spec.clifford_generators
        ]
        clifford_domain = _clifford_action_closed_domain(generators, domain)
        return _enumerate_symmetry_group(generators, clifford_domain)
    end

    throw(ArgumentError("`SymmetrySpec` does not contain any finite symmetry generators."))
end

function _enumerate_symmetry_group(generators::Vector{G}, domain::Vector{T}) where {G,T<:Integer}
    isempty(generators) && throw(ArgumentError("Symmetry reduction needs at least one finite generator."))

    for generator in generators
        _validate_finite_action(generator, domain)
    end

    identity = _identity_finite_action(G, domain)
    seen = Dict{Any,G}(_symmetry_key(identity, domain) => identity)
    queue = G[identity]

    while !isempty(queue)
        current = popfirst!(queue)
        for generator in generators
            candidate = _compose_finite_actions(generator, current, domain)
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

@inline _sw_group_value(g::_FiniteSymmetryGroupElement{G,T}) where {G,T<:Integer} = g.group.elements[g.index]

Base.eltype(::Type{_FiniteSymmetryGroup{G,T}}) where {G,T<:Integer} = _FiniteSymmetryGroupElement{G,T}
Base.IteratorSize(::Type{<:_FiniteSymmetryGroup}) = Base.HasLength()
Base.length(G::_FiniteSymmetryGroup) = length(G.elements)

function Base.iterate(G::_FiniteSymmetryGroup{GAction,T}, state::Int=1) where {GAction,T<:Integer}
    state > length(G.elements) && return nothing
    return _FiniteSymmetryGroupElement{GAction,T}(G, state), state + 1
end

Base.one(G::_FiniteSymmetryGroup{GAction,T}) where {GAction,T<:Integer} =
    _FiniteSymmetryGroupElement{GAction,T}(G, 1)
Base.parent(g::_FiniteSymmetryGroupElement) = g.group
Base.copy(g::_FiniteSymmetryGroupElement) = g
Base.isone(g::_FiniteSymmetryGroupElement) = g.index == 1

function Base.:(==)(a::_FiniteSymmetryGroupElement{G,T}, b::_FiniteSymmetryGroupElement{G,T}) where {G,T<:Integer}
    return a.group.domain == b.group.domain && _sw_group_value(a) == _sw_group_value(b)
end

function Base.hash(g::_FiniteSymmetryGroupElement, h::UInt)
    return hash((g.group.domain, _sw_group_value(g)), h)
end

function GroupsCore.order(::Type{I}, G::_FiniteSymmetryGroup) where {I<:Integer}
    return convert(I, length(G.elements))
end

function GroupsCore.order(::Type{I}, g::_FiniteSymmetryGroupElement) where {I<:Integer}
    return convert(I, g.group.element_orders[g.index])
end

function GroupsCore.gens(G::_FiniteSymmetryGroup{GAction,T}) where {GAction,T<:Integer}
    return [_FiniteSymmetryGroupElement{GAction,T}(G, idx) for idx in G.generator_indices]
end

function Base.inv(g::_FiniteSymmetryGroupElement{GAction,T}) where {GAction,T<:Integer}
    return _FiniteSymmetryGroupElement{GAction,T}(g.group, g.group.inv_table[g.index])
end

function Base.:(*)(a::_FiniteSymmetryGroupElement{GAction,T}, b::_FiniteSymmetryGroupElement{GAction,T}) where {GAction,T<:Integer}
    a.group.domain == b.group.domain && a.group.elements == b.group.elements || throw(ArgumentError(
        "Cannot multiply elements from different finite symmetry groups."
    ))
    return _FiniteSymmetryGroupElement{GAction,T}(a.group, a.group.mult_table[a.index, b.index])
end

function _sw_finite_action_group(
    generators::Vector{G},
    group::Vector{G},
    domain::Vector{T},
) where {G,T<:Integer}
    identity = _identity_finite_action(G, domain)
    identity_key = _symmetry_key(identity, domain)
    elements = G[identity]
    for g in group
        _symmetry_key(g, domain) == identity_key && continue
        push!(elements, g)
    end

    lookup = Dict{Any,Int}(_symmetry_key(g, domain) => i for (i, g) in enumerate(elements))
    generator_indices = Int[]
    for generator in generators
        idx = get(lookup, _symmetry_key(generator, domain), 0)
        idx == 0 && throw(ArgumentError(
            "Internal symmetry error: failed to locate a generator in the enumerated symmetry group."
        ))
        push!(generator_indices, idx)
    end

    n = length(elements)
    mult_table = Matrix{Int}(undef, n, n)
    for j in 1:n, i in 1:n
        product = _compose_finite_actions(elements[j], elements[i], domain)
        idx = get(lookup, _symmetry_key(product, domain), 0)
        idx == 0 && throw(ArgumentError(
            "Internal symmetry error: enumerated symmetry group is not closed under multiplication."
        ))
        mult_table[i, j] = idx
    end

    inv_table = Vector{Int}(undef, n)
    for i in 1:n
        inverse_idx = findfirst(j -> mult_table[i, j] == 1 && mult_table[j, i] == 1, 1:n)
        isnothing(inverse_idx) && throw(ArgumentError(
            "Internal symmetry error: enumerated symmetry group is missing an inverse element."
        ))
        inv_table[i] = inverse_idx
    end

    element_orders = ones(Int, n)
    for i in 2:n
        order_i = 1
        power = i
        while power != 1
            power = mult_table[power, i]
            order_i += 1
            order_i <= n || throw(ArgumentError(
                "Internal symmetry error: failed to compute an element order for the finite symmetry group."
            ))
        end
        element_orders[i] = order_i
    end

    return _FiniteSymmetryGroup(domain, elements, generator_indices, mult_table, inv_table, element_orders)
end

function _sw_signed_permutation_group(
    spec::SymmetrySpec,
    group::Vector{SignedPermutation{T}},
    domain::Vector{T},
) where {T<:Integer}
    generators = SignedPermutation{T}[_convert_signed_permutation(T, g) for g in spec.generators]
    return _sw_finite_action_group(generators, group, domain)
end

function _sw_fermionic_mode_permutation_group(
    spec::SymmetrySpec,
    group::Vector{FermionicModePermutation{T}},
    domain::Vector{T},
) where {T<:Integer}
    generators = FermionicModePermutation{T}[
        _convert_fermionic_mode_permutation(T, g) for g in spec.fermionic_generators
    ]
    return _sw_finite_action_group(generators, group, domain)
end

function SymbolicWedderburn.action(
    ::NCWordSignedPermutationAction{A,T},
    g::_SignedPermutationGroupElement{T},
    mono::NormalMonomial{A,T},
) where {A<:MonoidAlgebra,T<:Integer}
    sign, image = _act_monomial(_sw_group_value(g), mono)
    return image, sign
end

function SymbolicWedderburn.action(
    ::NCFermionicModePermutationAction{T},
    g::_FermionicModePermutationGroupElement{Int},
    mono::NormalMonomial{FermionicAlgebra,T},
) where {T<:Integer}
    sign, image = _act_monomial(_sw_group_value(g), mono)
    return image, sign
end

function SymbolicWedderburn.action(
    ::NCPauliCliffordAction{T},
    g::CliffordSymmetryGroupElement{T},
    mono::NormalMonomial{PauliAlgebra,T},
) where {T<:Integer}
    sign, image = _act_monomial(_clifford_group_value(g), mono)
    return image, sign
end

function _sw_decompose_half_basis(
    basis::Vector{M},
    group::_SignedPermutationGroup{T},
) where {A<:MonoidAlgebra,T<:Integer,M<:NormalMonomial{A,T}}
    action = NCWordSignedPermutationAction{A,T}(A)
    return SymbolicWedderburn.symmetry_adapted_basis(Float64, group, action, basis)
end

function _sw_decompose_half_basis(
    basis::Vector{NormalMonomial{FermionicAlgebra,T}},
    group::_FermionicModePermutationGroup{Int},
) where {T<:Integer}
    action = NCFermionicModePermutationAction{T}()
    return SymbolicWedderburn.symmetry_adapted_basis(Float64, group, action, basis)
end

function _sw_decompose_half_basis(
    basis::Vector{NormalMonomial{PauliAlgebra,T}},
    group::CliffordSymmetryGroup{T},
) where {T<:Integer}
    action = NCPauliCliffordAction{T}()
    return SymbolicWedderburn.symmetry_adapted_basis(Float64, group, action, basis)
end

@inline _sw_action_elements(group::_FiniteSymmetryGroup) = group.elements
@inline _sw_action_elements(group::CliffordSymmetryGroup) = group.inner.elements

function _basis_action_is_trivial(basis::AbstractVector{M}, sw_group) where {M<:NormalMonomial}
    for g in _sw_action_elements(sw_group), mono in basis
        sign, image = _act_monomial(g, mono)
        (sign == 1 && image == mono) || return false
    end
    return true
end

function _symmetry_support(poly::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    used = Set{T}()
    for (_, mono) in poly.terms
        union!(used, _symmetry_support(mono))
    end
    return used
end

function _symmetry_support(mono::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    used = Set{T}()
    for idx in mono.word
        push!(used, idx)
    end
    return used
end

function _symmetry_domain(
    pop::PolyOpt{A,T,P},
    corr_sparsity::CorrelativeSparsity{A,T,P,M,Nothing},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C},M<:NormalMonomial{A,T}}
    used = Set{T}()

    for poly in (pop.objective,)
        union!(used, _symmetry_support(poly))
    end
    for poly in pop.eq_constraints
        union!(used, _symmetry_support(poly))
    end
    for poly in pop.ineq_constraints
        union!(used, _symmetry_support(poly))
    end
    for poly in pop.moment_eq_constraints
        union!(used, _symmetry_support(poly))
    end

    for basis in corr_sparsity.clq_mom_mtx_bases
        for mono in basis
            union!(used, _symmetry_support(mono))
        end
    end
    for term_sparsities in cliques_term_sparsities, term_sparsity in term_sparsities, block_basis in term_sparsity.block_bases, mono in block_basis
        union!(used, _symmetry_support(mono))
    end

    return sort!(collect(used))
end

function _mode_symmetry_domain(
    pop::PolyOpt{A,T,P},
    corr_sparsity::CorrelativeSparsity{A,T,P,M,Nothing},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C},M<:NormalMonomial{A,T}}
    raw_domain = _symmetry_domain(pop, corr_sparsity, cliques_term_sparsities)
    used = unique!(sort!([Int(abs(idx)) for idx in raw_domain]))
    return used
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

function _permute_fermionic_word(
    g::FermionicModePermutation,
    word::Vector{T},
) where {T<:Integer}
    acted = similar(word)
    for (i, idx) in pairs(word)
        mode = Int(abs(idx))
        dst = _mode_image(g, mode)
        acted[i] = idx < 0 ? -T(dst) : T(dst)
    end
    return acted
end

function _act_monomial(
    g::FermionicModePermutation,
    mono::NormalMonomial{FermionicAlgebra,T},
) where {T<:Integer}
    terms = _simplified_to_terms(
        FermionicAlgebra,
        simplify(FermionicAlgebra, _permute_fermionic_word(g, mono.word)),
        T,
    )
    length(terms) == 1 || throw(ArgumentError(
        "FermionicModePermutation must preserve monomials after normal ordering; got $(length(terms)) terms for basis element $mono."
    ))

    coef, image = only(terms)
    coef == one(coef) && return 1, image
    coef == -one(coef) && return -1, image
    throw(ArgumentError(
        "FermionicModePermutation should induce only a sign on a basis monomial, got coefficient $coef for $mono."
    ))
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

function _act_polynomial(
    g::FermionicModePermutation,
    poly::Polynomial{FermionicAlgebra,T,C},
) where {T<:Integer,C<:Number}
    acted_terms = Tuple{C,NormalMonomial{FermionicAlgebra,T}}[]
    for (coef, mono) in poly.terms
        simplified = simplify(FermionicAlgebra, _permute_fermionic_word(g, mono.word))
        for (term_coef, image) in _simplified_to_terms(FermionicAlgebra, simplified, T)
            push!(acted_terms, (coef * term_coef, image))
        end
    end
    return Polynomial(acted_terms)
end

function _check_polynomial_invariance(
    label::AbstractString,
    poly::Polynomial{A,T,C},
    group::AbstractVector,
) where {A<:AlgebraType,T<:Integer,C<:Number}
    for g in group
        _act_polynomial(g, poly) == poly || throw(ArgumentError(
            "Symmetry reduction requires $label to be invariant under the supplied symmetry."
        ))
    end
    return nothing
end

function _check_symmetry_invariance(
    pop::PolyOpt{A,T,P},
    group::AbstractVector,
) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C}}
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
    group::AbstractVector,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T}}
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

@inline function _symmetric_monomial(mono::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    return NormalMonomial{A,T}(symmetric_canon(mono))
end

function _build_orbit_reducer(
    basis::AbstractVector{<:NormalMonomial{A,T}},
    group::AbstractVector,
) where {A<:AlgebraType,T<:Integer}
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
) where {A<:AlgebraType,T<:Integer,C<:Number}
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
) where {A<:AlgebraType,T<:Integer,C<:Number}
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
) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C}}
    reduced = Matrix{P}(undef, size(mat, 1), size(mat, 2))
    for j in axes(mat, 2), i in axes(mat, 1)
        reduced[i, j] = _orbit_reduce_polynomial(_chop_polynomial(mat[i, j]; atol), reducer; atol)
    end
    return reduced
end

function _poly_is_small(
    poly::Polynomial{A,T,C};
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:AlgebraType,T<:Integer,C<:Number}
    return iszero(_chop_polynomial(poly; atol))
end

function _assert_poly_small(
    poly::Polynomial{A,T,C},
    context::AbstractString;
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:AlgebraType,T<:Integer,C<:Number}
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
) where {A<:AlgebraType,T<:Integer,C<:Number}
    diff = _chop_polynomial(lhs - rhs; atol)
    iszero(diff) || throw(ArgumentError(
        "$context should match after symmetry reduction, but differs by $diff."
    ))
    return nothing
end

struct _BasisPartitionBlock{L}
    indices::Vector{Int}
    label::L
    provenance::Symbol
end

struct _BasisTransformBlock{L,M<:AbstractMatrix}
    row_basis::M
    label::L
    provenance::Symbol
end

function _chop_polynomial_matrix(mat::Matrix{P}; atol::Float64=_SYMMETRY_ATOL) where {P<:Polynomial}
    chopped = similar(mat)
    for j in axes(mat, 2), i in axes(mat, 1)
        chopped[i, j] = _chop_polynomial(mat[i, j]; atol)
    end
    return chopped
end

@inline _maybe_orbit_reduce_matrix(mat, ::Nothing; atol::Float64=_SYMMETRY_ATOL) = _chop_polynomial_matrix(mat; atol)
@inline _maybe_orbit_reduce_matrix(mat, reducer::_OrbitReducer; atol::Float64=_SYMMETRY_ATOL) =
    _orbit_reduce_matrix(mat, reducer; atol)

function _transform_polynomial_block(
    mat::Matrix{P},
    U_left::AbstractMatrix{<:Number},
    U_right::AbstractMatrix{<:Number};
    scale::Number=1,
    atol::Float64=_SYMMETRY_ATOL,
) where {P<:Polynomial}
    transformed = Matrix{P}(undef, size(U_left, 1), size(U_right, 1))

    for j in axes(transformed, 2), i in axes(transformed, 1)
        acc = zero(P)
        for b in axes(mat, 2), a in axes(mat, 1)
            coeff = U_left[i, a] * conj(U_right[j, b])
            abs(coeff) <= atol && continue
            acc += coeff * mat[a, b]
        end
        transformed[i, j] = _chop_polynomial(scale == 1 ? acc : scale * acc; atol)
    end

    return transformed
end

function _diagonal_transformed_blocks(
    mat::Matrix{P},
    row_bases::Vector{<:AbstractMatrix},
    reducer::Union{Nothing,_OrbitReducer};
    atol::Float64=_SYMMETRY_ATOL,
) where {P<:Polynomial}
    diagonal_blocks = Matrix{P}[]
    for row_basis in row_bases
        diagonal = _maybe_orbit_reduce_matrix(
            _transform_polynomial_block(mat, row_basis, row_basis; atol),
            reducer;
            atol,
        )
        push!(diagonal_blocks, _chop_polynomial_matrix(diagonal; atol))
    end
    return diagonal_blocks
end

function _append_transformed_offblock_zero_constraints!(
    constraints::Vector{Tuple{Symbol, Matrix{P}}},
    mat::Matrix{P},
    row_bases::Vector{<:AbstractMatrix},
    reducer::Union{Nothing,_OrbitReducer};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="constraint matrix",
) where {P<:Polynomial}
    for i in 1:length(row_bases), j in (i+1):length(row_bases)
        for (left, right) in ((i, j), (j, i))
            off_block = _maybe_orbit_reduce_matrix(
                _transform_polynomial_block(mat, row_bases[left], row_bases[right]; atol),
                reducer;
                atol,
            )
            for col in axes(off_block, 2), row in axes(off_block, 1)
                poly = _chop_polynomial(off_block[row, col]; atol)
                iszero(poly) || _append_constraint!(constraints, :Zero, reshape([poly], 1, 1), P)
            end
        end
    end
    return nothing
end

# SplitMix64 mixer: deterministic, dependency-free source of generic weights
# for the randomized off-block certificate. Pure integer ops, so the weights
# are reproducible across runs and platforms.
@inline function _splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

# Deterministic pseudo-random weight in (-1, 1) for variable `idx` under `seed`.
@inline function _pseudorandom_weight(idx::Int, seed::UInt64)
    bits = _splitmix64(seed ⊻ _splitmix64(UInt64(idx)))
    return Float64(bits >> 11) / 9.007199254740992e15 * 2.0 - 1.0
end

"""
    _verify_offblocks_randomized(mat, row_bases, reducer; atol, context, ntrials)

Certify that all off-diagonal symmetry-adapted blocks `Uᵢ mat Uⱼ†` (`i ≠ j`)
vanish after orbit reduction, without expanding them symbolically.

Each transformed entry is a *linear form* in the orbit-reduced moment
variables. A linear form that evaluates to ~0 at independent generic points
has ~0 coefficients, which is exactly what the `:full` check asserts
entry-by-entry. So: replace every canonical moment variable with a
deterministic pseudo-random weight, collapse `mat` to a numeric matrix, do
the congruence transform in BLAS, and bound the off-diagonal blocks. Repeats
with `ntrials` independent weight draws.

Throws the same kind of `ArgumentError` as the `:full` check on failure.
"""
function _verify_offblocks_randomized(
    mat::Matrix{P},
    row_bases::Vector{<:AbstractMatrix},
    reducer::Union{Nothing,_OrbitReducer{A,T}};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="constraint matrix",
    ntrials::Int=2,
) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C}}
    length(row_bases) <= 1 && return nothing

    n = size(mat, 1)
    size(mat, 2) == n || throw(ArgumentError("$context must be square for the off-block certificate."))

    # Sparse linear-form representation: mat = Σ_k A_k y_k with COO triplets.
    var_index = Dict{NormalMonomial{A,T},Int}()
    coo_rows = Int[]
    coo_cols = Int[]
    coo_vars = Int[]
    coo_coefs = ComplexF64[]

    for col in axes(mat, 2), row in axes(mat, 1)
        for (coef, mono) in mat[row, col].terms
            local key
            phase = 1
            if reducer === nothing
                key = mono
            else
                sym_mono = _symmetric_monomial(mono)
                haskey(reducer.representative, sym_mono) || throw(ArgumentError(
                    "Symmetry reducer is missing the monomial $sym_mono."
                ))
                phase = reducer.phase[sym_mono]
                phase == 0 && continue
                key = reducer.representative[sym_mono]
            end
            push!(coo_rows, row)
            push!(coo_cols, col)
            push!(coo_vars, get!(var_index, key, length(var_index) + 1))
            push!(coo_coefs, ComplexF64(coef) * phase)
        end
    end

    nvars = length(var_index)
    nvars == 0 && return nothing

    U = Matrix{ComplexF64}(reduce(vcat, (Matrix(b) for b in row_bases)))
    block_sizes = Int[size(b, 1) for b in row_bases]
    stops = cumsum(block_sizes)

    cmax = maximum(abs, coo_coefs)
    unorm2 = maximum(sum(abs2, U; dims=2))
    tol = atol * max(1.0, cmax) * max(1.0, unorm2) * nvars

    A_num = Matrix{ComplexF64}(undef, n, n)
    for trial in 1:ntrials
        seed = _splitmix64(0x4e435453536f5342 + UInt64(trial))  # "NCTSSoSB"
        weights = [_pseudorandom_weight(k, seed) for k in 1:nvars]

        Base.fill!(A_num, zero(ComplexF64))
        @inbounds for t in eachindex(coo_rows)
            A_num[coo_rows[t], coo_cols[t]] += coo_coefs[t] * weights[coo_vars[t]]
        end

        transformed = U * A_num * U'

        for bj in eachindex(block_sizes), bi in eachindex(block_sizes)
            bi == bj && continue
            (block_sizes[bi] == 0 || block_sizes[bj] == 0) && continue
            rows = (bi == 1 ? 1 : stops[bi-1] + 1):stops[bi]
            cols = (bj == 1 ? 1 : stops[bj-1] + 1):stops[bj]
            violation = maximum(abs, view(transformed, rows, cols))
            violation <= tol || throw(ArgumentError(
                "$context off-diagonal block ($bi, $bj) should vanish after symmetry reduction, " *
                "but the randomized certificate found magnitude $violation (tolerance $tol, trial $trial). " *
                "The supplied symmetry generators are likely not exact symmetries of this problem. " *
                "Rerun with `offblock_check = :full` in `SymmetrySpec` for an entry-level diagnosis."
            ))
        end
    end

    return nothing
end

function _reduce_transformed_blocks(
    mat::Matrix{P},
    row_bases::Vector{<:AbstractMatrix},
    reducer::Union{Nothing,_OrbitReducer};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="constraint matrix",
    offblock_check::Symbol=:full,
) where {P<:Polynomial}
    if offblock_check === :randomized
        _verify_offblocks_randomized(mat, row_bases, reducer; atol, context)
        return _diagonal_transformed_blocks(mat, row_bases, reducer; atol)
    elseif offblock_check === :off
        return _diagonal_transformed_blocks(mat, row_bases, reducer; atol)
    elseif offblock_check !== :full
        throw(ArgumentError("`offblock_check` must be `:full`, `:randomized`, or `:off`; got `$(repr(offblock_check))`."))
    end

    for i in 1:length(row_bases), j in (i+1):length(row_bases)
        off_block = _maybe_orbit_reduce_matrix(
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

        off_block_t = _maybe_orbit_reduce_matrix(
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

    return _diagonal_transformed_blocks(mat, row_bases, reducer; atol)
end

function _reduce_constraint_matrix_symmetric(
    mat::Matrix{P},
    basis::Vector{M},
    sw_group,
    reducer::Union{Nothing,_OrbitReducer{A,T}};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="constraint matrix",
    offblock_check::Symbol=:full,
) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C},M<:NormalMonomial{A,T}}
    if length(basis) == 1 || _basis_action_is_trivial(basis, sw_group)
        return [_maybe_orbit_reduce_matrix(mat, reducer; atol)]
    end

    blocks = _sw_decompose_half_basis(basis, sw_group)
    row_bases = [Matrix(block) for block in blocks]
    return _reduce_transformed_blocks(mat, row_bases, reducer; atol, context, offblock_check)
end

function _validate_fermionic_sector_spec(sector::FermionicSectorSpec)
    (sector.split_parity || sector.split_number || sector.split_spin || sector.split_irrep) ||
        throw(ArgumentError("`FermionicSectorSpec` must enable at least one sector split."))
    return nothing
end

function _neutral_fermionic_sector_label(sector::FermionicSectorSpec)
    return FermionicSectorLabel(
        sector.split_parity ? 0 : nothing,
        sector.split_number ? 0 : nothing,
        sector.split_spin ? 0 : nothing,
        sector.split_irrep ? _fermionic_irrep_identity(sector) : nothing,
    )
end

function _fermionic_mode_spin2(sector::FermionicSectorSpec, mode::Int)
    haskey(sector.mode_layout.spin2_of, mode) || throw(ArgumentError(
        "Fermionic spin-sector splitting needs `spin2_of[$mode]` in the supplied `FermionicModeLayout`."
    ))
    return sector.mode_layout.spin2_of[mode]
end

function _fermionic_sector_label(
    mono::NormalMonomial{FermionicAlgebra,T},
    sector::FermionicSectorSpec,
) where {T<:Integer}
    parity = sector.split_parity ? mod(length(mono.word), 2) : nothing
    number = nothing
    spin2 = sector.split_spin ? 0 : nothing
    irrep = sector.split_irrep ? _fermionic_irrep_identity(sector) : nothing

    if sector.split_number
        n_creators = count(<(zero(T)), mono.word)
        number = n_creators - (length(mono.word) - n_creators)
    end

    for idx in mono.word
        mode = Int(abs(idx))

        if sector.split_spin
            contribution = _fermionic_mode_spin2(sector, mode)
            spin2 += idx < 0 ? contribution : -contribution
        end

        if sector.split_irrep
            irrep = _fermionic_irrep_multiply(
                sector,
                irrep,
                idx < 0 ? _fermionic_mode_irrep(sector, mode) :
                    _fermionic_irrep_dual(sector, _fermionic_mode_irrep(sector, mode)),
            )
        end
    end

    return FermionicSectorLabel(parity, number, spin2, irrep)
end

function _fermionic_polynomial_sector_label(
    poly::Polynomial{FermionicAlgebra,T,C},
    sector::FermionicSectorSpec,
    context::AbstractString,
) where {T<:Integer,C<:Number}
    isempty(poly.terms) && return _neutral_fermionic_sector_label(sector)

    label = nothing
    for (_, mono) in poly.terms
        mono_label = _fermionic_sector_label(mono, sector)
        if isnothing(label)
            label = mono_label
        elseif mono_label != label
            throw(ArgumentError(
                "Fermionic sector splitting currently requires each $context to be homogeneous in the enabled sector labels. Mixed labels found in $poly."
            ))
        end
    end

    return something(label, _neutral_fermionic_sector_label(sector))
end

function _sector_partition_blocks(
    basis::AbstractVector{<:NormalMonomial{FermionicAlgebra,T}},
    sector::FermionicSectorSpec,
) where {T<:Integer}
    labels = FermionicSectorLabel[_fermionic_sector_label(mono, sector) for mono in basis]
    by_label = Dict{FermionicSectorLabel,Vector{Int}}()
    order = FermionicSectorLabel[]

    for (idx, label) in enumerate(labels)
        if !haskey(by_label, label)
            by_label[label] = Int[]
            push!(order, label)
        end
        push!(by_label[label], idx)
    end

    blocks = [_BasisPartitionBlock(by_label[label], label, :sector_split) for label in order]
    return blocks, labels
end

function _check_sectorwise_basis_closure(
    label::AbstractString,
    basis::AbstractVector{<:NormalMonomial{FermionicAlgebra,T}},
    group::AbstractVector,
    sector::FermionicSectorSpec,
) where {T<:Integer}
    sector_blocks, _ = _sector_partition_blocks(basis, sector)
    for block in sector_blocks
        _check_basis_closure("$label sector $(block.label)", basis[block.indices], group)
    end
    return nothing
end

function _append_sector_zero_constraints!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    mat::Matrix{MP},
    sector_blocks::Vector{<:_BasisPartitionBlock},
    reducer::Union{Nothing,_OrbitReducer};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="constraint matrix",
) where {MP<:Polynomial}
    for i in 1:length(sector_blocks), j in (i+1):length(sector_blocks)
        rows = sector_blocks[i].indices
        cols = sector_blocks[j].indices

        for row in rows, col in cols
            poly = reducer === nothing ? _chop_polynomial(mat[row, col]; atol) :
                _orbit_reduce_polynomial(_chop_polynomial(mat[row, col]; atol), reducer; atol)
            iszero(poly) || _append_constraint!(constraints, :Zero, reshape([poly], 1, 1), MP)
        end

        for row in cols, col in rows
            poly = reducer === nothing ? _chop_polynomial(mat[row, col]; atol) :
                _orbit_reduce_polynomial(_chop_polynomial(mat[row, col]; atol), reducer; atol)
            iszero(poly) || _append_constraint!(constraints, :Zero, reshape([poly], 1, 1), MP)
        end
    end

    return nothing
end

function _collect_reduced_basis(
    objective::Polynomial{A,T,C},
    constraints::Vector{Tuple{Symbol, Matrix{Polynomial{A,T,C}}}},
) where {A<:AlgebraType,T<:Integer,C<:Number}
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
    reducer::Union{Nothing,_OrbitReducer{A,T}},
    ::Type{MP};
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:AlgebraType,T<:Integer,CMP<:Number,CPP<:Number,M<:NormalMonomial{A,T},MP<:Polynomial{A,T,CMP},PP<:Polynomial{A,T,CPP}}
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
            poly = reducer === nothing ? _chop_polynomial(poly; atol) :
                _orbit_reduce_polynomial(_chop_polynomial(poly; atol), reducer; atol)
            iszero(poly) && continue
            _append_constraint!(constraints, :Zero, reshape([poly], 1, 1), MP)
        end
    end

    return nothing
end


"""
    moment_relax_symmetric(pop, corr_sparsity, cliques_term_sparsities, symmetry)

Construct a symmetry-reduced symbolic moment problem for dense ordinary polynomial
relaxations. The output contract is still the usual `MomentProblem`.
"""
function moment_relax_symmetric(
    pop::PolyOpt{A,T,P},
    corr_sparsity::CorrelativeSparsity{A,T,P,M,Nothing},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
    symmetry::SymmetrySpec;
    atol::Float64=_SYMMETRY_ATOL,
) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C},M<:NormalMonomial{A,T}}
    sector = symmetry.sector
    spin_adaptation = symmetry.spin_adaptation
    if !isnothing(sector)
        A === FermionicAlgebra || throw(ArgumentError(
            "Fermionic sector splitting is only supported for `FermionicAlgebra`; got `$(nameof(A))`."
        ))
        _validate_fermionic_sector_spec(sector)
    end
    if !isnothing(spin_adaptation)
        A === FermionicAlgebra || throw(ArgumentError(
            "Fermionic spin adaptation is only supported for `FermionicAlgebra`; got `$(nameof(A))`."
        ))
    end
    if !isempty(symmetry.clifford_generators)
        A === PauliAlgebra || throw(ArgumentError(
            "Clifford symmetry is only supported for `PauliAlgebra`; got `$(nameof(A))`."
        ))
    end

    moment_matrix_basis = _moment_matrix_basis(cliques_term_sparsities)
    total_basis, moment_eq_row_bases, moment_eq_row_basis_degrees =
        _polynomial_total_basis(pop, corr_sparsity, cliques_term_sparsities)
    _validate_polynomial_relaxation_support(pop, total_basis; source="Constructed relaxation basis")

    active_modes = if A === FermionicAlgebra && (!isnothing(sector) || !isnothing(spin_adaptation) || !isempty(symmetry.fermionic_generators))
        _mode_symmetry_domain(pop, corr_sparsity, cliques_term_sparsities)
    else
        Int[]
    end
    !isnothing(sector) && _validate_active_fermionic_sector_metadata(sector, active_modes)
    !isnothing(spin_adaptation) && _validate_fermionic_spin_adaptation_spec(
        symmetry,
        active_modes,
    )

    group = nothing
    sw_group = nothing
    reducer = nothing
    group_order = 1

    if !isempty(symmetry.generators)
        spec = _convert_symmetry_spec(T, symmetry)
        domain = _symmetry_domain(pop, corr_sparsity, cliques_term_sparsities)
        group = _enumerate_symmetry_group(spec, domain)
        sw_group = _sw_signed_permutation_group(spec, group, domain)
    elseif !isempty(symmetry.fermionic_generators)
        domain = _mode_symmetry_domain(pop, corr_sparsity, cliques_term_sparsities)
        group = _enumerate_symmetry_group(symmetry, domain)
        sw_group = _sw_fermionic_mode_permutation_group(symmetry, group, domain)
    elseif !isempty(symmetry.clifford_generators)
        support_domain = _symmetry_domain(pop, corr_sparsity, cliques_term_sparsities)
        nqubits = max(
            _clifford_nqubits_from_domain(support_domain),
            _clifford_max_generator_nqubits(symmetry.clifford_generators),
        )
        sw_group = CliffordSymmetryGroup(symmetry.clifford_generators; nqubits, integer_type=T, domain=support_domain)
        group = CliffordSymmetry{T}[_clifford_group_value(element) for element in sw_group]
    end

    if !isnothing(group)
        symmetry.check_invariance && _check_symmetry_invariance(pop, group)

        for (clique_idx, basis) in enumerate(corr_sparsity.clq_mom_mtx_bases)
            if !isnothing(sector) && A === FermionicAlgebra
                _check_sectorwise_basis_closure("moment basis of clique $clique_idx", basis, group, sector)
            else
                _check_basis_closure("moment basis of clique $clique_idx", basis, group)
            end
        end

        for (clique_idx, term_sparsities) in enumerate(cliques_term_sparsities)
            for (poly_idx, term_sparsity) in enumerate(term_sparsities)
                for (block_idx, basis) in enumerate(term_sparsity.block_bases)
                    if !isnothing(sector) && A === FermionicAlgebra
                        _check_sectorwise_basis_closure(
                            "constraint block basis $(poly_idx).$(block_idx) of clique $clique_idx",
                            basis,
                            group,
                            sector,
                        )
                    else
                        _check_basis_closure(
                            "constraint block basis $(poly_idx).$(block_idx) of clique $clique_idx",
                            basis,
                            group,
                        )
                    end
                end
            end
        end

        reducer = _build_orbit_reducer(total_basis, group)
        group_order = length(group)
    end

    psd_cone = _is_complex_problem(A) ? :HPSD : :PSD
    MP_C = _moment_problem_coeff_type(A, C)
    MP_P = Polynomial{A,T,MP_C}

    objective_mp = convert(MP_P, pop.objective)
    objective_mp = reducer === nothing ? _chop_polynomial(objective_mp; atol) :
        _orbit_reduce_polynomial(objective_mp, reducer; atol)

    if !isnothing(sector)
        neutral_label = _neutral_fermionic_sector_label(sector)
        objective_label = _fermionic_polynomial_sector_label(objective_mp, sector, "objective")
        objective_label == neutral_label || throw(ArgumentError(
            "Fermionic sector splitting currently supports only sector-neutral objectives. Got label $objective_label for the objective."
        ))
    end

    constraints = Tuple{Symbol, Matrix{MP_P}}[]
    reduced_moment_terms = NormalMonomial{A,T}[]
    reduced_moment_block_sizes = Int[]
    moment_block_labels = Any[]
    moment_block_provenance = Symbol[]

    for (clique_idx, (term_sparsities, cons_idx)) in enumerate(zip(cliques_term_sparsities, corr_sparsity.clq_cons))
        polys = [one(pop.objective); corr_sparsity.cons[cons_idx]...]

        for (poly_idx, (term_sparsity, poly)) in enumerate(zip(term_sparsities, polys))
            for (block_idx, basis) in enumerate(term_sparsity.block_bases)
                cone = poly in pop.eq_constraints ? :Zero : psd_cone
                _, raw_mat = _build_constraint_matrix(poly, basis, cone)
                mat = _convert_polynomial_matrix(MP_P, raw_mat)

                diagonal_blocks = Matrix{MP_P}[]
                diagonal_labels = Any[]
                diagonal_provenance = Symbol[]
                context = "clique $clique_idx block $block_idx for polynomial $poly_idx"

                if !isnothing(sector)
                    neutral_label = _neutral_fermionic_sector_label(sector)
                    poly_label = _fermionic_polynomial_sector_label(convert(MP_P, poly), sector, context)
                    poly_label == neutral_label || throw(ArgumentError(
                        "Fermionic sector splitting currently supports only sector-neutral multipliers. Got label $poly_label for $context."
                    ))

                    sector_blocks, _ = _sector_partition_blocks(basis, sector)
                    _append_sector_zero_constraints!(constraints, mat, sector_blocks, reducer; atol, context)

                    for sector_block in sector_blocks
                        sector_basis = basis[sector_block.indices]
                        sector_mat = mat[sector_block.indices, sector_block.indices]
                        transform_blocks = _fermionic_sector_transform_blocks(
                            sector_basis,
                            sector_block.label,
                            sw_group,
                            spin_adaptation,
                            active_modes;
                            atol,
                            context="$context sector $(sector_block.label)",
                        )
                        row_bases = [block.row_basis for block in transform_blocks]
                        if isnothing(spin_adaptation)
                            sector_diag_blocks = _reduce_transformed_blocks(
                                sector_mat,
                                row_bases,
                                reducer;
                                atol,
                                context="$context sector $(sector_block.label)",
                                offblock_check=symmetry.offblock_check,
                            )
                        else
                            _append_transformed_offblock_zero_constraints!(
                                constraints,
                                sector_mat,
                                row_bases,
                                reducer;
                                atol,
                                context="$context sector $(sector_block.label)",
                            )
                            sector_diag_blocks = _diagonal_transformed_blocks(
                                sector_mat,
                                row_bases,
                                reducer;
                                atol,
                            )
                        end

                        append!(diagonal_blocks, sector_diag_blocks)
                        append!(diagonal_labels, [block.label for block in transform_blocks])
                        append!(diagonal_provenance, [block.provenance for block in transform_blocks])
                    end
                else
                    diagonal_blocks = isnothing(sw_group) ?
                        [_maybe_orbit_reduce_matrix(mat, reducer; atol)] :
                        _reduce_constraint_matrix_symmetric(
                            mat,
                            basis,
                            sw_group,
                            reducer;
                            atol,
                            context=context,
                            offblock_check=symmetry.offblock_check,
                        )
                    diagonal_labels = fill(nothing, length(diagonal_blocks))
                    diagonal_provenance = fill(isnothing(sw_group) ? :identity : :wedderburn, length(diagonal_blocks))
                end

                for (diag_block, label, provenance) in zip(diagonal_blocks, diagonal_labels, diagonal_provenance)
                    _append_constraint!(constraints, cone, diag_block, MP_P)
                    if poly_idx == 1
                        for poly_entry in diag_block
                            append!(reduced_moment_terms, monomials(poly_entry))
                        end
                        push!(reduced_moment_block_sizes, size(diag_block, 1))
                        push!(moment_block_labels, label)
                        push!(moment_block_provenance, provenance)
                    end
                end
            end
        end
    end

    for global_con in corr_sparsity.global_cons
        poly = corr_sparsity.cons[global_con]
        cone = poly in pop.eq_constraints ? :Zero : psd_cone
        reduced_poly = convert(MP_P, poly)
        reduced_poly = reducer === nothing ? _chop_polynomial(reduced_poly; atol) :
            _orbit_reduce_polynomial(reduced_poly, reducer; atol)

        if !isnothing(sector)
            neutral_label = _neutral_fermionic_sector_label(sector)
            poly_label = _fermionic_polynomial_sector_label(reduced_poly, sector, "global constraint $global_con")
            poly_label == neutral_label || throw(ArgumentError(
                "Fermionic sector splitting currently supports only sector-neutral global constraints. Got label $poly_label for global constraint $global_con."
            ))
        end

        _append_constraint!(constraints, cone, reshape([reduced_poly], 1, 1), MP_P)
    end

    if !isnothing(sector)
        neutral_label = _neutral_fermionic_sector_label(sector)
        for (i, g) in pairs(pop.moment_eq_constraints)
            label = _fermionic_polynomial_sector_label(convert(MP_P, g), sector, "moment equality constraint $i")
            label == neutral_label || throw(ArgumentError(
                "Fermionic sector splitting currently supports only sector-neutral moment equality constraints. Got label $label for moment equality constraint $i."
            ))
        end
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
    _add_parity_constraints!(mp)

    report = SymmetryReport(
        group_order,
        count(mono -> !isone(mono), reduced_moment_basis),
        reduced_moment_block_sizes,
        length(only(corr_sparsity.clq_mom_mtx_bases)),
        length(moment_matrix_basis),
        moment_block_labels,
        moment_block_provenance,
    )

    return mp, report
end
