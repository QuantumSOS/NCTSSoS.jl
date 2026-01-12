"""
    ComposedMonomial{As<:Tuple,Ts<:Tuple} <: AbstractTensorMonomial{As}

Represents a product of monomials from DIFFERENT algebra types.

Each component is a `NormalMonomial{A,T}` with its own algebra type in the type parameter.
This enables tensor products of operators from different algebraic structures.

# Fields
- `components::Ts`: Tuple of monomials, e.g., `(NormalMonomial{PauliAlgebra}, NormalMonomial{FermionicAlgebra})`

# Type Parameters
- `As<:Tuple`: Ordered algebra signature, e.g. `Tuple{PauliAlgebra,FermionicAlgebra}`
- `Ts<:Tuple`: Tuple type of the component monomials

# Design
Algebra types are compile-time information (in each monomial's type parameter),
so ComposedMonomial has zero runtime overhead for algebra type storage.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m_pauli = NormalMonomial{PauliAlgebra}(UInt16[1, 4, 7]);  # σx on sites 1,2,3

julia> m_fermi = NormalMonomial{FermionicAlgebra}(Int32[-1, 2]);

julia> cm = ComposedMonomial((m_pauli, m_fermi));

julia> length(cm)
2

julia> degree(cm)
5

julia> cm[1] === m_pauli
true
```

Simplification dispatches to each component's algebra and returns a vector of `(coefficient, ComposedMonomial)` pairs:
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> m_pauli = NormalMonomial{PauliAlgebra}(UInt16[1, 4]);  # σx on sites 1,2

julia> m_unip = NormalMonomial{UnipotentAlgebra}(UInt16[encode_index(UInt16, 1, 1), encode_index(UInt16, 1, 1), encode_index(UInt16, 2, 1)]);

julia> cm = ComposedMonomial((m_pauli, m_unip));

julia> terms = simplify(cm);

julia> terms[1][1]
1.0 + 0.0im

julia> terms[1][2][1].word == UInt16[1, 4]  # Pauli component stays canonical
true

julia> terms[1][2][2].word == [encode_index(UInt16, 2, 1)]  # Unipotent: [1,1,2] -> [2]
true
```
"""
@inline function _unwrap_unionall(T)
    while T isa UnionAll
        T = T.body
    end
    return T
end

@generated function _composed_signature(::Type{Ts}) where {Ts<:Tuple}
    comps = Ts.parameters
    isempty(comps) && return :(Tuple{})

    algs = Any[]
    for C in comps
        sup = _unwrap_unionall(supertype(C))
        if !(sup <: AbstractMonomial)
            error("ComposedMonomial components must be subtypes of AbstractMonomial{A,T}, got $(C)")
        end
        push!(algs, sup.parameters[1])
    end

    return :(Tuple{$(algs...)})
end

# -----------------------------------------------------------------------------
# Type definition + constructors
# -----------------------------------------------------------------------------

struct ComposedMonomial{As<:Tuple,Ts<:Tuple} <: AbstractTensorMonomial{As}
    components::Ts

    function ComposedMonomial{As,Ts}(components::Ts) where {As<:Tuple,Ts<:Tuple}
        if As !== _composed_signature(Ts)
            throw(ArgumentError(
                "ComposedMonomial signature mismatch: As=$As but components imply $(_composed_signature(Ts))"
            ))
        end
        return new{As,Ts}(components)
    end
end

function ComposedMonomial(components::Ts) where {Ts<:Tuple}
    return ComposedMonomial{_composed_signature(Ts),Ts}(components)
end

"""
    Base.:(==)(cm1::ComposedMonomial, cm2::ComposedMonomial) -> Bool

Equality check via component-wise comparison.
"""
function Base.:(==)(
    cm1::ComposedMonomial{As1,Ts1},
    cm2::ComposedMonomial{As2,Ts2},
) where {As1,Ts1,As2,Ts2}
    (As1 !== As2 || Ts1 !== Ts2) && return false
    return cm1.components == cm2.components
end

"""
    Base.hash(cm::ComposedMonomial, h::UInt) -> UInt

Hash function computed from components.
"""
function Base.hash(cm::ComposedMonomial, h::UInt)
    h = hash(length(cm.components), h)
    for mono in cm.components
        h = hash(mono, h)
    end
    return h
end

"""
    Base.isless(cm1::ComposedMonomial{As,Ts}, cm2::ComposedMonomial{As,Ts}) where {As<:Tuple,Ts<:Tuple}

Compare two ComposedMonomials using degree-first ordering, then component-wise lexicographic.

This ordering enables sorting ComposedMonomials for polynomial operations.

# Algorithm
1. Compare total degrees (degree-first ordering)
2. If degrees equal, compare components lexicographically using their isless

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}(UInt16[1]);

julia> m2 = NormalMonomial{PauliAlgebra}(UInt16[1, 2]);

julia> m3 = NormalMonomial{FermionicAlgebra}(Int32[1]);

julia> cm1 = ComposedMonomial((m1, m3));  # degree 2

julia> cm2 = ComposedMonomial((m2, m3));  # degree 3

julia> isless(cm1, cm2)  # degree 2 < degree 3
true
```
"""
function Base.isless(cm1::ComposedMonomial{As,Ts}, cm2::ComposedMonomial{As,Ts}) where {As<:Tuple,Ts<:Tuple}
    # Degree-first ordering
    d1, d2 = degree(cm1), degree(cm2)
    d1 != d2 && return d1 < d2

    # Component-wise lexicographic ordering
    for (c1, c2) in zip(cm1.components, cm2.components)
        c1 != c2 && return isless(c1, c2)
    end
    return false  # Equal
end

"""
    Base.length(cm::ComposedMonomial) -> Int

Number of component monomials.
"""
Base.length(cm::ComposedMonomial) = length(cm.components)

"""
    Base.getindex(cm::ComposedMonomial, i::Int) -> NormalMonomial

Access the i-th component monomial.
"""
Base.getindex(cm::ComposedMonomial, i::Int) = cm.components[i]

"""
    degree(cm::ComposedMonomial) -> Int

Total degree across all components.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}(UInt16[1, 2, 3]);

julia> m2 = NormalMonomial{FermionicAlgebra}(Int32[1, 2]);

julia> cm = ComposedMonomial((m1, m2));

julia> degree(cm)
5
```
"""
degree(cm::ComposedMonomial) = mapreduce(degree, +, cm.components)

"""
    Base.isone(cm::ComposedMonomial) -> Bool

Check if a ComposedMonomial is the identity (all components are identity monomials).
"""
function Base.isone(cm::ComposedMonomial)
    return all(isone, cm.components)
end

"""
    Base.one(::Type{ComposedMonomial{As,Ts}}) where {As<:Tuple,Ts<:Tuple}

Create the identity ComposedMonomial (all components are identity monomials).
"""
function Base.one(::Type{ComposedMonomial{As,Ts}}) where {As<:Tuple,Ts<:Tuple}
    identity_components = ntuple(i -> one(fieldtype(Ts, i)), fieldcount(Ts))
    return ComposedMonomial(identity_components)
end

"""
    Base.one(cm::ComposedMonomial{As,Ts}) where {As,Ts}

Create the identity ComposedMonomial for the same type as `cm`.
"""
function Base.one(cm::ComposedMonomial{As,Ts}) where {As,Ts}
    return one(ComposedMonomial{As,Ts})
end

"""
    Base.isone(t::Tuple{C,<:ComposedMonomial}) where {C} -> Bool

Check if a composed term is the identity (coefficient 1, all components are identity monomials).
"""
function Base.isone(t::Tuple{C,<:ComposedMonomial}) where {C<:Number}
    coef, cm = t
    return isone(coef) && isone(cm)
end

"""
    simplify(cm::ComposedMonomial) -> Vector{Tuple{<:Number,<:ComposedMonomial}}

Simplify each component according to its algebra type.

Always returns `Vector{Tuple{coefficient,ComposedMonomial}}` for consistent API.
Each entry contains a ComposedMonomial with the simplified component monomials.

# Algorithm
1. Simplify each component using its algebra's simplify function
2. Convert all results to (coefficient, monomial) pairs
3. Compute Cartesian product across all components
4. Filter zero terms

# Examples
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> m_pauli = NormalMonomial{PauliAlgebra}(UInt16[]);  # identity

julia> u2 = encode_index(UInt16, 2, 1);

julia> m_unip = NormalMonomial{UnipotentAlgebra}(UInt16[u2, u2]);

julia> cm = ComposedMonomial((m_pauli, m_unip));

julia> result = simplify(cm);

julia> result isa Vector{<:Tuple}
true

julia> result[1][1]
1.0 + 0.0im

julia> isempty(result[1][2][1].word)  # Pauli identity stays identity
true

julia> isempty(result[1][2][2].word)  # Unipotent [2,2] -> []
true
```
"""
function simplify(cm::ComposedMonomial{As,Ts}) where {As<:Tuple,Ts<:Tuple}
    simplified_components = map(simplify, cm.components)
    return _expand_simplified_components(ComposedMonomial{As,Ts}, simplified_components)
end

# Helper to expand Cartesian product of simplified components
function _expand_simplified_components(::Type{CM}, simplified::Tuple) where {CM<:ComposedMonomial}
    # Infer coefficient type from tuple element types at compile time
    CoefType = _infer_coef_type_from_types(simplified)

    # Compute Cartesian product using iteration protocol
    result = Tuple{CoefType,CM}[]

    _cartesian_product_iter!(result, simplified, 1, (), one(CoefType))

    # Filter zeros and return
    filter!(t -> !iszero(t[1]), result)

    if isempty(result)
        # Return zero term with identity monomials
        one_components = map(_get_identity_monomial, simplified)
        return Tuple{CoefType,CM}[(zero(CoefType), ComposedMonomial(one_components))]
    end

    return result
end

# Get identity monomial from a simplified component (for zero-result fallback)
function _get_identity_monomial(::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    return one(NormalMonomial{A,T})
end

# Compile-time coefficient type inference using coeff_type
@generated function _infer_coef_type_from_types(::T) where {T<:Tuple}
    types = T.parameters
    result_type = Float64
    for t in types
        ct = coeff_type(t)
        result_type = promote_type(result_type, ct)
    end
    return :($result_type)
end

# New helper using iteration protocol
function _cartesian_product_iter!(
    result::Vector{Tuple{CoefType,CM}},
    components::Tuple,
    idx::Int,
    current_monomials::Tuple,
    current_coef::CoefType,
) where {CoefType<:Number,CM<:ComposedMonomial}
    if idx > length(components)
        # Base case: all components processed
        cm = ComposedMonomial(current_monomials)
        push!(result, (current_coef, cm))
        return nothing
    end

    # Iterate over simplified component using new iteration protocol
    for (coef_internal, mono) in components[idx]
        new_coef = current_coef * _coeff_to_number(mono, coef_internal)
        new_monomials = (current_monomials..., mono)
        _cartesian_product_iter!(result, components, idx + 1, new_monomials, new_coef)
    end
end

# Show method for ComposedMonomial
function Base.show(io::IO, cm::ComposedMonomial)
    print(io, "ComposedMonomial(")
    for (i, mono) in enumerate(cm.components)
        i > 1 && print(io, ", ")
        show(io, mono)
    end
    return print(io, ")")
end

# Minimal helpers for `(coefficient, ComposedMonomial)` pairs (legacy `Term` replacement)
Base.iszero(t::Tuple{C,<:ComposedMonomial}) where {C<:Number} = iszero(t[1])

function Base.show(io::IO, t::Tuple{C,<:ComposedMonomial}) where {C<:Number}
    coef, mono = t
    if iszero(coef)
        print(io, "0")
    elseif isone(coef) && isone(mono)
        print(io, "1")
    elseif coef == one(C)
        print(io, mono)
    else
        print(io, coef, " * ", mono)
    end
end
