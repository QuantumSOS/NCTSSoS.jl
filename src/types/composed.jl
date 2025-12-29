"""
    ComposedMonomial{Ts<:Tuple} <: AbstractMonomial

Represents a product of monomials from DIFFERENT algebra types.

Each component is a `Monomial{A,T}` with its own algebra type in the type parameter.
This enables tensor products of operators from different algebraic structures.

# Fields
- `components::Ts`: Tuple of monomials, e.g., `(Monomial{PauliAlgebra}, Monomial{FermionicAlgebra})`

# Type Parameters
- `Ts<:Tuple`: Tuple type of the component monomials

# Design
Algebra types are compile-time information (in each monomial's type parameter),
so ComposedMonomial has zero runtime overhead for algebra type storage.

# Examples
```jldoctest
julia> m_pauli = Monomial{PauliAlgebra}([1, 2, 3]);

julia> m_fermi = Monomial{FermionicAlgebra}(Int32[-1, 2]);

julia> cm = ComposedMonomial((m_pauli, m_fermi));

julia> length(cm)
2

julia> degree(cm)
5

julia> cm[1] === m_pauli
true
```

Simplification dispatches to each component's algebra and returns `Vector{Term}`:
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> m_pauli = Monomial{PauliAlgebra}(UInt16[1, 1, 2]);

julia> m_unip = Monomial{UnipotentAlgebra}(UInt16[encode_index(UInt16, 1, 1), encode_index(UInt16, 1, 1), encode_index(UInt16, 2, 1)]);

julia> cm = ComposedMonomial((m_pauli, m_unip));

julia> terms = simplify(cm);

julia> terms[1].coefficient
1.0 + 0.0im

julia> terms[1].monomial[1].word  # Pauli: [1,1,2] -> [2] (sigma_x squared = I)
1-element Vector{UInt16}:
 0x0002

julia> terms[1].monomial[2].word == [encode_index(UInt16, 2, 1)]  # Unipotent: [1,1,2] -> [2]
true
```
"""
struct ComposedMonomial{Ts<:Tuple} <: AbstractMonomial
    components::Ts
end

# Note: Default constructor ComposedMonomial(components::Tuple) is provided by Julia
# since the struct only has one field.

"""
    Base.:(==)(cm1::ComposedMonomial, cm2::ComposedMonomial) -> Bool

Equality check via component-wise comparison.
"""
function Base.:(==)(cm1::ComposedMonomial{T1}, cm2::ComposedMonomial{T2}) where {T1,T2}
    T1 !== T2 && return false
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
    Base.isless(cm1::ComposedMonomial{Ts}, cm2::ComposedMonomial{Ts}) where {Ts<:Tuple}

Compare two ComposedMonomials using degree-first ordering, then component-wise lexicographic.

This ordering enables sorting ComposedMonomials for polynomial operations.

# Algorithm
1. Compare total degrees (degree-first ordering)
2. If degrees equal, compare components lexicographically using their isless

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}(UInt16[1]);

julia> m2 = Monomial{PauliAlgebra}(UInt16[1, 2]);

julia> m3 = Monomial{FermionicAlgebra}(Int32[1]);

julia> cm1 = ComposedMonomial((m1, m3));  # degree 2

julia> cm2 = ComposedMonomial((m2, m3));  # degree 3

julia> isless(cm1, cm2)  # degree 2 < degree 3
true
```
"""
function Base.isless(cm1::ComposedMonomial{Ts}, cm2::ComposedMonomial{Ts}) where {Ts<:Tuple}
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
    Base.getindex(cm::ComposedMonomial, i::Int) -> Monomial

Access the i-th component monomial.
"""
Base.getindex(cm::ComposedMonomial, i::Int) = cm.components[i]

"""
    degree(cm::ComposedMonomial) -> Int

Total degree across all components.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}(UInt16[1, 2, 3]);

julia> m2 = Monomial{FermionicAlgebra}(Int32[1, 2]);

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
    Base.one(::Type{ComposedMonomial{Ts}}) where {Ts<:Tuple}

Create the identity ComposedMonomial (all components are identity monomials).
"""
function Base.one(::Type{ComposedMonomial{Ts}}) where {Ts<:Tuple}
    identity_components = ntuple(i -> one(fieldtype(Ts, i)), fieldcount(Ts))
    return ComposedMonomial(identity_components)
end

"""
    Base.one(cm::ComposedMonomial{Ts}) where {Ts}

Create the identity ComposedMonomial for the same type as `cm`.
"""
function Base.one(cm::ComposedMonomial{Ts}) where {Ts}
    return one(ComposedMonomial{Ts})
end

"""
    Base.isone(t::Term{ComposedMonomial{Ts},C}) where {Ts,C} -> Bool

Check if a composed term is the identity (coefficient 1, all empty monomials).
"""
function Base.isone(t::Term{ComposedMonomial{Ts},C}) where {Ts,C}
    isone(t.coefficient) || return false
    for mono in t.monomial.components
        isempty(mono.word) || return false
    end
    return true
end

"""
    simplify(cm::ComposedMonomial) -> Vector{Term{<:ComposedMonomial}}

Simplify each component according to its algebra type.

Always returns `Vector{Term{ComposedMonomial}}` for consistent API. Each term
contains a ComposedMonomial with the simplified component monomials.

# Algorithm
1. Simplify each component using its algebra's simplify function
2. Convert all results to (coefficient, monomial) pairs
3. Compute Cartesian product across all components
4. Filter zero terms

# Examples
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> m_pauli = Monomial{PauliAlgebra}(UInt16[1, 1]);

julia> u2 = encode_index(UInt16, 2, 1);

julia> m_unip = Monomial{UnipotentAlgebra}([u2, u2]);

julia> cm = ComposedMonomial((m_pauli, m_unip));

julia> result = simplify(cm);

julia> result isa Vector{<:Term}
true

julia> result[1].coefficient
1.0 + 0.0im

julia> isempty(result[1].monomial[1].word)  # Pauli [1,1] -> [] (sigma_x squared = I)
true

julia> isempty(result[1].monomial[2].word)  # Unipotent [2,2] -> []
true
```
"""
function simplify(cm::ComposedMonomial)
    simplified_components = map(simplify, cm.components)
    return _expand_simplified_components(simplified_components)
end

# Helper to expand Cartesian product of simplified components
function _expand_simplified_components(simplified::Tuple)
    # Infer coefficient type from tuple element types at compile time
    CoefType = _infer_coef_type_from_types(simplified)

    # Compute Cartesian product using iteration protocol
    result = Term[]

    _cartesian_product_iter!(result, simplified, 1, (), one(CoefType))

    # Filter zeros and return
    filter!(!iszero, result)

    if isempty(result)
        # Return zero term with identity monomials
        one_components = map(_get_identity_monomial, simplified)
        return [Term(zero(CoefType), ComposedMonomial(one_components))]
    end

    return result
end

# Get identity monomial from a simplified component (avoids brittle type reflection)
function _get_identity_monomial(x::Monomial)
    return one(x)
end

function _get_identity_monomial(x::Term)
    return one(x.monomial)
end

function _get_identity_monomial(x::Polynomial)
    if isempty(x.terms)
        # Empty polynomial - create identity from the Polynomial's type
        # Use the concrete monomial type from the polynomial
        return one(eltype(x.terms).parameters[1])
    else
        return one(x.terms[1].monomial)
    end
end

function _get_identity_monomial(x::Vector{<:Term})
    return one(x[1].monomial)
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
    result::Vector{Term},
    components::Tuple,
    idx::Int,
    current_monomials::Tuple,
    current_coef::Number,
)
    if idx > length(components)
        # Base case: all components processed
        cm = ComposedMonomial(current_monomials)
        push!(result, Term(current_coef, cm))
        return nothing
    end

    # Iterate over simplified component using new iteration protocol
    for (coef, mono) in components[idx]
        new_coef = current_coef * coef
        new_monomials = (current_monomials..., mono)
        _cartesian_product_iter!(result, components, idx + 1, new_monomials, new_coef)
    end
end

# Show method for ComposedMonomial
function Base.show(io::IO, cm::ComposedMonomial)
    print(io, "ComposedMonomial(")
    for (i, mono) in enumerate(cm.components)
        i > 1 && print(io, ", ")
        print(io, mono.word)
    end
    return print(io, ")")
end

# Show method for Term with ComposedMonomial
function Base.show(io::IO, t::Term{<:ComposedMonomial,C}) where {C}
    if iszero(t)
        print(io, "0")
    elseif isone(t)
        print(io, "1")
    elseif t.coefficient == one(C)
        print(io, t.monomial)
    else
        print(io, t.coefficient, " * ", t.monomial)
    end
end
