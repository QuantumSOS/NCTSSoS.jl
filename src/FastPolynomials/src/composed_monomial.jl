"""
    ComposedMonomial{Ts<:Tuple} <: AbstractMonomial

Represents a product of monomials from DIFFERENT algebra types.

Each component is a `Monomial{A,T}` with its own algebra type in the type parameter.
This enables tensor products of operators from different algebraic structures.

# Fields
- `components::Ts`: Tuple of monomials, e.g., `(Monomial{PauliAlgebra}, Monomial{FermionicAlgebra})`
- `hash::UInt64`: Precomputed hash for fast equality checks

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
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

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
    hash::UInt64
end

"""
    ComposedMonomial(components::Tuple)

Construct a ComposedMonomial from a tuple of monomials.

# Examples
```jldoctest
julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{FermionicAlgebra}(Int32[-1, 2]);

julia> cm = ComposedMonomial((m1, m2));

julia> length(cm)
2
```
"""
function ComposedMonomial(components::Ts) where {Ts<:Tuple}
    h = _compute_composed_hash(components)
    return ComposedMonomial{Ts}(components, h)
end

"""
    _compute_composed_hash(components::Tuple) -> UInt64

Compute hash from component monomials.
"""
function _compute_composed_hash(components::Tuple)
    h = hash(length(components))
    for mono in components
        h = hash(mono, h)
    end
    return h
end

"""
    Base.:(==)(cm1::ComposedMonomial, cm2::ComposedMonomial) -> Bool

Equality check using precomputed hash for fast inequality detection.
"""
function Base.:(==)(cm1::ComposedMonomial{T1}, cm2::ComposedMonomial{T2}) where {T1,T2}
    # Different tuple types means different structure
    T1 !== T2 && return false
    # Fast path: hash mismatch
    cm1.hash == cm2.hash || return false
    # Slow path: component-wise comparison
    return cm1.components == cm2.components
end

"""
    Base.hash(cm::ComposedMonomial, h::UInt) -> UInt

Hash function using precomputed hash.
"""
Base.hash(cm::ComposedMonomial, h::UInt) = hash(cm.hash, h)

"""
    Base.isless(cm1::ComposedMonomial{Ts}, cm2::ComposedMonomial{Ts}) where {Ts<:Tuple}

Compare two ComposedMonomials using degree-first ordering, then component-wise lexicographic.

This ordering enables sorting ComposedMonomials for polynomial operations.

# Algorithm
1. Compare total degrees (degree-first ordering)
2. If degrees equal, compare components lexicographically using their isless

# Examples
```jldoctest
julia> using FastPolynomials

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
julia> using FastPolynomials

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
    simplify(cm::ComposedMonomial) -> Union{ComposedMonomial, Term, Polynomial}

Simplify each component according to its algebra type.

The return type depends on the component simplification results:
- If any component returns `Polynomial` (Bosonic, Fermionic) → returns `Polynomial`
- Else if any component returns `Term` (Pauli) → returns `Term`
- Else (all components return `Monomial`: NC, Projector, Unipotent) → returns `ComposedMonomial`

This type hierarchy ensures that more complex return types propagate upward.

# Algorithm
1. Simplify each component using its algebra's simplify function
2. Determine return type based on component results
3. For multi-term algebras (Bosonic, Fermionic), compute Cartesian product

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> m_pauli = Monomial{PauliAlgebra}(UInt16[1, 1]);

julia> u2 = encode_index(UInt16, 2, 1);

julia> m_unip = Monomial{UnipotentAlgebra}([u2, u2]);

julia> cm = ComposedMonomial((m_pauli, m_unip));

julia> result = simplify(cm);  # Pauli returns Term, so result is Term

julia> result isa Term
true

julia> result.coefficient
1.0 + 0.0im

julia> isempty(result.monomial[1].word)  # Pauli [1,1] -> [] (sigma_x squared = I)
true

julia> isempty(result.monomial[2].word)  # Unipotent [2,2] -> []
true
```
"""
function simplify(cm::ComposedMonomial)
    # Simplify each component - collect results
    simplified_components = map(simplify, cm.components)

    # Determine return type based on component types:
    # - If any component is Polynomial -> return Polynomial
    # - Else if any component is Term -> return Term
    # - Else (all Monomial) -> return Monomial
    has_polynomial = any(x -> x isa Polynomial, simplified_components)
    has_term = any(x -> x isa Term, simplified_components)

    if has_polynomial
        # Return Vector{Term} for ComposedMonomial (Polynomial type doesn't support ComposedMonomial)
        term_vec = _expand_simplified_components(simplified_components)
        return term_vec
    elseif has_term
        # Return Term (single term from Cartesian product)
        term_vec = _expand_simplified_components(simplified_components)
        # Should always be single term when no Polynomial components
        return length(term_vec) == 1 ? term_vec[1] : term_vec
    else
        # All Monomials - return composed Monomial
        # Combine coefficients (should all be 1.0 for Monomials)
        coef = 1.0
        return ComposedMonomial(simplified_components)
    end
end

# Helper to expand Cartesian product of simplified components
function _expand_simplified_components(simplified::Tuple)
    # Convert each component to a vector of (coef, monomial) pairs
    component_terms = map(_to_term_vector, simplified)

    # Infer coefficient type from components (promotes Float64 + ComplexF64 -> ComplexF64)
    CoefType = _infer_coef_type(component_terms)

    # Compute Cartesian product
    result = Term[]

    _cartesian_product!(result, component_terms, 1, (), one(CoefType))

    # Filter zeros and return
    filter!(!iszero, result)

    if isempty(result)
        # Return zero term with identity monomials
        one_components = map(simplified) do x
            if x isa Polynomial
                isempty(x.terms) ? one(Monomial{typeof(x).parameters[1],typeof(x).parameters[2]}) : one(x.terms[1].monomial)
            elseif x isa Term
                one(x.monomial)
            elseif x isa Monomial
                one(x)
            else
                # Vector{Term}
                one(x[1].monomial)
            end
        end
        return [Term(zero(CoefType), ComposedMonomial(one_components))]
    end

    return result
end

# Infer the coefficient type from component terms
function _infer_coef_type(component_terms::Tuple)
    T = Float64
    for terms in component_terms
        for (coef, _) in terms
            T = promote_type(T, typeof(coef))
        end
    end
    return T
end

# Convert a Monomial to vector format (coefficient 1.0)
function _to_term_vector(m::Monomial)
    return [(1.0, m)]
end

# Convert a single Term to vector format
function _to_term_vector(t::Term)
    return [(t.coefficient, t.monomial)]
end

# Convert Vector{Term} to vector format (already in right form)
function _to_term_vector(ts::Vector{<:Term})
    return [(t.coefficient, t.monomial) for t in ts]
end

# Convert a Polynomial to vector format
function _to_term_vector(p::Polynomial)
    return [(t.coefficient, t.monomial) for t in p.terms]
end

# Recursive Cartesian product builder
function _cartesian_product!(
    result::Vector{Term},
    component_terms::Tuple,
    idx::Int,
    current_monomials::Tuple,
    current_coef::Number,
)
    if idx > length(component_terms)
        # Base case: all components processed
        cm = ComposedMonomial(current_monomials)
        push!(result, Term(current_coef, cm))
        return nothing
    end

    # Recursive case: iterate over all terms in current component
    for (coef, mono) in component_terms[idx]
        new_coef = current_coef * coef
        new_monomials = (current_monomials..., mono)
        _cartesian_product!(result, component_terms, idx + 1, new_monomials, new_coef)
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
