"""
    StateWord{ST<:StateType, A<:AlgebraType, T<:Integer} <: AbstractMonomial

A product of state expectations <M1><M2>...<Mk>.

The state expectations commute, so StateWord maintains monomials in sorted order.
All monomials share the same algebra type A.

# Fields
- `state_monos::Vector{Monomial{A,T}}`: Sorted, canonicalized monomials (expectations)
- `hash::UInt64`: Precomputed hash for fast equality checks

# Type Parameters
- `ST`: State type (Arbitrary or MaxEntangled)
- `A`: Algebra type (all expectations use the same algebra)
- `T`: Integer type for monomial words

# Invariants
1. **Involution**: Each monomial is canonicalized to min(m, adjoint(m))
2. **Commutativity**: state_monos is sorted by the monomial ordering
3. Identity monomials are filtered out (unless all are identity)
4. Hash is precomputed for fast equality checks

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);  # XY

julia> m2 = Monomial{PauliAlgebra}([3]);     # Z

julia> sw = StateWord{Arbitrary}([m1, m2]);  # <XY><Z>

julia> length(sw.state_monos)
2

julia> degree(sw)
3
```

See also: [`NCStateWord`](@ref), [`StatePolynomial`](@ref)
"""
struct StateWord{ST<:StateType,A<:AlgebraType,T<:Integer} <: AbstractMonomial
    state_monos::Vector{Monomial{A,T}}
    hash::UInt64

    function StateWord{ST}(monos::Vector{Monomial{A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
        # Apply involution invariant: canonicalize each monomial to min(m, adjoint(m))
        canonicalized = [_involution_canon(m) for m in monos]

        # Filter out identity monomials (unless all are identity)
        filtered = filter(m -> !isone(m), canonicalized)
        sorted_monos = isempty(filtered) ? [one(Monomial{A,T})] : sort(filtered)
        h = hash(sorted_monos)
        new{ST,A,T}(sorted_monos, h)
    end
end

"""
    _involution_canon(m::Monomial{A,T}) where {A,T}

Apply the involution canonicalization: return min(m, adjoint(m)).
This ensures that m and adjoint(m) are treated as equivalent in state expectations.
"""
function _involution_canon(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    m_adj = adjoint(m)
    return isless(m, m_adj) ? m : (m == m_adj ? m : m_adj)
end

# Convenience constructors
"""
    StateWord{ST,A,T}(monos::Vector{Monomial{A,T}}) where {ST,A,T}

Construct a StateWord with explicit type parameters.
"""
function StateWord{ST,A,T}(monos::Vector{Monomial{A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST}(monos)
end

"""
    StateWord{ST}(m::Monomial{A,T}) where {ST,A,T}

Construct a StateWord from a single monomial.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> sw = StateWord{Arbitrary}(m);

julia> length(sw.state_monos)
1
```
"""
function StateWord{ST}(m::Monomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST}([m])
end

# =============================================================================
# Identity and Zero
# =============================================================================

"""
    Base.one(::Type{StateWord{ST,A,T}}) where {ST,A,T}

Create the identity StateWord (single identity monomial).

# Examples
```jldoctest
julia> using FastPolynomials

julia> sw_one = one(StateWord{Arbitrary,PauliAlgebra,Int64});

julia> isone(sw_one)
true
```
"""
function Base.one(::Type{StateWord{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST}([one(Monomial{A,T})])
end

"""
    Base.one(sw::StateWord{ST,A,T}) where {ST,A,T}

Create the identity StateWord for the same type as `sw`.
"""
function Base.one(::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    one(StateWord{ST,A,T})
end

"""
    Base.isone(sw::StateWord) -> Bool

Check if a StateWord is the identity (single identity monomial).

# Examples
```jldoctest
julia> using FastPolynomials

julia> sw_one = one(StateWord{Arbitrary,PauliAlgebra,Int64});

julia> isone(sw_one)
true

julia> m = Monomial{PauliAlgebra}([1]);

julia> sw = StateWord{Arbitrary}([m]);

julia> isone(sw)
false
```
"""
Base.isone(sw::StateWord) = length(sw.state_monos) == 1 && isone(sw.state_monos[1])

# =============================================================================
# Degree and Variables
# =============================================================================

"""
    degree(sw::StateWord) -> Int

Compute the total degree of a StateWord (sum of monomial degrees).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);  # degree 2

julia> m2 = Monomial{PauliAlgebra}([3]);     # degree 1

julia> sw = StateWord{Arbitrary}([m1, m2]);

julia> degree(sw)
3
```
"""
degree(sw::StateWord) = sum(degree, sw.state_monos; init=0)

"""
    variables(sw::StateWord) -> Set

Get the set of all variable indices used in the StateWord's monomials.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([2, 3]);

julia> sw = StateWord{Arbitrary}([m1, m2]);

julia> variables(sw)
Set{Int64} with 3 elements:
  2
  3
  1
```
"""
function variables(sw::StateWord{ST,A,T}) where {ST,A,T}
    result = Set{T}()
    for m in sw.state_monos
        for idx in m.word
            push!(result, abs(idx))  # abs for fermionic/bosonic (negative = creation)
        end
    end
    result
end

# =============================================================================
# Comparison Operators
# =============================================================================

"""
    Base.:(==)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST,A,T} -> Bool

Check if two StateWords are equal (same state monomials).
"""
function Base.:(==)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    a.state_monos == b.state_monos
end

"""
    Base.hash(sw::StateWord, h::UInt) -> UInt

Hash function for StateWord. Uses precomputed hash.
"""
Base.hash(sw::StateWord, h::UInt) = hash(sw.hash, h)

"""
    Base.isless(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST,A,T} -> Bool

Compare two StateWords using degree-first ordering, then lexicographic on state_monos.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([1, 2]);

julia> sw1 = StateWord{Arbitrary}([m1]);

julia> sw2 = StateWord{Arbitrary}([m2]);

julia> isless(sw1, sw2)  # degree 1 < degree 2
true
```
"""
function Base.isless(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    degree(a) != degree(b) && return degree(a) < degree(b)
    return a.state_monos < b.state_monos
end

# =============================================================================
# Multiplication - Returns Term for consistency with other AbstractMonomial types
# =============================================================================

"""
    Base.:(*)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST,A,T}

Multiply two StateWords by concatenating and re-sorting their monomials.
Returns a Term{StateWord, Float64} for consistency with other monomial algebras.

State expectations commute, so the result is sorted.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3]);

julia> sw1 = StateWord{Arbitrary}([m1]);

julia> sw2 = StateWord{Arbitrary}([m2]);

julia> result = sw1 * sw2;

julia> result isa Term
true

julia> length(result.monomial.state_monos)
2
```
"""
function Base.:(*)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    product_sw = StateWord{ST}(vcat(a.state_monos, b.state_monos))
    Term(1.0, product_sw)
end

# =============================================================================
# Adjoint / Star
# =============================================================================

"""
    Base.adjoint(sw::StateWord{ST,A,T}) where {ST,A,T}

Compute the adjoint of a StateWord by adjointing each monomial.
Due to the involution invariant, adjoint(sw) == sw for all StateWords.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2, 3]);

julia> sw = StateWord{Arbitrary}([m]);

julia> sw_adj = adjoint(sw);

julia> sw == sw_adj  # Due to involution canonicalization
true
```
"""
function Base.adjoint(sw::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Because of involution canonicalization, adjoint(m) maps to same canonical form as m
    # So adjoint of StateWord is the StateWord itself
    sw
end

"""
    star(sw::StateWord)

Alias for `adjoint`. Compute the star (dagger) of a StateWord.
"""
star(sw::StateWord) = adjoint(sw)

# =============================================================================
# Display
# =============================================================================

"""
    Base.show(io::IO, sw::StateWord{ST,A,T}) where {ST,A,T}

Display a StateWord with appropriate brackets based on state type.
"""
function Base.show(io::IO, sw::StateWord{ST,A,T}) where {ST,A,T}
    if ST == MaxEntangled
        prefix = "tr("
        suffix = ")"
    else  # Arbitrary or other
        prefix = "<"
        suffix = ">"
    end
    parts = [prefix * string(m.word) * suffix for m in sw.state_monos]
    print(io, join(parts, " "))
end

# =============================================================================
# NCStateWord
# =============================================================================

"""
    NCStateWord{ST<:StateType, A<:AlgebraType, T<:Integer}

A product of state expectations and a non-commutative operator: <M1><M2>...<Mk> * Onc

Combines a commutative StateWord (expectations) with a non-commutative monomial (operator).

# Fields
- `sw::StateWord{ST,A,T}`: Commutative state word part
- `nc_word::Monomial{A,T}`: Non-commutative operator part
- `hash::UInt64`: Precomputed hash

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3]);

julia> sw = StateWord{Arbitrary}([m1]);  # <XY>

julia> ncsw = NCStateWord(sw, m2);       # <XY> * Z

julia> degree(ncsw)
3
```

See also: [`StateWord`](@ref), [`NCStatePolynomial`](@ref), [`expval`](@ref)
"""
struct NCStateWord{ST<:StateType,A<:AlgebraType,T<:Integer}
    sw::StateWord{ST,A,T}
    nc_word::Monomial{A,T}
    hash::UInt64

    function NCStateWord(sw::StateWord{ST,A,T}, nc_word::Monomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
        h = hash((sw.hash, nc_word.hash))
        new{ST,A,T}(sw, nc_word, h)
    end
end

# Convenience constructors
"""
    NCStateWord{ST}(sw::StateWord{ST,A,T}, nc_word::Monomial{A,T}) where {ST,A,T}

Construct an NCStateWord with explicit state type.
"""
function NCStateWord{ST}(sw::StateWord{ST,A,T}, nc_word::Monomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStateWord(sw, nc_word)
end

"""
    NCStateWord(sw::StateWord{ST,A,T}) where {ST,A,T}

Construct an NCStateWord from a StateWord with identity nc_word.
"""
function NCStateWord(sw::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStateWord(sw, one(Monomial{A,T}))
end

"""
    NCStateWord(m::Monomial{A,T}) where {A,T}

Construct an NCStateWord from a monomial with identity StateWord.
"""
function NCStateWord(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    NCStateWord(one(StateWord{Arbitrary,A,T}), m)
end

# =============================================================================
# NCStateWord - Degree and Variables
# =============================================================================

"""
    degree(ncsw::NCStateWord) -> Int

Compute the total degree of an NCStateWord (sw degree + nc_word degree).
"""
degree(ncsw::NCStateWord) = degree(ncsw.sw) + degree(ncsw.nc_word)

"""
    variables(ncsw::NCStateWord) -> Set

Get the set of all variable indices used in the NCStateWord.
"""
function variables(ncsw::NCStateWord{ST,A,T}) where {ST,A,T}
    sw_vars = variables(ncsw.sw)
    for idx in ncsw.nc_word.word
        push!(sw_vars, abs(idx))
    end
    sw_vars
end

# =============================================================================
# NCStateWord - Comparison Operators
# =============================================================================

"""
    Base.:(==)(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST,A,T} -> Bool

Check if two NCStateWords are equal.
"""
function Base.:(==)(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    a.sw == b.sw && a.nc_word == b.nc_word
end

"""
    Base.hash(ncsw::NCStateWord, h::UInt) -> UInt

Hash function for NCStateWord.
"""
Base.hash(ncsw::NCStateWord, h::UInt) = hash(ncsw.hash, h)

"""
    Base.isless(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST,A,T} -> Bool

Compare two NCStateWords: degree first, then nc_word, then sw.
"""
function Base.isless(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    degree(a) != degree(b) && return degree(a) < degree(b)
    a.nc_word != b.nc_word && return isless(a.nc_word, b.nc_word)
    return isless(a.sw, b.sw)
end

# =============================================================================
# NCStateWord - Multiplication
# =============================================================================

"""
    Base.:(*)(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST,A,T}

Multiply two NCStateWords: multiply sw parts (commutative) and nc_word parts (non-commutative).
"""
function Base.:(*)(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # nc_word multiplication uses monomial multiplication (returns Term)
    nc_prod = a.nc_word * b.nc_word
    # StateWord multiplication now returns Term, extract the monomial
    sw_prod = a.sw * b.sw
    NCStateWord(sw_prod.monomial, nc_prod.monomial)
end

# =============================================================================
# NCStateWord - Adjoint / Star
# =============================================================================

"""
    Base.adjoint(ncsw::NCStateWord{ST,A,T}) where {ST,A,T}

Compute the adjoint of an NCStateWord by adjointing both parts.
"""
function Base.adjoint(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStateWord(adjoint(ncsw.sw), adjoint(ncsw.nc_word))
end

"""
    star(ncsw::NCStateWord)

Alias for `adjoint`. Compute the star (dagger) of an NCStateWord.
"""
star(ncsw::NCStateWord) = adjoint(ncsw)

# =============================================================================
# NCStateWord - Operations
# =============================================================================

"""
    neat_dot(x::NCStateWord, y::NCStateWord)

Compute adjoint(x) * y. Common operation in quantum optimization.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1]);

julia> sw = StateWord{Arbitrary}([m]);

julia> ncsw = NCStateWord(sw, m);

julia> nd = neat_dot(ncsw, ncsw);

julia> nd == adjoint(ncsw) * ncsw
true
```
"""
neat_dot(x::NCStateWord, y::NCStateWord) = adjoint(x) * y

"""
    expval(ncsw::NCStateWord{ST,A,T}) where {ST,A,T}

Collapse NCStateWord to StateWord by adding nc_word as an expectation.

Converts `<M1><M2>...<Mk> * Onc` to `<M1><M2>...<Mk><Onc>`.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3, 4]);

julia> sw = StateWord{Arbitrary}([m1]);

julia> ncsw = NCStateWord(sw, m2);

julia> ev = expval(ncsw);

julia> length(ev.state_monos)
2
```
"""
function expval(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST}(vcat(ncsw.sw.state_monos, [ncsw.nc_word]))
end

# =============================================================================
# NCStateWord - Display
# =============================================================================

"""
    Base.show(io::IO, ncsw::NCStateWord{ST,A,T}) where {ST,A,T}

Display an NCStateWord.
"""
function Base.show(io::IO, ncsw::NCStateWord{ST,A,T}) where {ST,A,T}
    print(io, ncsw.sw, " * ", ncsw.nc_word.word)
end

# =============================================================================
# NCStateWord - one
# =============================================================================

"""
    Base.one(::Type{NCStateWord{ST,A,T}}) where {ST,A,T}

Create the identity NCStateWord.
"""
function Base.one(::Type{NCStateWord{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStateWord(one(StateWord{ST,A,T}), one(Monomial{A,T}))
end

"""
    Base.one(ncsw::NCStateWord{ST,A,T}) where {ST,A,T}

Create the identity NCStateWord for the same type.
"""
function Base.one(::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    one(NCStateWord{ST,A,T})
end

"""
    Base.isone(ncsw::NCStateWord) -> Bool

Check if an NCStateWord is the identity.
"""
Base.isone(ncsw::NCStateWord) = isone(ncsw.sw) && isone(ncsw.nc_word)

# =============================================================================
# Convenience Constructors: ς and tr
# =============================================================================

"""
    ς(m::Monomial{A,T}) where {A,T}

Create a StateWord{Arbitrary} from a monomial.

This is a convenience function for creating state expectations in the
arbitrary state formalism. Equivalent to `StateWord{Arbitrary}(m)`.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> sw = ς(m);

julia> sw isa StateWord{Arbitrary}
true

julia> length(sw.state_monos)
1
```

See also: [`tr`](@ref), [`StateWord`](@ref)
"""
ς(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer} = StateWord{Arbitrary}(m)

"""
    tr(m::Monomial{A,T}) where {A,T}

Create a StateWord{MaxEntangled} from a monomial.

This is a convenience function for creating trace expressions in the
maximally entangled state formalism. Equivalent to `StateWord{MaxEntangled}(m)`.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> sw = tr(m);

julia> sw isa StateWord{MaxEntangled}
true

julia> length(sw.state_monos)
1
```

See also: [`ς`](@ref), [`StateWord`](@ref)
"""
tr(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer} = StateWord{MaxEntangled}(m)
