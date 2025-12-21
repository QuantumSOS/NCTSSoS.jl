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
        # Apply canonicalization based on state type:
        # - MaxEntangled (trace): cyclic-symmetric canonicalization since tr(ABC) = tr(BCA) = tr(C†B†A†)
        # - Arbitrary: involution canonicalization (symmetric only) since <M> = <M†>
        canonicalized = [_state_canon(ST, m) for m in monos]

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

"""
    _state_canon(::Type{ST}, m::Monomial{A,T}) where {ST<:StateType,A,T}

Apply state-type-specific canonicalization to a monomial.

- For `Arbitrary`: involution canonicalization (min of m and adjoint(m))
- For `MaxEntangled`: cyclic-symmetric canonicalization (trace invariance)

# Arguments
- `ST`: The state type (Arbitrary or MaxEntangled)
- `m`: The monomial to canonicalize

# Returns
A canonicalized monomial.
"""
function _state_canon(::Type{Arbitrary}, m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # For arbitrary states: <M> = <M†>, so use involution (symmetric) canonicalization
    _involution_canon(m)
end

function _state_canon(::Type{MaxEntangled}, m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # For trace states: tr(ABC) = tr(BCA) = tr(CBA)† = tr(A†B†C†)
    # Use cyclic-symmetric canonicalization
    cyclic_symmetric_canon(m)
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

"""
    StateWord{ST,A,T}() where {ST,A,T}

Construct the identity StateWord (representing the constant 1).
Contains only the identity monomial.

# Examples
```jldoctest
julia> using FastPolynomials

julia> sw = StateWord{Arbitrary,PauliAlgebra,Int}();

julia> isone(sw)
true
```
"""
function StateWord{ST,A,T}() where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST}([one(Monomial{A,T})])
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
# Multiplication - Returns StateWord (commutative, no phase)
# =============================================================================

"""
    Base.:(*)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST,A,T}

Multiply two StateWords by concatenating and re-sorting their monomials.

State expectations commute, so the result is a sorted StateWord containing
all expectations from both operands.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3]);

julia> sw1 = StateWord{Arbitrary}([m1]);

julia> sw2 = StateWord{Arbitrary}([m2]);

julia> result = sw1 * sw2;

julia> result isa StateWord
true

julia> length(result.state_monos)
2
```
"""
function Base.:(*)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST}(vcat(a.state_monos, b.state_monos))
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
Uses registry from IO context if available for human-readable symbols.
"""
function Base.show(io::IO, sw::StateWord{ST,A,T}) where {ST,A,T}
    if ST == MaxEntangled
        prefix = "tr("
        suffix = ")"
    else  # Arbitrary or other
        prefix = "⟨"
        suffix = "⟩"
    end

    parts = String[]
    for m in sw.state_monos
        # Capture monomial display using the same IO context (with registry)
        mono_str = sprint(show, m; context=io)
        push!(parts, prefix * mono_str * suffix)
    end
    print(io, join(parts, ""))
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
    # nc_word multiplication returns Monomial (word concatenation only)
    nc_prod = a.nc_word * b.nc_word
    # StateWord multiplication returns StateWord (commutative, no phase)
    sw_prod = a.sw * b.sw
    NCStateWord(sw_prod, nc_prod)
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
Uses registry from IO context if available for human-readable symbols.
"""
function Base.show(io::IO, ncsw::NCStateWord{ST,A,T}) where {ST,A,T}
    # Print StateWord part (will use registry from context)
    show(io, ncsw.sw)

    # Only print nc_word if it's not identity
    if !isone(ncsw.nc_word)
        print(io, "⊗")
        show(io, ncsw.nc_word)
    end
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

# =============================================================================
# State Basis Generation
# =============================================================================

"""
    get_state_basis(registry::VariableRegistry{A,T}, d::Int;
                    state_type::Type{ST}=Arbitrary) where {A<:AlgebraType, T<:Integer, ST<:StateType}

Generate a basis of NCStateWord elements up to degree d.

This function generates all unique NCStateWord basis elements that can be formed
from the variables in the registry up to the specified degree. The basis includes
all (StateWord, Monomial) combinations where `degree(sw) + degree(nc_word) <= d`.

This generates:
- `<I>*M` forms (identity StateWord, operator monomial)
- `<M>*I` forms (single expectation, identity operator)
- `<M1><M2>*I` forms (compound expectations, identity operator)
- `<M>*N` mixed forms (expectation with operator)

# Arguments
- `registry`: Variable registry containing the variable indices
- `d`: Maximum total degree (inclusive)

# Keyword Arguments
- `state_type`: The state type for the basis elements (default: `Arbitrary`)

# Returns
- `Vector{NCStateWord{ST,A,T}}`: Sorted unique NCStateWord basis elements

# Examples
```jldoctest
julia> using FastPolynomials

julia> reg, (x,) = create_unipotent_variables([("x", 1:2)]);

julia> basis = get_state_basis(reg, 1);

julia> length(basis)  # includes <I>*I, <I>*x1, <I>*x2, <x1>*I, <x2>*I
5

julia> all(b -> b isa NCStateWord{Arbitrary}, basis)
true
```

For projector algebra:
```jldoctest
julia> using FastPolynomials

julia> reg, (P,) = create_projector_variables([("P", 1:2)]);

julia> basis = get_state_basis(reg, 2; state_type=MaxEntangled);

julia> all(b -> b isa NCStateWord{MaxEntangled}, basis)
true
```

See also: [`NCStateWord`](@ref), [`get_ncbasis`](@ref)
"""
function get_state_basis(
    registry::VariableRegistry{A,T},
    d::Int;
    state_type::Type{ST}=Arbitrary
) where {A<:AlgebraType, T<:Integer, ST<:StateType}
    
    # Collect monomials by degree
    monos_by_deg = Vector{Vector{Monomial{A,T}}}(undef, d + 1)
    for deg in 0:d
        poly_basis = get_ncbasis(registry, deg)
        monos = Monomial{A,T}[]
        for poly in poly_basis
            for term in poly.terms
                push!(monos, term.monomial)
            end
        end
        unique!(sort!(monos))
        monos_by_deg[deg + 1] = monos
    end
    
    # All monomials up to degree d
    all_monos = reduce(vcat, monos_by_deg)
    unique!(sort!(all_monos))
    
    result = NCStateWord{ST,A,T}[]
    
    # Generate all (StateWord, Monomial) combinations with total degree <= d
    for nc_deg in 0:d
        nc_monos = monos_by_deg[nc_deg + 1]
        sw_max_deg = d - nc_deg
        
        # Get all StateWords up to sw_max_deg
        sw_basis = _generate_statewords_up_to_degree(all_monos, sw_max_deg, ST)
        
        # Create NCStateWords for all valid combinations
        for sw in sw_basis
            for nc_m in nc_monos
                if degree(sw) + degree(nc_m) <= d
                    push!(result, NCStateWord(sw, nc_m))
                end
            end
        end
    end
    
    unique!(sort!(result))
    return result
end

"""
    _generate_statewords_up_to_degree(monos, max_deg, ST) -> Vector{StateWord}

Generate all StateWords with total degree <= max_deg.

Includes:
- Identity StateWord <I> (degree 0)
- Single expectations <M>
- Compound expectations <M1><M2>, <M1><M2><M3>, etc.

# Arguments
- `monos`: Vector of all available monomials
- `max_deg`: Maximum total degree for the StateWord
- `ST`: StateType (Arbitrary or MaxEntangled)

# Returns
Vector of unique, sorted StateWords.
"""
function _generate_statewords_up_to_degree(
    monos::Vector{Monomial{A,T}},
    max_deg::Int,
    ::Type{ST}
) where {A<:AlgebraType, T<:Integer, ST<:StateType}
    
    result = StateWord{ST,A,T}[]
    
    # Always include identity StateWord
    push!(result, one(StateWord{ST,A,T}))
    
    if max_deg == 0
        return result
    end
    
    # Filter to non-identity monomials with degree <= max_deg
    valid_monos = filter(m -> !isone(m) && degree(m) <= max_deg, monos)
    
    isempty(valid_monos) && return result
    
    # Use BFS to build StateWords level by level
    # Each level contains StateWords with one more expectation than previous
    current_level = Tuple{StateWord{ST,A,T}, Int}[]
    
    # Start with single expectations
    for m in valid_monos
        if degree(m) <= max_deg
            sw = StateWord{ST}([m])
            push!(current_level, (sw, degree(m)))
            push!(result, sw)
        end
    end
    
    # Build compound expectations iteratively
    while !isempty(current_level)
        next_level = Tuple{StateWord{ST,A,T}, Int}[]
        seen = Set{StateWord{ST,A,T}}()
        
        for (sw, sw_deg) in current_level
            for m in valid_monos
                new_deg = sw_deg + degree(m)
                if new_deg <= max_deg
                    # Create compound StateWord by adding this expectation
                    new_sw = StateWord{ST}(vcat(sw.state_monos, [m]))
                    if !(new_sw in seen)
                        push!(seen, new_sw)
                        push!(next_level, (new_sw, new_deg))
                        push!(result, new_sw)
                    end
                end
            end
        end
        
        current_level = next_level
    end
    
    unique!(result)
    return result
end

# =============================================================================
# StateWord Canonicalization
# =============================================================================
#
# Import symmetric_canon from canonicalization.jl (included before this file)
# and extend it for StateWord types.

"""
    symmetric_canon(sw::StateWord{ST,A,T}) where {ST,A,T}

Return a new StateWord with symmetrically canonicalized monomials.

For StateWords, symmetric canonicalization applies `symmetric_canon` to each
state monomial individually (state monomials are already involution-canonicalized
during StateWord construction via `_involution_canon`).

Since StateWords represent products of expectations which commute, the overall
StateWord is already in a canonical sorted form. This function ensures each
constituent monomial is in its symmetric canonical form.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([3, 2, 1]);

julia> sw = StateWord{Arbitrary}([m]);

julia> sw_canon = symmetric_canon(sw);

julia> sw_canon.state_monos[1].word
3-element Vector{Int64}:
 1
 2
 3
```
"""
function symmetric_canon(sw::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Apply symmetric_canon to each monomial in the state word
    canon_monos = [symmetric_canon(m) for m in sw.state_monos]
    StateWord{ST}(canon_monos)
end

"""
    symmetric_canon(sw::StateWord{ST,ProjectorAlgebra,T}) where {ST,T}

Specialized symmetric canonicalization for ProjectorAlgebra StateWords.

For projector algebra, this applies the idempotency rule (P² = P) via `simplify`
before cyclic canonicalization. This ensures that equivalent StateWords like
`<P₁P₁P₂>` and `<P₁P₂>` are recognized as identical.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1 = encode_index(UInt8, 1, 1);

julia> m = Monomial{ProjectorAlgebra}([idx1, idx1]);  # P₁P₁

julia> sw = StateWord{MaxEntangled}([m]);

julia> sw_canon = symmetric_canon(sw);

julia> length(sw_canon.state_monos[1].word)  # P₁P₁ → P₁
1
```
"""
function symmetric_canon(sw::StateWord{ST,ProjectorAlgebra,T}) where {ST<:StateType,T<:Integer}
    # For projector algebra, apply P²=P simplification before cyclic canonicalization
    canon_monos = Monomial{ProjectorAlgebra,T}[]
    for m in sw.state_monos
        simplified = simplify(m)  # Returns Monomial for ProjectorAlgebra
        # Only keep non-identity monomials
        if !isone(simplified)
            push!(canon_monos, symmetric_canon(simplified))
        end
    end
    StateWord{ST}(canon_monos)
end

"""
    symmetric_canon(ncsw::NCStateWord{ST,A,T}) where {ST,A,T}

Return a new NCStateWord with symmetrically canonicalized components.

For NCStateWords, this canonicalizes both the StateWord part and the nc_word part.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([3, 2, 1]);

julia> m2 = Monomial{PauliAlgebra}([2, 1]);

julia> sw = StateWord{Arbitrary}([m1]);

julia> ncsw = NCStateWord(sw, m2);

julia> ncsw_canon = symmetric_canon(ncsw);

julia> ncsw_canon.sw.state_monos[1].word
3-element Vector{Int64}:
 1
 2
 3

julia> ncsw_canon.nc_word.word
2-element Vector{Int64}:
 1
 2
```
"""
function symmetric_canon(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStateWord(symmetric_canon(ncsw.sw), symmetric_canon(ncsw.nc_word))
end
