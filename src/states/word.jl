# =============================================================================
# StateSymbol - Canonical state expectation atom
# =============================================================================

"""
    StateSymbol{ST<:StateType, A<:AlgebraType, T<:Integer} <: AbstractMonomial

A single state expectation symbol wrapping a canonicalized monomial.

StateSymbol is the atomic unit of state expectations. It stores a single monomial
that has been canonicalized according to the state type:
- `Arbitrary`: involution canon (min(m, adjoint(m)))
- `MaxEntangled`: cyclic_symmetric_canon(m)

# Fields
- `mono::Monomial{A,T}`: The canonicalized monomial

# Type Parameters
- `ST`: State type (Arbitrary or MaxEntangled)
- `A`: Algebra type
- `T`: Integer type for monomial words

# Invariants
- The monomial is always in canonical form for the given state type
- Canonicalization happens automatically at construction

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> sym = StateSymbol{Arbitrary}(m);

julia> sym.mono == m  # Already canonical
true

julia> degree(sym)
2
```

See also: [`StateWord`](@ref), [`_state_canon`](@ref)
"""
struct StateSymbol{ST<:StateType,A<:AlgebraType,T<:Integer} <: AbstractMonomial
    mono::Monomial{A,T}

    # Inner constructor: canonicalize on construction
    function StateSymbol{ST}(m::Monomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
        canon_m = _state_canon(ST, m)
        new{ST,A,T}(canon_m)
    end
end

"""
    StateSymbol(::Type{ST}, m::Monomial{A,T}) where {ST<:StateType,A,T}

Construct a StateSymbol with explicit state type.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> sym = StateSymbol(Arbitrary, m);

julia> sym isa StateSymbol{Arbitrary}
true
```
"""
function StateSymbol(::Type{ST}, m::Monomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateSymbol{ST}(m)
end

# =============================================================================
# StateSymbol - Identity and One
# =============================================================================

"""
    Base.one(::Type{StateSymbol{ST,A,T}}) where {ST,A,T}

Create the identity StateSymbol (wrapping identity monomial).
"""
function Base.one(::Type{StateSymbol{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateSymbol{ST}(one(Monomial{A,T}))
end

"""
    Base.one(sym::StateSymbol{ST,A,T}) where {ST,A,T}

Create the identity StateSymbol for the same type.
"""
function Base.one(::StateSymbol{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    one(StateSymbol{ST,A,T})
end

"""
    Base.isone(sym::StateSymbol) -> Bool

Check if a StateSymbol wraps the identity monomial.
"""
Base.isone(sym::StateSymbol) = isone(sym.mono)

# =============================================================================
# StateSymbol - Degree and Variables
# =============================================================================

"""
    degree(sym::StateSymbol) -> Int

Compute the degree of a StateSymbol (degree of its monomial).
"""
degree(sym::StateSymbol) = degree(sym.mono)

"""
    variables(sym::StateSymbol{ST,A,T}) -> Set{T}

Get the set of variable indices used in the StateSymbol.
For signed algebras (fermionic/bosonic), returns absolute values.
"""
function variables(sym::StateSymbol{ST,A,T}) where {ST,A,T}
    result = Set{T}()
    for idx in sym.mono.word
        push!(result, abs(idx))
    end
    result
end

# =============================================================================
# StateSymbol - Comparison Operators
# =============================================================================

"""
    Base.:(==)(a::StateSymbol{ST,A,T}, b::StateSymbol{ST,A,T}) where {ST,A,T} -> Bool

Check if two StateSymbols are equal (same canonicalized monomial).
"""
function Base.:(==)(a::StateSymbol{ST,A,T}, b::StateSymbol{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    a.mono == b.mono
end

"""
    Base.hash(sym::StateSymbol, h::UInt) -> UInt

Hash function for StateSymbol.
"""
Base.hash(sym::StateSymbol, h::UInt) = hash(sym.mono, h)

"""
    Base.isless(a::StateSymbol{ST,A,T}, b::StateSymbol{ST,A,T}) where {ST,A,T} -> Bool

Compare two StateSymbols: degree-first, then by monomial ordering.
"""
function Base.isless(a::StateSymbol{ST,A,T}, b::StateSymbol{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    degree(a) != degree(b) && return degree(a) < degree(b)
    return isless(a.mono, b.mono)
end

# =============================================================================
# StateSymbol - Adjoint
# =============================================================================

"""
    Base.adjoint(sym::StateSymbol{ST,A,T}) where {ST,A,T}

Return the StateSymbol itself (due to canonicalization, adjoint equals self).
"""
function Base.adjoint(sym::StateSymbol{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    sym
end

# =============================================================================
# StateSymbol - Display
# =============================================================================

"""
    Base.show(io::IO, sym::StateSymbol{ST,A,T}) where {ST,A,T}

Display a StateSymbol with appropriate brackets based on state type.
- Arbitrary: ⟨mono⟩
- MaxEntangled: tr(mono)
"""
function Base.show(io::IO, sym::StateSymbol{ST,A,T}) where {ST,A,T}
    if ST == MaxEntangled
        prefix = "tr("
        suffix = ")"
    else  # Arbitrary
        prefix = "⟨"
        suffix = "⟩"
    end
    mono_str = sprint(show, sym.mono; context=io)
    print(io, prefix, mono_str, suffix)
end

# =============================================================================
# StateWord - Product of state expectations
# =============================================================================

"""
    StateWord{ST<:StateType, A<:AlgebraType, T<:Integer} <: AbstractMonomial

A product of state expectations <M1><M2>...<Mk>.

The state expectations commute, so StateWord maintains symbols in sorted order.
All symbols share the same algebra type A.

# Fields
- `state_syms::Vector{StateSymbol{ST,A,T}}`: Sorted, canonicalized state symbols

# Type Parameters
- `ST`: State type (Arbitrary or MaxEntangled)
- `A`: Algebra type (all expectations use the same algebra)
- `T`: Integer type for monomial words

# Invariants
1. **Canonicalization**: Each symbol is canonicalized per state type
2. **Commutativity**: state_syms is sorted by the symbol ordering
3. Identity symbols are filtered out (unless all are identity)

!!! warning "Real expectation values"
    The involution canonicalization enforces ⟨M⟩ = ⟨M†⟩, which means all expectation
    values are treated as real variables. This is appropriate for Hermitian moment
    optimization but restricts the variable space to real-valued expectations.

!!! warning "MaxEntangled with PauliAlgebra"
    The `MaxEntangled` (trace) state type uses cyclic-symmetric canonicalization,
    which assumes tr(M) = tr(reverse(M)). This is NOT valid for PauliAlgebra where
    transposition can introduce signs (e.g., tr(XYZ) ≠ tr(ZYX) for Pauli matrices).
    Use `MaxEntangled` only with algebras where trace is symmetric under reversal
    (e.g., NonCommutativeAlgebra, ProjectorAlgebra, UnipotentAlgebra).

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}([1, 2]);  # XY

julia> m2 = Monomial{PauliAlgebra}([3]);     # Z

julia> sw = StateWord{Arbitrary}([m1, m2]);  # <XY><Z>

julia> length(sw.state_syms)
2

julia> degree(sw)
3
```

See also: [`StateSymbol`](@ref), [`NCStateWord`](@ref), [`StatePolynomial`](@ref)
"""
struct StateWord{ST<:StateType,A<:AlgebraType,T<:Integer} <: AbstractMonomial
    state_syms::Vector{StateSymbol{ST,A,T}}

    # Inner constructor from symbols (already canonicalized)
    function StateWord{ST,A,T}(syms::Vector{StateSymbol{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
        # Filter out identity symbols (unless all are identity)
        filtered = filter(s -> !isone(s), syms)
        sorted_syms = isempty(filtered) ? [one(StateSymbol{ST,A,T})] : sort(filtered)
        new{ST,A,T}(sorted_syms)
    end
end

# Constructor from monomials (canonical path)
function StateWord{ST}(monos::Vector{Monomial{A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    syms = StateSymbol{ST,A,T}[StateSymbol{ST}(m) for m in monos]
    StateWord{ST,A,T}(syms)
end

# =============================================================================
# Canonicalization Helpers
# =============================================================================

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

# =============================================================================
# StateWord - Additional Constructors
# =============================================================================

# Convenience constructors

# Explicit type parameters from monomials
function StateWord{ST,A,T}(monos::Vector{Monomial{A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST}(monos)
end

# From vector of symbols (alternate path)
function StateWord{ST}(syms::Vector{StateSymbol{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST,A,T}(syms)
end

"""
    StateWord{ST}(m::Monomial{A,T}) where {ST,A,T}
    StateWord{ST}(monos::Vector{Monomial{A,T}}) where {ST,A,T}
    StateWord{ST}(sym::StateSymbol{ST,A,T}) where {ST,A,T}
    StateWord{ST}(syms::Vector{StateSymbol{ST,A,T}}) where {ST,A,T}

Construct a StateWord from monomials or StateSymbols.
Each monomial is lifted to a StateSymbol with canonicalization.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3]);

julia> sw = StateWord{Arbitrary}([m1, m2]);

julia> length(sw.state_syms)
2

julia> sw_single = StateWord{Arbitrary}(m1);

julia> length(sw_single.state_syms)
1
```
"""
function StateWord{ST}(m::Monomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST}([m])
end

# From single symbol
function StateWord{ST}(sym::StateSymbol{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST,A,T}([sym])
end

"""
    StateWord{ST,A,T}() where {ST,A,T}

Construct the identity StateWord (representing the constant 1).
Contains only the identity symbol.

# Examples
```jldoctest
julia> using NCTSSoS

julia> sw = StateWord{Arbitrary,PauliAlgebra,Int}();

julia> isone(sw)
true
```
"""
function StateWord{ST,A,T}() where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST,A,T}([one(StateSymbol{ST,A,T})])
end

# =============================================================================
# StateWord - Identity and Zero
# =============================================================================

"""
    Base.one(::Type{StateWord{ST,A,T}}) where {ST,A,T}

Create the identity StateWord (single identity symbol).

# Examples
```jldoctest
julia> using NCTSSoS

julia> sw_one = one(StateWord{Arbitrary,PauliAlgebra,Int64});

julia> isone(sw_one)
true
```
"""
function Base.one(::Type{StateWord{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST,A,T}([one(StateSymbol{ST,A,T})])
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

Check if a StateWord is the identity (single identity symbol).

# Examples
```jldoctest
julia> using NCTSSoS

julia> sw_one = one(StateWord{Arbitrary,PauliAlgebra,Int64});

julia> isone(sw_one)
true

julia> m = Monomial{PauliAlgebra}([1]);

julia> sw = StateWord{Arbitrary}([m]);

julia> isone(sw)
false
```
"""
Base.isone(sw::StateWord) = length(sw.state_syms) == 1 && isone(sw.state_syms[1])

# =============================================================================
# StateWord - Degree and Variables
# =============================================================================

"""
    degree(sw::StateWord) -> Int

Compute the total degree of a StateWord (sum of symbol degrees).

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}([1, 2]);  # degree 2

julia> m2 = Monomial{PauliAlgebra}([3]);     # degree 1

julia> sw = StateWord{Arbitrary}([m1, m2]);

julia> degree(sw)
3
```
"""
degree(sw::StateWord) = sum(degree, sw.state_syms; init=0)

"""
    variables(sw::StateWord) -> Set

Get the set of all variable indices used in the StateWord's symbols.

# Examples
```jldoctest
julia> using NCTSSoS

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
    for sym in sw.state_syms
        union!(result, variables(sym))
    end
    result
end

# =============================================================================
# StateWord - Comparison Operators
# =============================================================================

"""
    Base.:(==)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST,A,T} -> Bool

Check if two StateWords are equal (same state symbols).
"""
function Base.:(==)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    a.state_syms == b.state_syms
end

"""
    Base.hash(sw::StateWord, h::UInt) -> UInt

Hash function for StateWord. Computed from state_syms.
"""
Base.hash(sw::StateWord, h::UInt) = hash(sw.state_syms, h)

"""
    Base.isless(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST,A,T} -> Bool

Compare two StateWords using degree-first ordering, then lexicographic on state_syms.

# Examples
```jldoctest
julia> using NCTSSoS

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
    return a.state_syms < b.state_syms
end

# =============================================================================
# StateWord - Multiplication (commutative, no phase)
# =============================================================================

"""
    Base.:(*)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST,A,T}

Multiply two StateWords by concatenating and re-sorting their symbols.

State expectations commute, so the result is a sorted StateWord containing
all expectations from both operands.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3]);

julia> sw1 = StateWord{Arbitrary}([m1]);

julia> sw2 = StateWord{Arbitrary}([m2]);

julia> result = sw1 * sw2;

julia> result isa StateWord
true

julia> length(result.state_syms)
2
```
"""
function Base.:(*)(a::StateWord{ST,A,T}, b::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST,A,T}(vcat(a.state_syms, b.state_syms))
end

# =============================================================================
# Adjoint
# =============================================================================

"""
    Base.adjoint(sw::StateWord{ST,A,T}) where {ST,A,T}

Compute the adjoint (Hermitian conjugate) of a StateWord.
Due to the involution invariant, adjoint(sw) == sw for all StateWords.

!!! note "Physics notation"
    This is the dagger (†) or star (*) operation in physics notation.
    You can also use the Julia syntax `sw'` as shorthand for `adjoint(sw)`.

# Examples
```jldoctest
julia> using NCTSSoS

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

# =============================================================================
# StateWord - Display
# =============================================================================

"""
    Base.show(io::IO, sw::StateWord{ST,A,T}) where {ST,A,T}

Display a StateWord with appropriate brackets based on state type.
Uses registry from IO context if available for human-readable symbols.
"""
function Base.show(io::IO, sw::StateWord{ST,A,T}) where {ST,A,T}
    parts = String[]
    for sym in sw.state_syms
        # StateSymbol already handles display with proper brackets
        sym_str = sprint(show, sym; context=io)
        push!(parts, sym_str)
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

# Examples
```jldoctest
julia> using NCTSSoS

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

    function NCStateWord(sw::StateWord{ST,A,T}, nc_word::Monomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
        new{ST,A,T}(sw, nc_word)
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

Hash function for NCStateWord. Computed from sw and nc_word.
"""
Base.hash(ncsw::NCStateWord, h::UInt) = hash((ncsw.sw, ncsw.nc_word), h)

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

Compute the adjoint (Hermitian conjugate) of an NCStateWord by adjointing both parts.

!!! note "Physics notation"
    This is the dagger (†) or star (*) operation in physics notation.
    You can also use the Julia syntax `ncsw'` as shorthand for `adjoint(ncsw)`.
"""
function Base.adjoint(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStateWord(adjoint(ncsw.sw), adjoint(ncsw.nc_word))
end

# =============================================================================
# NCStateWord - Operations
# =============================================================================

"""
    neat_dot(x::NCStateWord, y::NCStateWord)

Compute adjoint(x) * y. Common operation in quantum optimization.

# Examples
```jldoctest
julia> using NCTSSoS

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
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3, 4]);

julia> sw = StateWord{Arbitrary}([m1]);

julia> ncsw = NCStateWord(sw, m2);

julia> ev = expval(ncsw);

julia> length(ev.state_syms)
2
```
"""
function expval(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Create a new symbol for the nc_word and concatenate with existing symbols
    new_sym = StateSymbol{ST}(ncsw.nc_word)
    StateWord{ST,A,T}(vcat(ncsw.sw.state_syms, [new_sym]))
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
# NCStateWord - Simplification
# =============================================================================

"""
    simplify(ncsw::NCStateWord) -> NCStatePolynomial

Simplify an NCStateWord by applying algebra-specific simplification rules to its
nc_word (Monomial) part. Returns an NCStatePolynomial since simplification may
produce multiple terms (e.g., Pauli algebra phase factors).

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = Monomial{UnipotentAlgebra}(UInt[1, 1]);  # x₁²

julia> sw = StateWord{Arbitrary}([one(Monomial{UnipotentAlgebra,UInt})]);

julia> ncsw = NCStateWord(sw, m);

julia> result = simplify(ncsw);

julia> result isa NCStatePolynomial
true
```
"""
function simplify(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Simplify the nc_word part (may produce multiple terms with phases)
    nc_poly = Polynomial(simplify(ncsw.nc_word))

    # Convert to NCStatePolynomial: each term gets the same StateWord
    coeffs = [t.coefficient for t in nc_poly.terms]
    ncsws = [NCStateWord(ncsw.sw, t.monomial) for t in nc_poly.terms]
    return NCStatePolynomial(coeffs, ncsws)
end

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
julia> using NCTSSoS

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> sw = ς(m);

julia> sw isa StateWord{Arbitrary}
true

julia> length(sw.state_syms)
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
julia> using NCTSSoS

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> sw = tr(m);

julia> sw isa StateWord{MaxEntangled}
true

julia> length(sw.state_syms)
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
julia> using NCTSSoS

julia> reg, (x,) = create_unipotent_variables([("x", 1:2)]);

julia> basis = get_state_basis(reg, 1);

julia> length(basis)  # includes <I>*I, <I>*x1, <I>*x2, <x1>*I, <x2>*I
5

julia> all(b -> b isa NCStateWord{Arbitrary}, basis)
true
```

For projector algebra:
```jldoctest
julia> using NCTSSoS

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
                    # Create compound StateWord by adding this expectation via symbol
                    new_sym = StateSymbol{ST}(m)
                    new_sw = StateWord{ST,A,T}(vcat(sw.state_syms, [new_sym]))
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
state symbol's monomial individually (symbols are already canonicalized
during StateSymbol construction).

Since StateWords represent products of expectations which commute, the overall
StateWord is already in a canonical sorted form. This function ensures each
constituent monomial is in its symmetric canonical form.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = Monomial{PauliAlgebra}([3, 2, 1]);

julia> sw = StateWord{Arbitrary}([m]);

julia> sw_canon = symmetric_canon(sw);

julia> sw_canon.state_syms[1].mono.word
3-element Vector{Int64}:
 1
 2
 3
```
"""
function symmetric_canon(sw::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Apply symmetric_canon to each monomial in the state symbols
    canon_monos = [symmetric_canon(sym.mono) for sym in sw.state_syms]
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
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> idx1 = encode_index(UInt8, 1, 1);

julia> m = Monomial{ProjectorAlgebra}([idx1, idx1]);  # P₁P₁

julia> sw = StateWord{MaxEntangled}([m]);

julia> sw_canon = symmetric_canon(sw);

julia> length(sw_canon.state_syms[1].mono.word)  # P₁P₁ → P₁
1
```
"""
function symmetric_canon(sw::StateWord{ST,ProjectorAlgebra,T}) where {ST<:StateType,T<:Integer}
    # For projector algebra, apply P²=P simplification before cyclic canonicalization
    canon_monos = Monomial{ProjectorAlgebra,T}[]
    for sym in sw.state_syms
        simplified = simplify(sym.mono)  # Returns Monomial for ProjectorAlgebra
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
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}([3, 2, 1]);

julia> m2 = Monomial{PauliAlgebra}([2, 1]);

julia> sw = StateWord{Arbitrary}([m1]);

julia> ncsw = NCStateWord(sw, m2);

julia> ncsw_canon = symmetric_canon(ncsw);

julia> ncsw_canon.sw.state_syms[1].mono.word
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
