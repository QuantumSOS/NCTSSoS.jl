"""
    StateSymbol{ST<:StateType, A<:AlgebraType, T<:Integer} <: AbstractMonomial{A,T}

A single state expectation symbol wrapping a canonicalized monomial.

StateSymbol is the atomic unit of state expectations. It stores a single monomial
that has been canonicalized according to the state type:
- `Arbitrary`: involution canon (min(m, adjoint(m)))
- `MaxEntangled`: cyclic_symmetric_canon(m)

# Fields
- `mono::Vector{T}`: The canonicalized monomial

# Type Parameters
- `ST`: State type (Arbitrary or MaxEntangled)
- `A`: Algebra type
- `T`: Integer type for monomial words

# Invariants
- The monomial is always in canonical form for the given state type
- Canonicalization happens automatically at construction

# Note on non-monoid algebras
For algebras where simplification can introduce phases (Pauli / `TwistedGroupAlgebra`) or
multi-term PBW expansions (Fermionic/Bosonic / `PBWAlgebra`), the expectation value of a
`Monomial` is not representable by a single `StateSymbol`. Use `expval(ST, m::Monomial)`
to obtain a `StatePolynomial` (a linear combination of `StateWord`s) instead.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{ProjectorAlgebra}(UInt8[1, 2]);

julia> sym = StateSymbol{Arbitrary}(m);

julia> degree(sym)
2
```

See also: [`StateWord`](@ref), [`_state_canon`](@ref)
"""
struct StateSymbol{ST<:StateType,A<:AlgebraType,T<:Integer} <: AbstractMonomial{A,T}
    mono::Vector{T}

    # Inner constructor: canonicalize on construction
    function StateSymbol{ST}(m::NormalMonomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
        canon_w = canonicalize(ST,m)
        new{ST,A,T}(canon_w)
    end
end

# =============================================================================
# StateSymbol - Identity and One
# =============================================================================

"""
    Base.one(::Type{StateSymbol{ST,A,T}}) where {ST,A,T}

Create the identity StateSymbol (wrapping identity monomial).
"""
function Base.one(::Type{StateSymbol{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateSymbol{ST}(one(NormalMonomial{A,T}))
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
Base.isone(sym::StateSymbol) = isempty(sym.mono)

# =============================================================================
# StateSymbol - Degree and Variables
# =============================================================================

"""
    degree(sym::StateSymbol) -> Int

Compute the degree of a StateSymbol (degree of its monomial).
"""
degree(sym::StateSymbol) = length(sym.mono)

"""
    variables(sym::StateSymbol{ST,A,T}) -> Set{T}

Get the set of variable indices used in the StateSymbol.
"""
function variables(sym::StateSymbol{ST,A,T}) where {ST,A,T}
    result = Set{T}()
    for idx in sym.mono
        push!(result, idx)
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
    return a.mono < b.mono
end

# =============================================================================
# StateSymbol - Adjoint
# =============================================================================

"""
    Base.adjoint(sym::StateSymbol{ST,A,T}) where {ST,A,T}

Return the StateSymbol itself (due to canonicalization, adjoint equals self).
"""
Base.adjoint(sym::StateSymbol{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer} = sym
Base.adjoint(sym::StateSymbol{ST,A,T}) where {ST<:StateType,A<:TwistedGroupAlgebra,T<:Integer} = sym
Base.adjoint(sym::StateSymbol{ST,A,T}) where {ST<:StateType,A<:PBWAlgebra,T<:Integer} = throw(error("Adjoint not defined for PBWAlgebra"))


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
    StateWord{ST<:StateType, A<:AlgebraType, T<:Integer} <: AbstractMonomial{A,T}

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

!!! warning "MaxEntangled (trace) assumptions"
    The `MaxEntangled` (trace) state type uses cyclic-symmetric canonicalization,
    which assumes tr(M) = tr(reverse(M)). This can fail for algebras where
    transposition introduces signs/phases (e.g. Pauli algebra). `NCTSSoS.jl`
    restricts state polynomials to `MonoidAlgebra` (NC/Projector/Unipotent),
    so `PauliAlgebra` is not supported here.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

julia> m2 = NormalMonomial{NonCommutativeAlgebra}([3]);

julia> sw = StateWord{Arbitrary}([m1, m2]);

julia> length(sw.state_syms)
2

julia> degree(sw)
3
```

See also: [`StateSymbol`](@ref), [`NCStateWord`](@ref), [`StatePolynomial`](@ref)
"""
struct StateWord{ST<:StateType,A<:AlgebraType,T<:Integer} <: AbstractMonomial{A,T}
    state_syms::Vector{StateSymbol{ST,A,T}}

    # Inner constructor from symbols (already canonicalized)
    function StateWord{ST,A,T}(syms::Vector{StateSymbol{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
        # Filter out identity symbols (unless all are identity)
        filtered = filter(s -> !isone(s), syms)
        sorted_syms = isempty(filtered) ? [one(StateSymbol{ST,A,T})] : sort!(filtered)
        new{ST,A,T}(sorted_syms)
    end
end

StateWord{ST,A,T}() where {ST<:StateType,A<:AlgebraType,T<:Integer} = one(StateWord{ST,A,T})

StateWord{ST}(sw::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer} = sw

StateWord{ST}(sym::StateSymbol{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer} =
    StateWord{ST,A,T}([sym])

StateWord{ST}(syms::AbstractVector{StateSymbol{ST,A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer} =
    StateWord{ST,A,T}(collect(syms))

StateWord{ST}(m::NormalMonomial{A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer} =
    StateWord{ST,A,T}([StateSymbol{ST}(m)])

StateWord{ST}(monos::AbstractVector{NormalMonomial{A,T}}) where {ST<:StateType,A<:AlgebraType,T<:Integer} =
    StateWord{ST,A,T}([StateSymbol{ST}(m) for m in monos])

# =============================================================================
# StateWord - Identity and Zero
# =============================================================================

"""
    Base.one(::Type{StateWord{ST,A,T}}) where {ST,A,T}

Create the identity StateWord (single identity symbol).

# Examples
```jldoctest
julia> using NCTSSoS

julia> sw_one = one(StateWord{Arbitrary,NonCommutativeAlgebra,Int64});

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

julia> sw_one = one(StateWord{Arbitrary,NonCommutativeAlgebra,Int64});

julia> isone(sw_one)
true

julia> m = NormalMonomial{NonCommutativeAlgebra}([1]);

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

julia> m1 = NormalMonomial{NonCommutativeAlgebra}([1, 2]);  # degree 2

julia> m2 = NormalMonomial{NonCommutativeAlgebra}([3]);     # degree 1

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

julia> m1 = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

julia> m2 = NormalMonomial{NonCommutativeAlgebra}([2, 3]);

julia> sw = StateWord{Arbitrary}([m1, m2]);

julia> sort(collect(variables(sw)))
3-element Vector{Int64}:
 1
 2
 3
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

julia> m1 = NormalMonomial{NonCommutativeAlgebra}([1]);

julia> m2 = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

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
# StateSymbol - Multiplication (lifts to StateWord)
# =============================================================================

function Base.:(*)(a::StateSymbol{ST,A,T}, b::StateSymbol{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST,A,T}([a, b])
end

function Base.:(*)(a::StateSymbol{ST,A,T}, b::StateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST,A,T}(vcat([a], b.state_syms))
end

function Base.:(*)(a::StateWord{ST,A,T}, b::StateSymbol{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    StateWord{ST,A,T}(vcat(a.state_syms, [b]))
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

julia> m1 = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

julia> m2 = NormalMonomial{NonCommutativeAlgebra}([3]);

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

julia> m = NormalMonomial{NonCommutativeAlgebra}([1, 2, 3]);

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
    NCStateWord{ST<:StateType, A<:MonoidAlgebra, T<:Integer}

A product of state expectations and a non-commutative operator: <M1><M2>...<Mk> * Onc

Combines a commutative StateWord (expectations) with a non-commutative monomial (operator).

# Fields
- `sw::StateWord{ST,A,T}`: Commutative state word part
- `nc_word::NormalMonomial{A,T}`: Non-commutative operator part

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

julia> m2 = NormalMonomial{NonCommutativeAlgebra}([3]);

julia> sw = StateWord{Arbitrary}([m1]);

julia> ncsw = NCStateWord(sw, m2);

julia> degree(ncsw)
3
```

See also: [`StateWord`](@ref), [`NCStatePolynomial`](@ref), [`expval`](@ref)
"""
struct NCStateWord{ST<:StateType,A<:MonoidAlgebra,T<:Integer}
    sw::StateWord{ST,A,T}
    nc_word::NormalMonomial{A,T}

    function NCStateWord(sw::StateWord{ST,A,T}, nc_word::NormalMonomial{A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer}
        new{ST,A,T}(sw, nc_word)
    end
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
function Base.:(==)(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer}
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
function Base.isless(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer}
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
function Base.:(*)(a::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer}
    # NormalMonomial * NormalMonomial returns a Polynomial; for MonoidAlgebra this must be a single term.
    nc_prod_poly = a.nc_word * b.nc_word
    length(nc_prod_poly.terms) == 1 || throw(ArgumentError("NCStateWord requires monoid-algebra monomial products to stay single-term"))
    (c, nc_prod) = nc_prod_poly.terms[1]
    isone(c) || throw(ArgumentError("NCStateWord cannot represent non-unit coefficients (got $c)"))
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
function Base.adjoint(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer}
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

julia> m = NormalMonomial{NonCommutativeAlgebra}([1]);

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

julia> m1 = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

julia> m2 = NormalMonomial{NonCommutativeAlgebra}([3, 4]);

julia> sw = StateWord{Arbitrary}([m1]);

julia> ncsw = NCStateWord(sw, m2);

julia> ev = expval(ncsw);

julia> length(ev.state_syms)
2
```
"""
function expval(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer}
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
function Base.one(::Type{NCStateWord{ST,A,T}}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer}
    NCStateWord(one(StateWord{ST,A,T}), one(NormalMonomial{A,T}))
end

"""
    Base.one(ncsw::NCStateWord{ST,A,T}) where {ST,A,T}

Create the identity NCStateWord for the same type.
"""
function Base.one(::NCStateWord{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer}
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
produce a new canonical word (even when the result stays single-term for
`MonoidAlgebra` types).

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{UnipotentAlgebra}(UInt[1, 1]);  # x₁²

julia> sw = StateWord{Arbitrary}([one(NormalMonomial{UnipotentAlgebra,UInt})]);

julia> ncsw = NCStateWord(sw, m);

julia> result = simplify(ncsw);

julia> result isa NCStatePolynomial
true
```
"""
simplify(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer} = simplify!(copy(ncsw))
simplify!(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:MonoidAlgebra,T<:Integer} = ncsw

# =============================================================================
# Convenience Constructors: ς and tr
# =============================================================================

"""
    ς(m::NormalMonomial{A,T}) where {A,T}

Create a StateSymbol{Arbitrary} from a monomial.

This is a convenience function for creating state expectations in the
arbitrary state formalism. Equivalent to `StateSymbol{Arbitrary}(m)`.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

julia> sym = ς(m);

julia> sym isa StateSymbol{Arbitrary}
true
```
ς(m::NormalMonomial{A,T}) where {A<:MonoidAlgebra,T<:Integer} = StateSymbol{Arbitrary}(m)
ς(pairs::Vector{Tuple{Val{1},NormalMonomial{A,T}}}) where {A<:MonoidAlgebra,T<:Integer} = StateSymbol{Arbitrary}(pairs)

varsigma(args...) = ς(args...)
"""
ς(m::NormalMonomial{A,T}) where {A<:MonoidAlgebra,T<:Integer} = StateWord{Arbitrary}(m)

"""
    tr(m::NormalMonomial{A,T}) where {A,T}

Create a StateSymbol{MaxEntangled} from a monomial.

This is a convenience function for creating trace expressions in the
maximally entangled state formalism. Equivalent to `StateSymbol{MaxEntangled}(m)`.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

julia> sym = tr(m);

julia> sym isa StateSymbol{MaxEntangled}
true
```
tr(m::NormalMonomial{A,T}) where {A<:MonoidAlgebra,T<:Integer} = StateSymbol{MaxEntangled}(m)
tr(pairs::Vector{Tuple{Val{1},NormalMonomial{A,T}}}) where {A<:MonoidAlgebra,T<:Integer} = StateSymbol{MaxEntangled}(pairs)
"""
tr(m::NormalMonomial{A,T}) where {A<:MonoidAlgebra,T<:Integer} = StateWord{MaxEntangled}(m)

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
) where {A<:MonoidAlgebra, T<:Integer, ST<:StateType}

    # Collect monomials by *actual* degree after simplification.
    # For monoid algebras (e.g. unipotent), a word of length `deg` can simplify
    # to a lower-degree monomial, so bucketing by requested `deg` misses valid
    # StateWord/Monomial combinations.
    monos_by_deg = [NormalMonomial{A,T}[] for _ in 0:d]
    for elem in get_ncbasis(registry, d)
        mono = elem
        deg = degree(mono)
        deg <= d && push!(monos_by_deg[deg + 1], mono)
    end
    foreach(monos_by_deg) do monos
        unique!(sort!(monos))
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
    monos::Vector{NormalMonomial{A,T}},
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
