# =============================================================================
# Contract 1: word operations
# =============================================================================
#
# Provider: the core layer (packed words for each algebra).
# User:     the pipeline layer (block family, strengthening, symmetry reduction).
#
# The pipeline calls only the functions in this file. It has no knowledge of bit
# layouts. A change to this file needs the agreement of the architect thread.
# `test/next/contract_checks.jl` holds the conformance checks for a word type.

"""
    AbstractWord

A basis word of one algebra in canonical (normal) form.

A concrete word type `W` must satisfy all of these:

- `isbitstype(W)` is `true`. No operation in this contract allocates.
- `one(W)` is the empty word and `isone(w)` tests for it.
- `==`, `hash` and `isless` are defined. `isless` is a total order with the
  degree as its first key: `degree(a) < degree(b)` implies `a < b`.
- Each value of `W` is in normal form. A constructor that receives a word that is
  not in normal form must throw; it must not store it.
"""
abstract type AbstractWord end

"""
    Phase(k)

The scalar `im^k`, with `k` stored modulo 4. Products of words return a `Phase`
so that Pauli arithmetic is exact and has no floating-point coefficient.
"""
struct Phase
    k::UInt8
    Phase(k::Integer) = new(UInt8(mod(k, 4)))
end

Base.one(::Type{Phase}) = Phase(0)
Base.isone(p::Phase) = p.k == 0x00
Base.:*(a::Phase, b::Phase) = Phase(a.k + b.k)
Base.conj(p::Phase) = Phase(4 - p.k)
Base.isreal(p::Phase) = iseven(p.k)

"""
    scale(p::Phase, c::Number)

Return `im^k * c` in the number type of `c`. A real `c` with an odd `k` throws
`InexactError`: an algebra with imaginary product phases needs complex
coefficients.
"""
function scale(p::Phase, c::T) where {T<:Complex}
    p.k == 0x00 && return c
    p.k == 0x01 && return T(-imag(c), real(c))
    p.k == 0x02 && return -c
    return T(imag(c), -real(c))
end
function scale(p::Phase, c::T) where {T<:Real}
    p.k == 0x00 && return c
    p.k == 0x02 && return -c
    throw(InexactError(:scale, T, p))
end

# -----------------------------------------------------------------------------
# Inspection
# -----------------------------------------------------------------------------

"""
    degree(w::AbstractWord) -> Int

The number of letters of `w`. `degree(one(W)) == 0`.
"""
function degree end

"""
    sites(w::AbstractWord)

An iterator of the site (or mode) numbers on which `w` acts, as `Int`, in
ascending order and without repeats. It must not allocate.
"""
function sites end

# -----------------------------------------------------------------------------
# Adjoint and product
# -----------------------------------------------------------------------------

"""
    star(w::W) -> (Phase, W)

The adjoint of `w` as `phase * word`. In each algebra of this package the
adjoint of a normal-form word is one normal-form word with a phase.
`star` is an involution: `star` of the result gives `w` back, with the
conjugate phase.
"""
function star end

"""
    ProductStyle(::Type{W})

`SingleTerm()` when the product of two words of `W` is one word with a phase
(monoid algebras, Pauli, Majorana). `MultiTerm()` when the normal form of a
product is a sum of words (normal-ordered fermionic and bosonic words).
"""
abstract type ProductStyle end
struct SingleTerm <: ProductStyle end
struct MultiTerm <: ProductStyle end

"""
    mul(a::W, b::W) -> (Phase, W)

The product `a * b` in normal form. Required when `ProductStyle(W)` is
`SingleTerm()`.
"""
function mul end

"""
    push_term!(sink, coefficient::Number, w::AbstractWord)

Add `coefficient * w` to `sink`. The sink is owned by the caller: a polynomial
buffer, a moment accumulator, or a test dictionary. `mul!` and `act!` must not
allocate when `push_term!` does not.
"""
function push_term! end

"""
    mul!(sink, a::W, b::W, c::Number = 1) -> sink

Add the normal form of `c * a * b` to `sink` with one `push_term!` call for each
term. Coefficients have the number type of `c`. Required for each word type; the
method below covers `SingleTerm()` types.
"""
function mul!(sink, a::W, b::W, c::Number=1) where {W<:AbstractWord}
    return _mul!(ProductStyle(W), sink, a, b, c)
end
function _mul!(::SingleTerm, sink, a::W, b::W, c::Number) where {W<:AbstractWord}
    p, w = mul(a, b)
    push_term!(sink, scale(p, c), w)
    return sink
end

# -----------------------------------------------------------------------------
# Moments
# -----------------------------------------------------------------------------

"""
    moment_key(w::W) -> (key::W, phase::Phase, conjugated::Bool)

The representative of the moment of `w`. With `L` a state,
`L(w) == phase * (conjugated ? conj(L(key)) : L(key))`.

The method below uses only `L(star(w)) == conj(L(w))`: the key is the smaller of
`w` and the word of `star(w)`. A state type with more relations (for example a
trace) adds a method for its word type.
"""
function moment_key(w::W) where {W<:AbstractWord}
    p, s = star(w)
    # L(w) = conj(L(star(w))) = conj(p * L(s)) = conj(p) * conj(L(s))
    return isless(s, w) ? (s, conj(p), true) : (w, one(Phase), false)
end

"""
    moment_reality(key::W) -> Symbol

`:real` when `star(key)` is `+key`, so the moment of `key` is a real number;
`:imaginary` when `star(key)` is `-key`; `:complex` in the other cases.
The pipeline uses this to give a moment one real variable or two.
"""
function moment_reality(key::W) where {W<:AbstractWord}
    p, s = star(key)
    s == key || return :complex
    p.k == 0x00 && return :real
    p.k == 0x02 && return :imaginary
    throw(ArgumentError("star(w) == ±im * w is not possible for a word; got phase im^$(Int(p.k))."))
end

# -----------------------------------------------------------------------------
# Symmetry actions
# -----------------------------------------------------------------------------

"""
    AbstractAction

One symmetry generator: a map on words that keeps products and adjoints.

A concrete action type must define `act` (or `act!` when an image is a sum of words),
`compose(g, h)` (first `h`, then `g`), `inv`, `one(g)`, `isone`, `==` and
`hash`. The pipeline makes groups and orbits from these functions only.
"""
abstract type AbstractAction end

"""
    is_antiunitary(g::AbstractAction) -> Bool

`true` when `g` also conjugates scalars. The default is `false`.
"""
is_antiunitary(::AbstractAction) = false

"""
    act(g::AbstractAction, w::W) -> (Phase, W)

The image of `w` as `phase * word`, for the label, permutation and antiunitary
action types. The image is in normal form. For an antiunitary `g` this is the
image of the word; the caller conjugates the coefficient.
"""
function act end

"""
    act!(sink, g::AbstractAction, w::W, c::Number = 1) -> sink

Add `c * g(w)` to `sink`. Required for an action whose image of a word is a sum
of words. The method below covers each action with an `act` method. A continuous
symmetry is given by its generators, which are `AbstractDerivation`s.
"""
function act!(sink, g::AbstractAction, w::AbstractWord, c::Number=1)
    p, image = act(g, w)
    push_term!(sink, scale(p, c), image)
    return sink
end

"""
    compose(g::AbstractAction, h::AbstractAction)

The action `w -> g(h(w))`.
"""
function compose end

"""
    AbstractDerivation

One generator of a continuous symmetry group: a linear map `d` on polynomials
with `d(a * b) == d(a) * b + a * d(b)` and `d(star(a)) == star(d(a))`. The group
element is `exp(t * d)`. A derivation is not an `AbstractAction`: it does not
keep products, and it has no `compose` and no `inv`. A concrete type defines
`derive!`, `==` and `hash`.
"""
abstract type AbstractDerivation end

"""
    derive!(sink, d::AbstractDerivation, w::W, c::Number = 1) -> sink

Add `c * d(w)` to `sink` in normal form, with one `push_term!` call for each
term. It must not allocate when `push_term!` does not.
"""
function derive! end

# -----------------------------------------------------------------------------
# Width selection
# -----------------------------------------------------------------------------

"""
    word_type(algebra, nsites::Integer, max_letters::Integer) -> Type{<:AbstractWord}

The smallest word type of `algebra` that holds each word on `nsites` sites (or
modes) with at most `max_letters` letters. Throws an `ArgumentError` that names
the hard maximum when no type is large enough. The pipeline calls this one time
for a problem; `max_letters` includes the letters that strengthening adds.

`algebra` is a descriptor value of the core layer, not only a tag. It carries
the data that the width needs and that the two integers do not give, for example
the number of different letters on a site for a site-word algebra. The pipeline
gets the descriptor from the registry and passes it on; it does not read it.
"""
function word_type end

# -----------------------------------------------------------------------------
# Bridge to the current package. Temporary: removed with the old code.
# -----------------------------------------------------------------------------

"""
    legacy_letters(w::AbstractWord) -> Vector{Int}

The index vector that `NCTSSoS.NormalMonomial` uses for the same word. Define it
for a word type whose words map one to one to the normal-form words of the
current package. Used only by the tests that compare the new code with the
current package.
"""
function legacy_letters end

"""
    word_from_legacy(::Type{W}, letters::AbstractVector{<:Integer}) -> W

The inverse of `legacy_letters`. `letters` is a word in the normal form of the
current package.
"""
function word_from_legacy end

"""
    legacy_terms(w::AbstractWord) -> Vector{Tuple{ComplexF64,Vector{Int}}}

The expansion of `w` in normal-form words of the current package, as
`(coefficient, letters)` terms. The comparison tests use this function. The
method below covers each word type with `legacy_letters`. A word type whose
words are sums in the current normal form (for example Majorana words against
normal-ordered fermionic words) defines `legacy_terms` and
`terms_from_legacy!` and does not define the two functions above.
"""
legacy_terms(w::AbstractWord) = [(one(ComplexF64), legacy_letters(w))]

"""
    terms_from_legacy!(sink, ::Type{W}, letters::AbstractVector{<:Integer}, c::Number = 1) -> sink

Add the expansion of `c` times the legacy word `letters` in words of `W` to
`sink`, with `push_term!`. The inverse of `legacy_terms` on polynomials.
"""
function terms_from_legacy!(sink, ::Type{W}, letters::AbstractVector{<:Integer}, c::Number=1) where {W<:AbstractWord}
    push_term!(sink, c, word_from_legacy(W, letters))
    return sink
end
