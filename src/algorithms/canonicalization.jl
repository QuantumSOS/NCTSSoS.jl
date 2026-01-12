canonicalize(::Type{Arbitrary},m::NormalMonomial) = symmetric_canon(m)
canonicalize(::Type{MaxEntangled},m::NormalMonomial) = cyclic_symmetric_canon(m)

"""
    symmetric_canon(m::NormalMonomial{A}) where {A<:AlgebraType} -> Vector

Return a *word-level* canonical representative for `m` under the involution
symmetry `w ~ reverse(w)` (equivalently, `m ~ adjoint(m)` for real monoid
algebras).

This is used to identify monomials that should be treated as equivalent in
moment/SOS constructions where expectations are real-valued.

Notes
- The return value is a fresh copy of the word; mutating it will not affect `m`.
- For `PauliAlgebra`, `FermionicAlgebra`, and `BosonicAlgebra`, this currently
  defaults to the identity representative (`copy(m.word)`).
"""
symmetric_canon(m::NormalMonomial{PauliAlgebra}) = copy(m.word)
symmetric_canon(m::NormalMonomial{FermionicAlgebra}) = copy(m.word)
symmetric_canon(m::NormalMonomial{BosonicAlgebra}) = copy(m.word)

function symmetric_canon(m::NormalMonomial{A}) where {A<:MonoidAlgebra}
    degree(m) <= 1 && return copy(m.word)

    word = copy(m.word)
    word_rev = simplify!(A, reverse(word))
    word_rev < word && return word_rev
    return word
end


"""
    cyclic_canon(m::NormalMonomial{A}) where {A<:AlgebraType} -> Vector

Return a *word-level* canonical representative for `m` under cyclic rotation of
the word (trace invariance: `tr(ABC) = tr(BCA) = tr(CAB)`).

Notes
- The return value is a fresh copy of the word.
- For `PauliAlgebra`, `FermionicAlgebra`, and `BosonicAlgebra`, this currently
  defaults to the identity representative (`copy(m.word)`).
"""

cyclic_canon(m::NormalMonomial{PauliAlgebra}) = copy(m.word)
cyclic_canon(m::NormalMonomial{FermionicAlgebra}) = copy(m.word)
cyclic_canon(m::NormalMonomial{BosonicAlgebra}) = copy(m.word)

function cyclic_canon(m::NormalMonomial{A}) where {A<:MonoidAlgebra}
    word = copy(m.word)
    n = length(word)
    n <= 1 && return word

    best = copy(word)
    for offset in 1:(n-1)
        rotation = [word[mod1(i + offset, n)] for i in 1:n]
        simplify!(A,rotation)
        if rotation < best
            best = rotation
        end
    end
    return best
end


"""
    cyclic_symmetric_canon(m::NormalMonomial{A}) where {A<:AlgebraType} -> Vector

Return a *word-level* canonical representative for `m` under both cyclic
rotation and involution (`w ~ reverse(w)`), implemented as:

`min(cyclic_canon(m), cyclic_canon(adjoint(m)))`.

Notes
- The return value is a fresh copy of the word.
- For `PauliAlgebra`, `FermionicAlgebra`, and `BosonicAlgebra`, this currently
  defaults to the identity representative (`copy(m.word)`).
"""
cyclic_symmetric_canon(m::NormalMonomial{PauliAlgebra}) = copy(m.word)
cyclic_symmetric_canon(m::NormalMonomial{FermionicAlgebra}) = copy(m.word)
cyclic_symmetric_canon(m::NormalMonomial{BosonicAlgebra}) = copy(m.word)

function cyclic_symmetric_canon(m::NormalMonomial{A}) where {A<:MonoidAlgebra}
    return min(cyclic_canon(m), cyclic_canon(adjoint(m)))
end
