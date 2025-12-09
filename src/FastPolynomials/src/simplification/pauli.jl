"""
    Pauli Algebra Simplification

Implements simplification for Pauli spin operators satisfying:
- σᵢ² = I (involution: each Pauli squares to identity)
- Different sites commute: [σᵢⱼ, σₖₗ] = 0 for j ≠ l
- Same site cyclic products: σₓσᵧ = iσz, σᵧσz = iσₓ, σzσₓ = iσᵧ

# Variable Encoding Convention
Variables are ordered by site first, then by type (x, y, z):
- Index 1: σx₁, Index 2: σy₁, Index 3: σz₁
- Index 4: σx₂, Index 5: σy₂, Index 6: σz₂

For index `idx`:
- site = (idx - 1) ÷ 3 + 1
- pauli_type = (idx - 1) % 3 (0=X, 1=Y, 2=Z)

# Algorithm
1. Sort by site using stable sort (sites commute, stable sort gives deterministic ordering)
2. Linear pass: reduce each site group to at most one Pauli operator
3. Track complex phase coefficient through reductions

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);  # σx₁ σy₁ = iσz₁

julia> t = simplify(m);

julia> t.coefficient
0.0 + 1.0im

julia> t.monomial.word
1-element Vector{Int64}:
 3
```
"""

# Encoding helper functions (exported for testing)

"""
    _pauli_site(idx::Integer) -> Int

Extract site number from Pauli variable index.
Site = (idx - 1) ÷ 3 + 1
"""
@inline _pauli_site(idx::Integer) = (idx - 1) ÷ 3 + 1

"""
    _pauli_type(idx::Integer) -> Int

Extract Pauli type from variable index.
Type = (idx - 1) % 3, where 0=X, 1=Y, 2=Z
"""
@inline _pauli_type(idx::Integer) = (idx - 1) % 3

"""
    _pauli_index(site::Integer, type::Integer) -> Int

Create variable index from site and Pauli type.
Index = (site - 1) * 3 + type + 1
"""
@inline _pauli_index(site::Integer, type::Integer) = (site - 1) * 3 + type + 1

# Cyclic product table:
# σₐσᵦ = i * ε(a,b) * σ_c where c = (a+b) mod 3 if a≠b (sort of)
# More precisely:
#   XY = iZ, YZ = iX, ZX = iY (cyclic: phase = +i)
#   YX = -iZ, ZY = -iX, XZ = -iY (anti-cyclic: phase = -i)

"""
    _pauli_product(type1::Int, type2::Int) -> Tuple{ComplexF64, Int}

Compute the product of two Pauli operators on the SAME site.
Returns (phase, result_type) where result_type is 0, 1, or 2.

If type1 == type2, returns (1.0, -1) to indicate identity (σᵢ² = I).
"""
@inline function _pauli_product(type1::Int, type2::Int)
    # Same type: σᵢ² = I
    type1 == type2 && return (ComplexF64(1.0), -1)

    # Different types: cyclic or anti-cyclic product
    # XY→Z, YZ→X, ZX→Y (cyclic, +i)
    # YX→Z, ZY→X, XZ→Y (anti-cyclic, -i)

    # Result type: the one that's neither type1 nor type2
    # For {0,1,2}, the third one is 3 - type1 - type2
    result_type = 3 - type1 - type2

    # Determine phase: +i for cyclic order, -i for anti-cyclic
    # Cyclic order: 0→1→2→0 (X→Y→Z→X)
    # (type2 - type1 + 3) % 3 == 1 means cyclic
    phase = if (type2 - type1 + 3) % 3 == 1
        ComplexF64(0.0, 1.0)   # +i
    else
        ComplexF64(0.0, -1.0)  # -i
    end

    return (phase, result_type)
end

"""
    simplify!(m::Monomial{PauliAlgebra,T}) where T -> Term{Monomial{PauliAlgebra,T},ComplexF64}

In-place simplification of a Pauli algebra monomial.
Mutates the monomial's word vector and returns a Term.

Applies:
- Site sorting (operators on different sites commute)
- Cyclic product rules (σₓσᵧ = iσz, etc.)
- Idempotency (σᵢ² = I)

# Warning
This mutates the input monomial. Use `simplify` for a non-mutating version.

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2]);  # σx₁ σy₁

julia> t = simplify!(m);

julia> t.coefficient
0.0 + 1.0im

julia> t.monomial.word  # Result: σz₁
1-element Vector{Int64}:
 3

julia> m.word  # Original was mutated
1-element Vector{Int64}:
 3
```
"""
function simplify!(m::Monomial{PauliAlgebra,T}) where {T}
    word = m.word
    coef = ComplexF64(1.0)

    # Empty or single: nothing to simplify
    length(word) <= 1 && return Term(coef, m)

    # Stage 1: Sort by site (stable sort preserves within-site order for determinism)
    sort!(word, alg=InsertionSort, by=_pauli_site)

    # Stage 2: Reduce same-site groups in a single linear pass
    # After sorting, all operators on the same site are adjacent
    result = T[]
    i = 1
    while i <= length(word)
        current_site = _pauli_site(word[i])

        # Find extent of this site's group: word[i:j-1]
        j = i + 1
        while j <= length(word) && _pauli_site(word[j]) == current_site
            j += 1
        end

        # Reduce all operators in word[i:j-1] to at most one operator
        # current_type: -1 = identity accumulated, 0/1/2 = X/Y/Z accumulated
        current_type = _pauli_type(word[i])
        for k in i+1:j-1
            next_type = _pauli_type(word[k])
            if current_type == -1
                # Identity * next = next (no phase change)
                current_type = next_type
            else
                phase, result_type = _pauli_product(current_type, next_type)
                coef *= phase
                current_type = result_type  # Could be -1 (identity)
            end
        end

        # Append result for this site (skip if identity)
        if current_type != -1
            push!(result, T(_pauli_index(current_site, current_type)))
        end

        i = j  # Move to next site group
    end

    # Update word in-place
    resize!(word, length(result))
    copyto!(word, result)

    return Term(coef, m)
end

"""
    simplify(m::Monomial{PauliAlgebra,T}) where T -> Term{Monomial{PauliAlgebra,T},ComplexF64}

Simplify a Pauli algebra monomial.

Non-mutating version - creates a copy and simplifies it.

# Algebraic Rules
- σᵢ² = I (idempotency)
- Different sites commute
- Same site cyclic products: σₓσᵧ = iσz, σᵧσz = iσₓ, σzσₓ = iσᵧ

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2]);  # σx₁ σy₁

julia> t = simplify(m);

julia> t.coefficient
0.0 + 1.0im

julia> t.monomial.word
1-element Vector{Int64}:
 3

julia> m.word  # Original unchanged
2-element Vector{Int64}:
 1
 2
```
"""
function simplify(m::Monomial{PauliAlgebra,T}) where {T}
    # Copy and delegate to simplify!
    m_copy = Monomial{PauliAlgebra,T}(copy(m.word), m.hash)
    simplify!(m_copy)
end
