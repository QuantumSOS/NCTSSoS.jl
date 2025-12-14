# Implementation Plan: Mutable Monomial

## Approach
Change `Monomial` from immutable `struct` to `mutable struct`, add `update_hash!` helper, and update simplify! functions in simple algebras to mutate in-place.

## Steps

### 1. [DONE] ✓ Update monomial.jl - Core struct definition
**File**: `src/FastPolynomials/src/monomial.jl`

Change line 78-81:
```julia
# BEFORE
struct Monomial{A<:AlgebraType,T<:Integer} <: AbstractMonomial
    word::Vector{T}
    hash::UInt64
end

# AFTER
mutable struct Monomial{A<:AlgebraType,T<:Integer} <: AbstractMonomial
    word::Vector{T}
    hash::UInt64
end
```

Add `update_hash!` helper after the struct definition:
```julia
"""
    update_hash!(m::Monomial{A,T}) where {A,T} -> Monomial{A,T}

Recompute the hash field from the current word vector.
Call this after mutating the word to maintain hash consistency.
"""
@inline function update_hash!(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    m.hash = hash(m.word)
    return m
end
```

### 2. [DONE] ✓ Update noncommutative.jl - simplify!
**File**: `src/FastPolynomials/src/simplification/noncommutative.jl`

```julia
# BEFORE (creates new Monomial)
function simplify!(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    word = m.word
    length(word) <= 1 && return m
    sort!(word, alg=Base.Sort.InsertionSort, by=decode_site)
    simplified_mono = Monomial{NonCommutativeAlgebra}(word)
    return simplified_mono
end

# AFTER (mutates in-place)
function simplify!(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    word = m.word
    length(word) <= 1 && return m
    sort!(word, alg=Base.Sort.InsertionSort, by=decode_site)
    update_hash!(m)
    return m
end
```

### 3. [DONE] ✓ Update unipotent.jl - simplify!
**File**: `src/FastPolynomials/src/simplification/unipotent.jl`

Replace the final `Monomial{UnipotentAlgebra}(word)` construction with:
```julia
update_hash!(m)
return m
```

### 4. [DONE] ✓ Update projector.jl - simplify!
**File**: `src/FastPolynomials/src/simplification/projector.jl`

Replace the final `Monomial{ProjectorAlgebra}(word)` construction with:
```julia
update_hash!(m)
return m
```

### 5. [DONE] ✓ Clean up pauli.jl
**File**: `src/FastPolynomials/src/simplification/pauli.jl`

Remove outdated comment about "can't update hash in-place since Monomial is immutable".
No functional change needed - Pauli builds a NEW result vector and returns Term.

### 6. [SKIP] Clean up fermionic.jl (no outdated comments found)
**File**: `src/FastPolynomials/src/simplification/fermionic.jl`

Remove any outdated comments about immutability.
No functional change - returns Polynomial.

### 7. [SKIP] Clean up bosonic.jl (no outdated comments found)
**File**: `src/FastPolynomials/src/simplification/bosonic.jl`

Remove any outdated comments about immutability.
No functional change - returns Polynomial.

### 8. [DONE] ✓ Add tests for hash consistency
**File**: `test/fastpoly_test/simplify.jl` or new file

```julia
@testset "Mutable Monomial Hash Consistency" begin
    # NonCommutative: verify in-place mutation
    m_nc = Monomial{NonCommutativeAlgebra}(UInt16[2, 1])
    original_id = objectid(m_nc)
    result_nc = simplify!(m_nc)
    @test objectid(result_nc) == original_id  # Same object
    @test m_nc.hash == hash(m_nc.word)  # Hash correct

    # Unipotent: with pair cancellation
    m_uni = Monomial{UnipotentAlgebra}(UInt16[1, 1])
    result_uni = simplify!(m_uni)
    @test result_uni === m_uni
    @test m_uni.hash == hash(m_uni.word)
    @test isempty(m_uni.word)  # U^2 = I = empty

    # Projector: with duplicate removal
    m_proj = Monomial{ProjectorAlgebra}(UInt16[1, 1, 1])
    result_proj = simplify!(m_proj)
    @test result_proj === m_proj
    @test m_proj.hash == hash(m_proj.word)
    @test length(m_proj.word) == 1  # P^3 = P
end
```

### 9. [DONE] ✓ Run test suite
```bash
make test-FastPoly
LOCAL_TESTING=true make test
```

## Testing Strategy
1. Run existing tests to ensure no regressions
2. Add specific tests for hash consistency after mutation
3. Verify `objectid` is same before/after simplify! for simple algebras
4. Test edge cases: empty word, single element, pair cancellation

## Code to Remove
- Comments mentioning "can't update hash in-place since Monomial is immutable"
- Any workaround code for immutability

## Citations
- Research from Explore agent (agentId: acff63d)
- Plan from Plan agent (agentId: ac0936e)
