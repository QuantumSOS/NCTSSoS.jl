# Task: mutable-monomial

## Request
Make `Monomial` a mutable struct to allow `simplify!` to directly mutate `word` and update `hash`. Remove legacy code that is no longer needed.

## Project State
- NCTSSoS.jl is a Julia package for sparse noncommutative polynomial optimization
- Core types: `Monomial{A,T}` → `Term{M,C}` → `Polynomial{A,T,C}`
- Current Monomial is immutable but `simplify!` already mutates the word vector
- Six algebra types with different simplification patterns

## Research Findings

### Current Implementation
```julia
# src/FastPolynomials/src/monomial.jl (lines 78-81)
struct Monomial{A<:AlgebraType,T<:Integer} <: AbstractMonomial
    word::Vector{T}
    hash::UInt64
end
```

### Current simplify! Pattern
All `simplify!` functions follow this pattern:
```julia
function simplify!(m::Monomial{SomeAlgebra,T})
    word = m.word
    # ... mutate word in-place (sort!, deleteat!, etc.)
    # Create NEW Monomial because hash can't be updated
    return Monomial{SomeAlgebra}(word)
end
```

### Key Findings
1. Monomial IS NOT used as Dict/Set keys in most places (word vectors are used instead)
2. GNS code uses monomials as Dict keys but only after simplification (safe)
3. Term is already mutable - there's precedent
4. Three algebra patterns:
   - Simple (NC, Unipotent, Projector): Return same Monomial, just update hash
   - Pauli: Returns Term (builds new result vector)
   - Fermionic/Bosonic: Return Polynomial (multiple terms)

## Decisions
- Use `update_hash!(m)` helper function for DRY principle
- For simple algebras: mutate in-place and return same monomial
- For Pauli/Fermionic/Bosonic: keep return type semantics (Term/Polynomial)
- No legacy support code needed

## Progress
- [x] Research complete
- [x] Implementation
- [x] Testing (all 1425 tests pass)
- [x] Documentation update

## Handoff Summary
- **Completed**: Full implementation of mutable Monomial refactor
- **Key Finding**: Only simple algebras (NC, Unipotent, Projector) benefit from in-place mutation
- **Decision Made**: Add `update_hash!` helper, change struct to mutable
- **Git Commit**: `1c40eb4` - feat(FastPolynomials): make Monomial mutable for in-place simplification
- **Next Step**: Push to remote, consider merging to main
