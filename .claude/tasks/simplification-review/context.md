# Task: Simplification Scheme Review

## Request
Create a comprehensive review plan for the simplification schemes across all algebra types in the NCTSSoS codebase. The user will execute this review in pair programming mode to understand the landscape before making changes.

## Project State
- **Branch**: fix/simplification-review
- **Codebase**: NCTSSoS - Polynomial optimization with sum-of-squares (SoS) methods
- **Core Module**: FastPolynomials (located in src/FastPolynomials/src/)
- **Recent Work**: Recent commits focus on fastpoly refactoring, basis generation, and canonicalization improvements

## Architecture Overview

### Algebra Type Hierarchy
Six algebra types are defined in `src/FastPolynomials/src/algebra_types.jl`:
1. **NonCommutativeAlgebra** - Generic non-commutative (no simplification rules)
2. **PauliAlgebra** - Pauli spin matrices (σᵢ² = I, cyclic products)
3. **FermionicAlgebra** - Fermionic operators (anticommutation, Wick's theorem)
4. **BosonicAlgebra** - Bosonic operators (commutation with delta terms)
5. **ProjectorAlgebra** - Projectors (Pᵢ² = Pᵢ, idempotent)
6. **UnipotentAlgebra** - Unipotent operators (U² = I, simpler than Pauli)

### Key Data Structures
- **Monomial{A,T}** - Single term with algebra type A and index type T
- **Polynomial{A,T,C}** - Sum of Term{Monomial{A,T}, C}
- **Variable** - Legacy compatibility struct (being phased out)
- **SimplifyAlgorithm** - Legacy config struct (compatibility layer in utils.jl)

## Simplification Landscape

### Simplification Files
Located in `src/FastPolynomials/src/simplification/`:
- `noncommutative.jl` (6.6k) - Site-based commutation, no algebraic rules
- `pauli.jl` (7.2k) - Cyclic products, site sorting, idempotency
- `fermionic.jl` (12k) - Generalized Wick's theorem, normal ordering
- `bosonic.jl` (16k) - Rook numbers on Ferrers boards, delta term generation
- `projector.jl` (5.2k) - Site-based idempotency (P² = P)
- `unipotent.jl` (7.9k) - Site-based pair removal (U² = I)

### Canonicalization
Located in `src/FastPolynomials/src/canonicalization.jl`:
- **symmetric_canon** - min(word, reverse(word)) for eigenvalue optimization
- **cyclic_canon** - Smallest cyclic rotation for trace optimization
- **cyclic_symmetric_canon** - Combined cyclic + symmetric
- **canonicalize** - Unified interface with `cyclic::Bool` parameter

### Common Patterns Across Algebras

#### Pattern 1: Return Types
- **Single Term**: NonCommutativeAlgebra, PauliAlgebra, ProjectorAlgebra, UnipotentAlgebra
  - `simplify(m) -> Term{Monomial,C}` where C is Float64 or ComplexF64
- **Multiple Terms**: FermionicAlgebra, BosonicAlgebra
  - `simplify(m) -> Vector{Term{Monomial,Float64}}` due to delta corrections

#### Pattern 2: Site-Based Encoding
Algebras using site-based commutation (Unsigned types only):
- NonCommutativeAlgebra, ProjectorAlgebra, UnipotentAlgebra
- Use bit-packed encoding: `index = (operator_id << site_bits) | site`
- Functions: `encode_index`, `decode_site`, `decode_operator_id` in algebra_types.jl
- Operators on different sites commute, same site preserves order

#### Pattern 3: Index Type Conventions
- **UInt16**: PauliAlgebra, ProjectorAlgebra, UnipotentAlgebra (self-adjoint operators)
- **Int32**: FermionicAlgebra, BosonicAlgebra (signed for creation/annihilation)
- **Flexible**: NonCommutativeAlgebra (supports both Signed and Unsigned)

#### Pattern 4: Multiplication Overloading
Each algebra implements `Base.:*(m1::Monomial{A}, m2::Monomial{A})`:
- Concatenates words: `vcat(m1.word, m2.word)`
- Calls `simplify!` on result
- Returns: Term or Vector{Term} depending on algebra

#### Pattern 5: Adjoint Operations
- Default in monomial.jl: `adjoint(m) = reverse(m.word)` with sign handling for signed types
- Specialized: UnipotentAlgebra (signed version) in unipotent.jl

## Known Issues and TODOs

### From Code Comments
1. **projector.jl:86** - "TODO: need to look at NCTSSOS this algorithm seems bad"
2. **pauli.jl:270** - "TODO: I would like to perform simplification during multiplication process itself"
3. **canonicalization.jl:143** - "TODO: do I need to simplify first before comparing?"
4. **bosonic.jl:53** - "TODO: Verify complexity bounds against literature"

### Inconsistencies Identified
1. **Return type heterogeneity**: Some algebras return single Term, others Vector{Term}
2. **Site-based encoding**: Only Unsigned types support it, Signed fallback is simple concatenation
3. **Legacy compatibility**: SimplifyAlgorithm in utils.jl provides compatibility but adds complexity
4. **Canonicalization integration**: Unclear when to canonicalize vs simplify

## Decisions
- Focus on MAPPING the current state, not proposing changes
- Document what exists, note patterns and inconsistencies
- Create actionable review items for pair programming sessions

## Progress
- [x] Explore codebase structure
- [x] Map all algebra types and their simplification schemes
- [x] Create context.md
- [x] Create plan.md
- [x] Create features.json

## Implementation Summary

### Task: Refactor Monomial Multiplication
Refactored all algebra-specific `*` (multiplication) functions for `Monomial` types to simply concatenate `word` fields instead of doing algebra-specific simplification.

### Files Modified
1. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/simplification/noncommutative.jl`
2. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/simplification/pauli.jl`
3. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/simplification/fermionic.jl`
4. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/simplification/bosonic.jl`
5. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/simplification/projector.jl`
6. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/simplification/unipotent.jl` (both Unsigned and Signed versions)
7. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/utils.jl` (neat_dot, _neat_dot3, Variable multiplication, Variable power)
8. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/polynomial.jl` (_add_simplified_terms! to handle Monomial returns)
9. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/state_word.jl` (NCStateWord multiplication)
10. `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/monomial.jl` (equality and hash handling for zero-hash monomials)

### Tests Updated
1. `/Users/yushengzhao/projects/NCTSSoS-main/test/fastpoly_test/monomials.jl`
2. `/Users/yushengzhao/projects/NCTSSoS-main/test/fastpoly_test/simplify.jl`
3. `/Users/yushengzhao/projects/NCTSSoS-main/test/fastpoly_test/arithmetic.jl`
4. `/Users/yushengzhao/projects/NCTSSoS-main/test/fastpoly_test/variables.jl`
5. `/Users/yushengzhao/projects/NCTSSoS-main/test/fastpoly_test/utils.jl`

### Key Changes
- **Monomial * Monomial** now returns `Monomial` (word concatenation only) instead of `Term`
- Callers must explicitly call `simplify!` if algebra-specific simplification is needed
- `_add_simplified_terms!` in polynomial.jl handles the new return type by calling `simplify!`
- `neat_dot` and `_neat_dot3` now directly use `adjoint(a) * b` which returns Monomial
- Fixed hash handling: equality and hash functions now handle zero-hash monomials correctly

### Test Status
- **FastPolynomials**: All 1060 tests pass
- Pre-existing failures remain in Aqua.jl, DocTest, Stale Imports, and Correlative Sparsity (solver-related)

## Handoff Summary
- All 6 algebra types identified with their simplification files
- Return type split: 4 algebras return Term, 2 return Vector{Term}
- Site-based encoding pattern documented for 3 algebras
- 4 TODOs found in code comments suggesting areas for improvement
- Legacy compatibility layer in utils.jl bridges old SimplifyAlgorithm API to new AlgebraType dispatch
- **COMPLETED**: Monomial multiplication refactored to return Monomial (word concatenation only)
