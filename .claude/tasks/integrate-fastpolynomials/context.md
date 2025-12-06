# Task: integrate-fastpolynomials

## Request
Refactor the FastPolynomials folder in NCTSSoS by integrating the new custom implementation from `/Users/yushengzhao/projects/FastPolynomials.jl`, addressing performance, memory, and API issues.

## Clarified Specification

### What it does
Replace the existing `src/FastPolynomials/` module in NCTSSoS with the new implementation featuring:
- Variable registry system instead of `Variable` struct
- New `Polynomial{M,C}` type signature (vs old `Polynomial{A,T,C}`)
- Consolidated `star`/`star!` API (removed `adjoint` for monomials)
- Support for ComposedMonomial (mixed algebra types)
- Better simplification algorithms for Pauli, Fermionic, Bosonic, Projector, and Unipotent algebras

### Source Locations
- **New implementation:** `/Users/yushengzhao/projects/FastPolynomials.jl/src/` (17 files)
- **New tests:** `/Users/yushengzhao/projects/FastPolynomials.jl/test/` (20+ files)
- **Old tests to migrate:** `/Users/yushengzhao/projects/NCTSSoS-main/test/fastpoly_test/` (11 files)
- **Old implementation:** `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/` (10 files)

### Acceptance Criteria
1. All new FastPolynomials.jl tests pass
2. Migrate old NCTSSoS fastpoly tests (11 files) to new API - logic must still pass
3. API changes propagated to NCTSSoS dependent code

### Constraints
- Term struct must remain mutable
- Must support both `Monomial{A,T}` and `ComposedMonomial`
- `Base.adjoint` kept for Polynomial (Julia ecosystem compatibility)

### Out of Scope
- Adding new features beyond what's in FastPolynomials.jl
- Changing NCTSSoS core algorithms

## Project State

### Current NCTSSoS Structure
- Main module: `src/NCTSSoS.jl`
- FastPolynomials submodule: `src/FastPolynomials/src/` with 10 files
- Test files: `test/fastpoly_test/` with 11 test files

### New FastPolynomials.jl Structure
- 17 source files including simplification algorithms
- 20+ comprehensive test files
- Detailed refactoring plan in `plan.md`

## Decisions
- **API Migration Strategy**: **Direct migration—no adapter layers, shims, or compatibility wrappers**. All test files and NCTSSoS source files will be rewritten to use new API directly.
- **Type Signature Change**: Accept breaking change from `Polynomial{T}` to `Polynomial{A,T,C}` (algebra-aware)
- **Simplification Returns**: New system returns `Term{Monomial, Coefficient}` instead of bare `Monomial` (captures phases)
- **Variable System**: Replace `Variable` struct with `VariableRegistry{T}` system (type-safe index management)
- **Migration Approach Rationale**: User explicitly rejected adapter/shim approach in favor of clean, direct adoption of new API to avoid technical debt

## Progress
- [x] Requirements clarified with user
- [x] Research complete
- [x] Implementation plan created
- [x] Phase 1: Infrastructure Setup (Steps 1-2)
- [x] Phase 2: Direct Test Migration (Steps 3-13)
- [ ] Phase 3: NCTSSoS Source Migration (Steps 14-17)
- [ ] Phase 4: Validation & Documentation (Steps 18-19)
- [ ] QA Review

## Research Summary (DO NOT RE-RESEARCH)

### Key API Differences Discovered

**Type System Changes:**
1. **Variable**: Old uses `Variable` struct with `name::Symbol` and `iscomplex::Bool`. New uses `VariableRegistry{T<:Integer}` with bidirectional Dict mappings (symbol ↔ index). Registry creation is algebra-specific via factory functions (`create_pauli_variables`, `create_fermionic_variables`, etc.).

2. **Monomial**: Old uses `Monomial(vars::Vector{Variable}, z::Vector{Int})` with exponent representation. New uses `Monomial{A<:AlgebraType, T<:Integer}` with word representation `word::Vector{T}` and precomputed hash. Algebra type is type parameter for zero-cost dispatch.

3. **Polynomial**: Old is `Polynomial{T}(coeffs::Vector{T}, monos::Vector{Monomial})`. New is `Polynomial{A,T,C}(terms::Vector{Term{Monomial{A,T}, C}})`. Added algebra type `A` and coefficient type `C` as parameters. Terms are sorted/deduplicated automatically.

4. **SimplifyAlgorithm**: Old uses configuration struct `SimplifyAlgorithm(comm_gps, is_unipotent, is_projective)`. New uses singleton AlgebraType dispatch (`PauliAlgebra`, `FermionicAlgebra`, `BosonicAlgebra`, `ProjectorAlgebra`, `UnipotentAlgebra`, `NonCommutativeAlgebra`).

5. **Term**: New type `Term{M<:AbstractMonomial, C<:Number}` (mutable) is result of simplification. Pairs coefficient with monomial. Does not exist in old implementation.

**Operation Changes:**
- `simplify(m::Monomial, sa::SimplifyAlgorithm)` → `simplify(m::Monomial{A,T})` returns `Term` with phase coefficient (e.g., Pauli products return complex phases like `iσz`)
- `star(m::Monomial)` unchanged semantically, but `Base.adjoint` removed for Monomial per plan (kept for Polynomial)
- `canonicalize(m, sa)` built into Polynomial constructor, not separate function
- `get_basis(vars, d, sa)` → `get_ncbasis(A, registry, degree)` with algebra-specific dispatch

**Variable Creation Examples:**
```julia
# Old API
@ncpolyvar x[1:3]  # Creates Variable array

# New API
reg, (x,) = create_noncommutative_variables([("x", 1:3)])
# Returns: (VariableRegistry{UInt8}, (Vector{Monomial{NonCommutativeAlgebra,UInt8}},))

reg, (σx, σy, σz) = create_pauli_variables(1:2)
# Returns: (VariableRegistry{UInt8}, (Vector{Monomial{PauliAlgebra,UInt8}}, ...))
```

**Simplification Example:**
```julia
# Old: Returns Monomial
simplify(σx₁ * σy₁, SimplifyAlgorithm(comm_gps=...)) → σz₁

# New: Returns Term with phase
simplify(σx₁ * σy₁) → Term(0.0 + 1.0im, σz₁::Monomial{PauliAlgebra})
```

### Files Using FastPolynomials in NCTSSoS
**Source files (2):**
- `src/NCTSSoS.jl`: Imports Variable, Monomial, utility functions
- `src/gns.jl`: Uses `Variable`, `get_basis`, `monomial` for GNS reconstruction

**Test files (23):**
- Core fastpoly tests (11): `test/fastpoly_test/*.jl` - MUST be migrated
- Integration tests (12): Various NCTSSoS tests importing FastPolynomials

### Migration Complexity Assessment
- **High complexity**: State polynomial system (StateWord, NCStateWord, StatePolynomial) - needs careful porting
- **Medium complexity**: Simplification system (SimplifyAlgorithm → AlgebraType requires adapter)
- **Low complexity**: Basic arithmetic (addition, multiplication largely compatible via adapters)

### Citations
- [Old Variable implementation]: `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/variables.jl`
- [New VariableRegistry]: `/Users/yushengzhao/projects/FastPolynomials.jl/src/variable_registry.jl`
- [Old Monomial]: `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/monomials.jl`
- [New Monomial]: `/Users/yushengzhao/projects/FastPolynomials.jl/src/monomial.jl`
- [Old Polynomial]: `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/polynomial.jl`
- [New Polynomial]: `/Users/yushengzhao/projects/FastPolynomials.jl/src/polynomial.jl`
- [Old SimplifyAlgorithm]: `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/simplify.jl`
- [New Pauli simplification]: `/Users/yushengzhao/projects/FastPolynomials.jl/src/simplification/pauli.jl`
- [AlgebraTypes]: `/Users/yushengzhao/projects/FastPolynomials.jl/src/algebra_types.jl`
- [Term struct]: `/Users/yushengzhao/projects/FastPolynomials.jl/src/term.jl`

## Implementation Summary

### Phase 1 Complete (2025-12-06)

**Files deleted** (old implementation):
- `src/FastPolynomials/src/`: arithmetic.jl, compare.jl, monomials.jl, polynomial.jl, simplify.jl, state_word.jl, statepolynomial.jl, utils.jl, variables.jl, FastPolynomials.jl (10 files)

**Files created** (new implementation):
- `src/FastPolynomials/src/`: FastPolynomials.jl, algebra_types.jl, variable_registry.jl, monomial.jl, term.jl, polynomial.jl, composed_monomial.jl, canonicalization.jl, basis.jl, state_types.jl, state_word.jl, state_polynomial.jl (12 files)
- `src/FastPolynomials/src/simplification/`: pauli.jl, fermionic.jl, bosonic.jl, projector.jl, unipotent.jl, noncommutative.jl (6 files)

**Test status**: Module loads and exports verified (30+ symbols accessible)

**Commits**:
- `178cfbe` - feat(fastpoly): replace with new FastPolynomials.jl implementation

### Phase 2 Complete (2025-12-06)

**Files modified** (test migration):
- `test/fastpoly_test/setup.jl` - NEW: Loads FastPolynomials directly for testing
- `test/fastpoly_test/runtests.jl` - Updated to use setup.jl
- `test/fastpoly_test/variables.jl` - Rewritten for VariableRegistry API
- `test/fastpoly_test/monomials.jl` - Rewritten for word-based Monomial{A,T}
- `test/fastpoly_test/polynomial.jl` - Rewritten for Term-based construction
- `test/fastpoly_test/arithmetic.jl` - Updated for new Polynomial/Term arithmetic
- `test/fastpoly_test/compare.jl` - Updated for new comparison API
- `test/fastpoly_test/simplify.jl` - Rewritten for AlgebraType dispatch
- `test/fastpoly_test/state_word.jl` - Rewritten for new StateWord/NCStateWord
- `test/fastpoly_test/statepolynomial.jl` - Rewritten for new StatePolynomial
- `test/fastpoly_test/utils.jl` - Updated for encode_index, decode_*, basis generation
- `test/fastpoly_test/allocations.jl` - Updated allocation tests

**Test status**: 350 tests pass in fastpoly_test suite

## Handoff Summary (for successor agent)

**Completed**: Phase 2 Direct Test Migration (Steps 3-13) - all 11 test files migrated

**Key Actions**:
1. Created `setup.jl` to load FastPolynomials directly (bypasses NCTSSoS module during test migration)
2. Migrated all test files to new API:
   - `variables.jl` - VariableRegistry tests, create_*_variables functions
   - `monomials.jl` - Word-based Monomial{A,T} construction, star operation
   - `polynomial.jl` - Term-based construction, accessor functions
   - `arithmetic.jl` - Polynomial/Term arithmetic operations
   - `compare.jl` - Equality, hashing, ordering comparisons
   - `simplify.jl` - AlgebraType dispatch, Term return values
   - `state_word.jl` - StateWord, NCStateWord with new Monomial types
   - `statepolynomial.jl` - StatePolynomial, NCStatePolynomial operations
   - `utils.jl` - encode_index, decode_*, basis generation, has_consecutive_repeats
   - `allocations.jl` - Performance/allocation tests
   - `runtests.jl` - Test runner using setup.jl
3. Updated imports to use `.FastPolynomials` (relative to setup.jl include)
4. Fixed API naming differences (decode_operator -> decode_operator_id, etc.)
5. Fixed integer type constraints (UInt8 max sites = 3)
6. All 350 tests pass

**Decision Made**: Use setup.jl to load FastPolynomials directly for testing (NCTSSoS module migration is Phase 3)

**Blocker/Question**: None - Phase 2 complete

**Next Step**: Phase 3, Step 14 - Update src/NCTSSoS.jl imports and exports for new FastPolynomials API

**Git Commit**: `68638f2` - test(fastpoly): migrate all test files to new FastPolynomials API
