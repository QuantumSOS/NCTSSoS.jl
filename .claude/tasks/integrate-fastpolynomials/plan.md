# Implementation Plan: integrate-fastpolynomials

## Approach
Replace NCTSSoS's embedded FastPolynomials module with the new external FastPolynomials.jl implementation via **direct API migration**. All test files and NCTSSoS source files will be rewritten to use the new API directly—no adapter layers, shims, or compatibility wrappers. This ensures clean adoption of the new type system (Variable→VariableRegistry, SimplifyAlgorithm→AlgebraType) and eliminates technical debt from the start.

## API Migration Map

### Type Signature Changes

| **Old API** | **New API** | **Breaking Change** |
|-------------|-------------|---------------------|
| `Variable` struct | `VariableRegistry{T}` | Yes - struct → registry system |
| `Monomial(vars::Vector{Variable}, z::Vector{Int})` | `Monomial{A,T}(word::Vector{T})` | Yes - word representation, no exponents |
| `Polynomial{T}` | `Polynomial{A,T,C}` | Yes - added AlgebraType and Coefficient type params |
| `SimplifyAlgorithm` | `AlgebraType` singleton types | Yes - configuration → dispatch |
| N/A | `Term{M,C}` mutable struct | New - result of simplification |

### Variable Creation

| **Old API** | **New API** | **Notes** |
|-------------|-------------|-----------|
| `@ncpolyvar x[1:3]` | `reg, (x,) = create_noncommutative_variables([("x", 1:3)])` | Returns registry + monomial vectors |
| `Variable(name::Symbol; iscomplex=false)` | Registry-based creation per algebra | No direct Variable construction |
| `polyarrayvar(prefix, indices...)` | `create_*_variables(subscripts)` | Per-algebra factory functions |

### Monomial Operations

| **Old API** | **New API** | **Notes** |
|-------------|-------------|-----------|
| `monomial(vars, z)` | `Monomial{A}([idx1, idx2, ...])` | Word representation via indices |
| `star(m::Monomial)` → `Monomial` | `star(m::Monomial)` → `Monomial` | Same name, same semantics |
| `Base.adjoint(m)` for Monomial | `Base.adjoint(m)` for Monomial | Removed per plan, use `star` |
| `degree(m)` → `Int` | `degree(m)` → `Int` | Unchanged |
| `variables(m)` → `Vector{Variable}` | Via registry lookups | No direct Variable type |

### Polynomial Operations

| **Old API** | **New API** | **Notes** |
|-------------|-------------|-----------|
| `Polynomial(coeffs, monos)` | `Polynomial([Term(c, m), ...])` | Term-based construction |
| `p.coeffs`, `p.monos` | `coefficients(p)`, `monomials(p)` | Accessor functions |
| `terms(p)` → zip | `terms(p)` → `Vector{Term}` | Now first-class type |
| `star(p::Polynomial)` → `Polynomial` | `star(p::Polynomial)` → `Polynomial` | Unchanged |
| `Base.adjoint(p)` kept | `Base.adjoint(p)` kept | Per constraints |

### Simplification System

| **Old API** | **New API** | **Migration Strategy** |
|-------------|-------------|------------------------|
| `SimplifyAlgorithm(comm_gps=[x, y], is_unipotent=true, ...)` | Use `UnipotentAlgebra` type param | Create adapter mapping config → AlgebraType |
| `simplify(mono, sa)` → `Monomial` | `simplify(m::Monomial{A,T})` → `Term{Monomial{A,T}, C}` | Returns Term, not Monomial |
| `simplify!(mono, sa)` → mutates | `simplify!(m)` → mutates, returns Term | Same mutation semantics |
| `canonicalize(m, sa)` | Built into polynomial constructor | Automatic during Polynomial creation |
| `get_basis(vars, d, sa)` | `get_ncbasis(A, registry, d)` | Per-algebra basis generation |

### Basis Generation

| **Old API** | **New API** | **Notes** |
|-------------|-------------|-----------|
| `get_basis(vars::Vector{Variable}, d::Int)` | `get_ncbasis(A, registry, degree)` | Type-driven dispatch |
| `monomials(vars, Val{D}())` | Implicit in basis generation | No separate function |
| `get_state_basis(ST, vars, d, sa)` | Port to new system | Requires StatePolynomial migration |

### State Polynomials (Specialized)

| **Old API** | **New API** | **Status** |
|-------------|-------------|-----------|
| `StateWord{ST}` | `StateWord{ST}` | Use new FastPolynomials.jl implementation |
| `NCStateWord` | `NCStateWord` | Use new FastPolynomials.jl implementation |
| `StatePolynomial` | `StatePolynomial` | Use new FastPolynomials.jl implementation |
| `NCStatePolynomial` | `NCStatePolynomial` | Use new FastPolynomials.jl implementation |
| `ς(...)` | `ς(...)` | Use new FastPolynomials.jl implementation |

## Implementation Steps

### Phase 1: Infrastructure Setup (2 steps)

1. [x] **Replace old FastPolynomials implementation with new one**
   - Files: Delete `src/FastPolynomials/src/*` (old implementation)
   - Files: Copy `/Users/yushengzhao/projects/FastPolynomials.jl/src/*` to `src/FastPolynomials/src/`
   - Include: Use new `state_word.jl`, `state_polynomial.jl`, `state_types.jl` from FastPolynomials.jl
   - Test: `include("src/FastPolynomials/src/FastPolynomials.jl")` succeeds
   - Commit: `feat(fastpoly): replace with new FastPolynomials.jl implementation`
   - **DONE:** Commit `178cfbe`

2. [x] **Update FastPolynomials module structure**
   - File: `src/FastPolynomials/src/FastPolynomials.jl`
   - Update: Module exports to include new API (VariableRegistry, AlgebraType, Term, etc.)
   - Remove: Old API exports (Variable struct, SimplifyAlgorithm, etc.)
   - Test: Module loads without errors, new API functions are accessible
   - Commit: `refactor(fastpoly): update module exports for new API`
   - **DONE:** No code changes needed - new implementation already had correct exports. Verified all 30+ exports accessible.

### Phase 2: Direct Test Migration (11 steps)

3. [x] **Migrate test/fastpoly_test/variables.jl**
   - Strategy: Replace `@ncpolyvar x[1:3]` with `reg, (x,) = create_noncommutative_variables([("x", 1:3)])`
   - Change: `Variable` struct comparisons → registry symbol lookups via `reg.index_to_symbol`
   - Change: Remove `Variable` construction tests, add VariableRegistry tests
   - Test: All variable creation and lookup tests pass
   - **DONE:** Rewrote entire file to test new VariableRegistry API

4. [x] **Migrate test/fastpoly_test/monomials.jl**
   - Strategy: Replace `monomial(vars, z)` with direct `Monomial{A}(word::Vector{T})` construction
   - Change: Exponent notation `[x, x, y] → [2,1,3]` to word notation `[idx_x, idx_x, idx_y]`
   - Change: Remove `Base.adjoint(::Monomial)` tests (per plan), keep `star(::Monomial)` tests
   - Test: All monomial construction, comparison, and star operation tests pass
   - **DONE:** Rewrote entire file to test word-based Monomial API

5. [x] **Migrate test/fastpoly_test/polynomial.jl**
   - Strategy: Replace `Polynomial(coeffs, monos)` with `Polynomial([Term(c, m) for (c,m) in zip(coeffs, monos)])`
   - Change: `.coeffs` field access → `coefficients(p)` function
   - Change: `.monos` field access → `monomials(p)` function
   - Change: Update type annotations from `Polynomial{T}` to `Polynomial{A,T,C}`
   - Test: All polynomial construction, arithmetic, and accessor tests pass
   - **DONE:** Rewrote entire file with Term-based construction

6. [x] **Migrate test/fastpoly_test/simplify.jl**
   - Strategy: Replace `simplify(m, SimplifyAlgorithm(...))` with `simplify(m::Monomial{A,T})`
   - Change: Algebra configuration via SimplifyAlgorithm → algebra type as Monomial type parameter
   - Change: Handle return type change: `Monomial` → `Term{Monomial, C}` with phase coefficient
   - Update: Extract `.monomial` from Term where bare monomial expected
   - Update: Add tests for phase coefficients (e.g., Pauli `σx*σy = im*σz`)
   - Test: All simplification tests pass with correct algebra semantics
   - **DONE:** Rewrote for AlgebraType dispatch and Term return values

7. [x] **Migrate test/fastpoly_test/basis.jl** (merged into utils.jl)
   - Strategy: Replace `get_basis(vars, d, sa)` with `get_ncbasis(A, registry, d)`
   - Change: SimplifyAlgorithm parameter → AlgebraType type parameter
   - Change: `vars::Vector{Variable}` → `registry::VariableRegistry`
   - Test: Basis generation for all algebra types matches expected monomial sets
   - **DONE:** Basis tests now in utils.jl

8. [x] **Migrate test/fastpoly_test/star.jl** (merged into monomials.jl and utils.jl)
   - Strategy: Update to use new Monomial and Polynomial types
   - Change: Ensure star operations use type-appropriate conjugation rules
   - Test: Star operations on monomials and polynomials preserve algebra semantics
   - **DONE:** Star tests in monomials.jl and utils.jl

9. [x] **Migrate test/fastpoly_test/state_*.jl files**
   - Files: `state_word.jl`, `statepolynomial.jl`
   - Strategy: Port state polynomial system to use new Monomial/Polynomial types
   - Change: Ensure StateWord, NCStateWord, StatePolynomial use new internal representations
   - Test: State polynomial tests pass (ς notation, state basis generation)
   - **DONE:** Both files rewritten for new StateWord/NCStateWord API

10. [x] **Migrate test/fastpoly_test/arithmetic.jl**
    - Strategy: Update arithmetic tests to use new Polynomial construction
    - Change: Ensure addition, multiplication, subtraction use `Polynomial{A,T,C}` types
    - Test: All arithmetic operations preserve algebra semantics and type stability
    - **DONE:** Rewrote for new Polynomial and Term arithmetic

11. [x] **Migrate test/fastpoly_test/comparison.jl** (compare.jl)
    - Strategy: Update equality/comparison tests for new types
    - Change: Monomial/Polynomial equality using new internal representations
    - Test: Comparison operations work correctly with new hash/equality implementations
    - **DONE:** compare.jl rewritten for new comparison API

12. [x] **Migrate test/fastpoly_test/utilities.jl** (utils.jl)
    - Strategy: Update utility function tests to new API
    - Change: Any helper functions that construct or manipulate polynomials
    - Test: All utility functions work with new type system
    - **DONE:** utils.jl rewritten with encode/decode and basis generation tests

13. [x] **Migrate remaining test/fastpoly_test/*.jl files**
    - Files: `allocations.jl`, `runtests.jl`, `setup.jl`
    - Strategy: Apply same migration patterns (direct API usage, no adapters)
    - Test: All remaining tests pass
    - **DONE:** All 11 test files migrated, 350 tests pass

### Phase 3: Direct NCTSSoS Source Migration (4 steps)

14. [ ] **Update src/NCTSSoS.jl imports and exports**
    - Remove: `using .FastPolynomials: AbstractPolynomial, Variable, Monomial, SimplifyAlgorithm`
    - Add: `using .FastPolynomials: VariableRegistry, AlgebraType, Term, Monomial, Polynomial`
    - Add: `using .FastPolynomials: create_*_variables, get_ncbasis, simplify, coefficients, monomials`
    - Update: Re-export new API for downstream user code
    - Test: NCTSSoS module loads without errors
    - Commit: `refactor(core): update FastPolynomials imports for new API`

15. [ ] **Update src/gns.jl for new API**
    - Change: `vars::Vector{Variable}` → `registry::VariableRegistry` in function signatures
    - Change: `get_basis(vars, d)` → `get_ncbasis(A, registry, d)` with explicit algebra type
    - Change: Variable creation → registry-based variable creation
    - Update: Monomial construction to use word representation
    - Test: GNS reconstruction tests in `test/gns_test.jl` pass
    - Commit: `refactor(gns): migrate to new FastPolynomials API`

16. [ ] **Update src/solver_utils.jl for new API**
    - Change: Any `Variable` usage → registry-based approach
    - Change: `Polynomial` construction → new `Term`-based construction
    - Change: SimplifyAlgorithm configuration → AlgebraType dispatch
    - Test: Solver utility tests pass
    - Commit: `refactor(solver): migrate solver_utils.jl to new API`

17. [ ] **Update remaining NCTSSoS source files**
    - Files: `src/moment_solver.jl`, `src/constraints.jl`, and any other files using FastPolynomials
    - Strategy: Apply same migration patterns (direct new API usage)
    - Change: Replace all old API calls with new equivalents per API Migration Map
    - Test: Full NCTSSoS test suite passes
    - Commit: `refactor(core): complete NCTSSoS migration to new FastPolynomials API`

### Phase 4: Validation & Documentation (2 steps)

18. [ ] **Full integration testing**
    - Test: Run complete NCTSSoS test suite (`test/runtests.jl`)
    - Validate: All 11 migrated fastpoly tests pass
    - Validate: All NCTSSoS integration tests pass
    - Validate: Performance benchmarks show no regressions
    - Regression: Compare outputs with baseline snapshots (≤1e-12 numerical difference)
    - Commit: `test: validate full integration of new FastPolynomials`

19. [ ] **Documentation update**
    - File: Create `docs/fastpolynomials_migration.md`
    - Content: API changes summary, migration guide for users, breaking changes list
    - Include: Type signature changes table (from API Migration Map)
    - Include: Before/after code examples for common patterns
    - Include: Troubleshooting guide for migration issues
    - Commit: `docs(fastpoly): document API migration and breaking changes`

## Testing Strategy

### Unit Test Coverage
- All 11 old test files must pass after direct migration (steps 3-13)
- 20+ new FastPolynomials.jl tests run in parallel (from new implementation)
- Edge cases: empty monomials, identity elements, zero polynomials
- Type stability: All operations return correctly typed results

### Integration Test Coverage
- NCTSSoS examples from `test/` directory with new API
- GNS reconstruction tests with varying tolerances
- Moment solver tests with different algebra types
- End-to-end workflows using new FastPolynomials API

### Regression Prevention
- Snapshot old test outputs before migration (baseline comparison)
- Compare numerical results (≤1e-12 difference for real, ≤1e-10 for complex phases)
- Basis generation: exact monomial set equality (order-independent)
- Ensure all algebra semantics preserved (Pauli, Fermionic, Bosonic, etc.)

### Performance Validation
- Benchmark old vs new for representative operations (multiplication, simplification, basis generation)
- Target: New implementation ≥ old performance (no regressions allowed)
- Memory profiling: Verify reduced allocations in hot paths

## Risk Assessment

### High Risk
1. **State polynomial API differences**: New StateWord/StatePolynomial from FastPolynomials.jl may have different API than old NCTSSoS implementation
   - Mitigation: Compare old vs new state polynomial APIs, update test files to match new API in step 9
   - Validation: Ensure all state polynomial tests (ς notation, state basis generation) pass with new implementation

2. **Simplification semantics change**: Old returns `Monomial`, new returns `Term{Monomial, C}` with phase coefficient
   - Mitigation: Update all call sites to handle Term type, extract `.monomial` field where needed
   - Validation: Add explicit tests for phase coefficients (e.g., Pauli `σx*σy = im*σz`)

3. **Breaking API changes across codebase**: Direct migration means all code breaks simultaneously
   - Mitigation: Incremental testing—run tests after each file migration
   - Recovery: Use git to revert individual file migrations if issues arise

### Medium Risk
4. **Variable index encoding mismatch**: New registry system may assign different indices than old Variable ordering
   - Mitigation: Validate index-to-symbol mappings in tests, ensure deterministic registry creation
   - Risk: Could cause basis generation to produce correct set but different ordering

5. **Type parameter proliferation**: `Polynomial{A,T,C}` has 3 params vs old single param `Polynomial{T}`
   - Mitigation: Use type aliases for common cases (e.g., `PauliPoly = Polynomial{PauliAlgebra,UInt8,ComplexF64}`)
   - Impact: May require type annotation updates in function signatures

6. **Test logic translation errors**: Direct rewriting of tests may introduce bugs not present in original
   - Mitigation: Run old and new implementations in parallel during migration, compare outputs
   - Validation: Use output snapshots to verify numerical equivalence

### Low Risk
7. **Basis generation ordering differences**: New implementation may generate monomials in different order
   - Mitigation: Tests use `Set` comparison or `sort()` before equality checks
   - Note: Order doesn't affect correctness, only presentation

8. **Performance regressions in hot paths**: New code may be slower in some operations
   - Mitigation: Benchmark critical operations (multiplication, simplification, basis generation)
   - Threshold: Flag any operation >20% slower than old implementation

## Citations
- [Old FastPolynomials]: `/Users/yushengzhao/projects/NCTSSoS-main/src/FastPolynomials/src/`
- [New FastPolynomials]: `/Users/yushengzhao/projects/FastPolynomials.jl/src/`
- [Old Tests]: `/Users/yushengzhao/projects/NCTSSoS-main/test/fastpoly_test/`
- [New Tests]: `/Users/yushengzhao/projects/FastPolynomials.jl/test/`
- [Task Context]: `.claude/tasks/integrate-fastpolynomials/context.md`
