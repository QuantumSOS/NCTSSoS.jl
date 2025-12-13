# Task: solver-numerical-accuracy

## Request
Debug and fix numerical accuracy issues in `cs_nctssos()` solver after FastPolynomials migration.

## Background
The fastpoly-integration task refactored NCTSSoS to use the new FastPolynomials API with:
- Registry-based variable creation (`create_*_variables()`)
- Algebra type parameters (`PolyOpt{A,P}`, `CorrelativeSparsity{A,T,P,M}`)
- Unified `MomentProblem{A,T,M,P}` type

The refactoring compiles and passes unit tests (1357 tests), but **integration tests that run the solver produce incorrect numerical results**.

## Completion Summary

**Completed**: 2025-12-13

### Root Causes Identified & Fixed

**1. Missing Simplification in `_neat_dot3`**
- **Problem**: `_neat_dot3(a, m, b)` computed `adjoint(a) * m * b` but returned raw concatenated monomial without simplification
- **Impact**: Pauli algebra rules (σ²=I, σₓσᵧ=iσz) were not applied in moment matrix construction
- **Fix**: Changed `_neat_dot3` to call `simplify()` and return a `Polynomial` instead of `Monomial`

**2. Hash Staleness After Simplification**
- **Problem**: `simplify!()` modified monomial words in-place but returned Term with original monomial containing stale hash
- **Impact**: `unique()` in `sorted_union()` failed to deduplicate identical monomials, causing basis construction errors
- **Fix**: Modified all `simplify!()` functions to create new Monomial with freshly computed hash

### Files Modified

#### src/FastPolynomials/src/utils.jl
- `_neat_dot3(Monomial, Monomial, Monomial)`: Now returns `Polynomial` with simplified terms

#### src/FastPolynomials/src/simplification/*.jl
- `pauli.jl`: Create new Monomial with correct hash after simplification
- `noncommutative.jl`: Create new Monomial with correct hash after simplification
- `projector.jl`: Create new Monomial with correct hash after simplification
- `unipotent.jl`: Create new Monomial with correct hash after simplification

#### src/moment_solver.jl
- Updated `total_basis` computation to extract monomials from Polynomial return

#### src/sparse.jl
- Updated `get_term_sparsity_graph()` to handle Polynomial return from `_neat_dot3`
- Updated `term_sparsity_graph_supp()` to extract monomials from Polynomial results

### Testing Status

- [x] All 1357 unit tests passing
- [x] Minimal reproduction verified (N=2 Heisenberg: -0.75 expected, -0.7500 actual)
- [x] Error within numerical tolerance (~1e-10)

### Key Technical Insights

1. **Monomial hashing is critical**: Julia's `unique()` uses `hash()` for deduplication. Monomials with same word but different hash are treated as distinct.

2. **Simplification must return new objects**: Since `Monomial` has immutable `hash` field, `simplify!()` cannot truly be "in-place" - must create new Monomial with correct hash.

3. **Return type changes propagate**: Changing `_neat_dot3` return from `Monomial` to `Polynomial` required updates to all callers.

## Progress
- [x] Create minimal reproduction case
- [x] Trace moment matrix construction
- [x] Identify root cause
- [x] Implement fix
- [x] All tests passing

## Decisions
- Changed `_neat_dot3` return type from `Monomial` to `Polynomial` for correctness
- All `simplify!` functions now create new Monomial objects (semantic change from mutation)
