# StatePolyOpt Solver Support

## Task Status: ✅ COMPLETE

## Summary

StatePolyOpt solver support is fully implemented and working. All tests pass including the previously broken test 7.2.1 (squared expectations).

## Working Features

1. **StatePolyOpt type**: `polyopt(ncsp, registry)` creates state polynomial optimization problems
2. **Solver interface**: `cs_nctssos(spop, solver_config)` solves state polynomial optimization
3. **Basic state polynomials**: Objectives with single expectations like `<xy>` work correctly
4. **Covariance expressions**: Test 7.2.2 with `cov(a,b) = <xy> - <x><y>` works
5. **Squared expectations**: Test 7.2.1 with `<xy>²` terms now works correctly

## Test Status

- Test 7.2.0: PASS (-2.828...)
- Test 7.2.0 Sparse (with MMD term sparsity): SKIP (term sparsity needs separate fix)
- Test 7.2.1: PASS (-4.0) ✅ Fixed!
- Test 7.2.2: PASS (-5.0)

## The Fix: Complete Basis Generation

### Problem
The original `get_state_basis` only generated `<I>*M` forms (identity StateWord, operator monomial). This couldn't represent compound expectations like `<xy>²` needed for test 7.2.1.

### Solution
Modified `get_state_basis` to generate **all (StateWord, Monomial) combinations** where `degree(sw) + degree(nc_word) <= d`:

| Form | Example | Purpose |
|------|---------|---------|
| `<I>*M` | `<I>*xy` | Operator constraints (existing) |
| `<M>*I` | `<xy>*I` | Single expectation with identity operator |
| `<M1><M2>*I` | `<x><y>*I` | Compound expectations |
| `<M>*N` | `<x>*y` | Mixed forms (bridges both requirements) |

The key insight from analyzing the main branch was that **mixed forms** like `<M>*N` are essential for bridging the constraint structures needed by different test cases.

### Why Previous Attempts Failed

1. **Adding only `<M>*I` forms**: Caused double-counting because `<M>*I` and `<I>*M` have the same `expval` but produce different moment matrix entries.

2. **Using only `<M>*I` forms**: Broke test 7.2.0 which needs `<I>*M` forms for operator product constraints.

3. **The solution**: Generate ALL combinations including mixed forms. This provides a complete basis that can represent both operator constraints and compound expectations correctly.

## Implementation Details

### Files Modified

- `src/FastPolynomials/src/state_word.jl`: 
  - Rewrote `get_state_basis` to generate all (StateWord, Monomial) combinations
  - Added `_generate_statewords_up_to_degree` helper function
  - Removed `required_state_words` parameter (no longer needed)

- `src/sparse.jl`: 
  - Updated `correlative_sparsity` for StatePolyOpt to use simplified `get_state_basis`

- `test/state_poly_opt.jl`: 
  - Changed `@test_skip` to `@test` for test 7.2.1

### Known Issue: Term Sparsity (MMD) with State Polynomials

The MMD term sparsity algorithm creates incorrect blocks for state polynomial optimization.
When enabled, it returns 0 instead of the correct value. This is because the term sparsity
graph construction doesn't properly handle NCStateWord/StateWord types.

Workaround: Use `ts_algo=NoElimination()` for state polynomial optimization.

## Original Files Modified for StatePolyOpt Support

- `src/FastPolynomials/src/state_polynomial.jl`: Added `variable_indices`, `maxdegree`
- `src/FastPolynomials/src/state_word.jl`: Added `symmetric_canon` for StateWord/NCStateWord
- `src/FastPolynomials/src/utils.jl`: Extended `_neat_dot3` for NCStateWord
- `src/sparse.jl`: Added `StateCorrelativeSparsity`, `correlative_sparsity(::StatePolyOpt)`
- `src/moment_solver.jl`: Added `StateMomentProblem`, `moment_relax(::StatePolyOpt)`
- `src/sos_solver.jl`: Added `sos_dualize(::StateMomentProblem)`, `_get_state_Cαj`
- `src/interface.jl`: Added `StatePolyOptResult`, `cs_nctssos(::StatePolyOpt)`

See `.claude/tasks/state-poly-fix/plan.md` for detailed fix documentation.
