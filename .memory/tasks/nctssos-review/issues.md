# Issues Found During Review

## Bugs

*(None found yet)*

---

## Performance Concerns

*(None found yet)*

---

## Maintainability Improvements

### MAINT-001: StatePolyOpt/PolyOpt Duplication
**Files**: `problem.jl`, `interface.jl`, `moment.jl`, `sparsity.jl`
**Impact**: ~200-300 lines duplicated per file
**Description**: `StatePolyOpt` and `PolyOpt` have nearly identical structs and parallel `cs_nctssos()` implementations.
**Root cause**: Extra `ST<:StateType` parameter + `NCStatePolynomial` vs `Polynomial`.
**Possible fix**: Unify into `PolyOpt{A, P<:AbstractPolynomial}` with dispatch on `P`.
**Trade-off**: Current = explicit types, easy dispatch. Unified = DRY but complex generics.

---

## Test Coverage Gaps

### TEST-001: compute_sparsity tests are structural only
**File**: `test/relaxations/interface.jl`
**Description**: Current tests only verify field existence and types, not concrete values.
**Action**: When debugging a failing polynomial optimization, add the problem as a concrete test case with expected `initial_activated_supps` and `cliques_term_sparsities` values.

---

## Design Questions

*(To be filled during review)*
