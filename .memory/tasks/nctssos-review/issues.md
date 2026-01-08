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

*(None found yet)*

---

## Design Questions

*(To be filled during review)*
