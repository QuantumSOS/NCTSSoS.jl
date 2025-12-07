# Progress Log: basis-registry

## Session: 2025-12-07 21:45

**Agent:** orchestrator
**Feature:** Registry-aware basis generation - test migration

### Actions
- Updated test/fastpoly_test/utils.jl to use registry-based API
- Updated test/fastpoly_test/simplify.jl - removed redundant basis tests
- Updated test/fastpoly_test/allocations.jl to use registry-based API
- Removed tests for deleted functionality (has_consecutive_repeats, filter_constraint)

### Outcome
All 1017 FastPolynomials tests passing.
Old API removed from basis.jl, fastpoly_test/ tests migrated.

### Current State
- Steps 1-3: COMPLETE (VariableRegistry{A,T}, registry-aware get_ncbasis)
- Step 4 (Remove old API): COMPLETE for fastpoly_test/
- Step 5 (Update tests): COMPLETE for fastpoly_test/
- External tests (sparse.jl, etc.): NOT UPDATED (per user request)

**Commit:** `319e10a` - refactor(fastpoly): migrate tests to registry-based basis API

---

## Session: 2025-12-07 18:30 - Checkpoint

**Agent:** orchestrator + polyglot-impl
**Feature:** Registry-aware basis generation

### Actions
- Validated original plan with lead-researcher (found algorithm issues)
- Revised plan: "generate all → simplify → return Terms" approach
- Implemented VariableRegistry{A,T} type parameter
- Added registry-aware get_ncbasis and get_ncbasis_deg returning Vector{Term}
- Added 52 new tests for registry-based API
- Started removing old API (in progress)

### Outcome
Core functionality complete. New registry-aware API working with all algebra types.
Old API removal in progress - basis.jl updated, tests partially updated.

### Current State
- Steps 1-3: COMPLETE (VariableRegistry{A,T}, registry-aware get_ncbasis)
- Step 4 (Remove old API): IN PROGRESS - 60%
  - basis.jl: old functions removed
  - test/fastpoly_test/basis.jl: updated
  - test/fastpoly_test/utils.jl: needs update
  - utils.jl get_basis wrappers: needs update
  - Other dependent files: needs review

### Blockers
- Many files use old get_ncbasis(Algebra, n, d) and get_basis APIs
- Full migration requires updating: sparse.jl, gns.jl, test files, docs

### Next Steps
1. Update test/fastpoly_test/utils.jl to use registry API
2. Update src/FastPolynomials/src/utils.jl (remove/update get_basis wrappers)
3. Update other dependent code (sparse.jl, gns.jl, etc.)
4. Run full test suite to verify

**Commit:** (WIP commit to follow)

---
