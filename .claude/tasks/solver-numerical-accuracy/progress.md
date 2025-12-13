# Progress Log: solver-numerical-accuracy

## Session: 2025-12-13 - Task Created

**Agent:** orchestrator
**Feature:** N/A (task initialization)

### Actions
- Created task from fastpoly-integration blocker
- Documented symptoms and likely root causes
- Identified files to investigate

### Outcome
Task initialized. Ready for investigation.

### Next Steps
- Create minimal reproduction case
- Trace moment matrix construction to identify where values diverge

---

## Session: 2025-12-13 - Bug Fixed

**Agent:** orchestrator

### Actions
1. Created minimal reproduction (N=2 Heisenberg, expected -0.75)
2. Found root cause: `_neat_dot3` returned unsimplified Monomials
3. Fixed `_neat_dot3` to return simplified Polynomials
4. Found second bug: `simplify!` returned monomials with stale hashes
5. Fixed all `simplify!` functions to create new Monomials with correct hashes
6. Verified fix: N=2 Heisenberg now returns -0.75 (error ~1e-10)
7. All 1357 tests passed

### Files Modified
- `src/FastPolynomials/src/utils.jl`: `_neat_dot3` now returns simplified Polynomial
- `src/FastPolynomials/src/simplification/pauli.jl`: Fixed hash staleness
- `src/FastPolynomials/src/simplification/noncommutative.jl`: Fixed hash staleness
- `src/FastPolynomials/src/simplification/projector.jl`: Fixed hash staleness
- `src/FastPolynomials/src/simplification/unipotent.jl`: Fixed hash staleness
- `src/moment_solver.jl`: Updated to handle Polynomial return from `_neat_dot3`
- `src/sparse.jl`: Updated to handle Polynomial return from `_neat_dot3`

### Root Causes
1. **Missing simplification**: `_neat_dot3(a, m, b)` computed `adjoint(a) * m * b` but didn't simplify, so Pauli rules (σ²=I) weren't applied in moment matrices
2. **Hash inconsistency**: `simplify!` modified monomial words in-place but returned original monomial with stale hash, causing `unique()` to fail deduplication

### Outcome
**FIXED** - Solver now returns correct numerical results:
- N=2 Heisenberg: -0.7500 (expected -0.75, error ~1e-10)

### Next Steps
- Clean up debug files
- Commit changes

---
