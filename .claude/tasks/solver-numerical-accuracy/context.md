# Task: solver-numerical-accuracy

## Request
Debug and fix numerical accuracy issues in `cs_nctssos()` solver after FastPolynomials migration.

## Background
The fastpoly-integration task refactored NCTSSoS to use the new FastPolynomials API with:
- Registry-based variable creation (`create_*_variables()`)
- Algebra type parameters (`PolyOpt{A,P}`, `CorrelativeSparsity{A,T,P,M}`)
- Unified `MomentProblem{A,T,M,P}` type

The refactoring compiles and passes unit tests (1357 tests), but **integration tests that run the solver produce incorrect numerical results**.

## Problem Statement

### Symptoms
The solver (`cs_nctssos`) returns wrong optimization values:

| Model | Expected | Got | Error |
|-------|----------|-----|-------|
| XXX Heisenberg (N=6) | -0.467129 | -0.480 | ~3% |
| J1-J2 Heisenberg (N=6) | -0.427 | -13.35 | ~3000% |
| Transverse Field Ising (N=3) | -1.017/-1.010 | -0.876 | ~14% |

### Key Observations
1. Solver **runs without errors** - no crashes or exceptions
2. SDP solver (Clarabel) converges normally
3. Results are consistently wrong across different physics models
4. The J1-J2 result is wildly off (order of magnitude), suggesting a fundamental issue

## Likely Root Causes

### 1. Moment Matrix Construction (moment_solver.jl)
- Basis monomials may be incorrectly mapped to JuMP variables
- `monomap::Dict{M, JuMP.Variable}` may have wrong entries
- Constraint matrices may be incorrectly constructed

### 2. Sparse.jl Algebra Handling
- `_neat_dot3()` may not correctly handle typed monomials
- `symmetric_canon()` may not work correctly with new types
- Term sparsity graph construction may be wrong

### 3. Index Type Mismatches
- Old API used `UInt64` indices, new uses algebra-specific types (UInt8, Int8)
- Basis generation (`get_ncbasis`) may have inconsistent index handling

### 4. Polynomial Simplification
- Pauli algebra simplification may not be triggered correctly
- Coefficients may be incorrect after simplification

## Files to Investigate

**Primary suspects:**
- `src/moment_solver.jl` - MomentProblem construction, monomap creation
- `src/sparse.jl` - Term sparsity, correlative sparsity
- `src/interface.jl` - cs_nctssos orchestration

**Supporting files:**
- `src/FastPolynomials/src/simplification/` - Algebra-specific simplification
- `src/FastPolynomials/src/basis.jl` - get_ncbasis implementation

## Debugging Strategy

1. **Compare old vs new API outputs:**
   - Run same problem with old API (if still available) and new API
   - Compare intermediate values: basis, monomap, constraint matrices

2. **Minimal reproduction:**
   - Create smallest possible test case that shows the bug
   - N=2 Pauli system might be sufficient

3. **Trace moment matrix construction:**
   - Print basis monomials
   - Print monomap dictionary
   - Verify constraint matrix structure

4. **Verify Pauli simplification:**
   - Check that σx² = I is correctly applied
   - Check that σx·σy = iσz produces correct coefficient

## Related Files

- `.claude/tasks/fastpoly-integration/` - Parent task (blocked by this issue)
- `test/algebra_constructors.jl` - Integration tests that fail (commented out)
- `test/heisenberg.jl` - Additional failing tests (disabled)

## Progress
- [ ] Create minimal reproduction case
- [ ] Trace moment matrix construction
- [ ] Identify root cause
- [ ] Implement fix
- [ ] Re-enable integration tests

## Decisions
(None yet)
