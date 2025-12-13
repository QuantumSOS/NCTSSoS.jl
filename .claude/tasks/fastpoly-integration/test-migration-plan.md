# Test Migration Plan

## Summary
- **Total files**: 11
- **EASY (mechanical migration)**: 4 files
- **MEDIUM (needs attention)**: 4 files
- **HARD (comment out - state poly related)**: 2 files
- **NO CHANGES NEEDED**: 1 file

## Key API Changes Reference

### Variable Creation
- **OLD**: `@ncpolyvar x[1:N]`
- **NEW**: `registry, (x,) = create_noncommutative_variables([("x", 1:N)])`
- **NEW (Pauli)**: `registry, (x, y, z) = create_pauli_variables(1:N)`

### Basis Generation
- **OLD**: `get_basis([x, y, z], order)`
- **NEW**: `get_ncbasis(registry, order)` → returns `Vector{Polynomial}`
- **Extract monomials**: `extract_monomials_from_basis(basis_polys)`

### Variable Extraction
- **OLD**: `variables(mono)` → `Vector{Variable}` (legacy type)
- **NEW**: `variable_indices(mono)` → `Set{T}` (just indices)

### Complex Polynomial API
- **OLD**: `cpolyopt(...)` and `ComplexPolyOpt` as separate type
- **NEW**: Both still work! `cpolyopt` is now an alias to `polyopt`, `ComplexPolyOpt` is alias to `PolyOpt`
- **No changes needed** for existing `cpolyopt` calls

### Problem Construction
- **OLD**: Manual `comm_gps`, `is_unipotent`, `is_projective` flags
- **NEW**: Algebra type determined by variable creation function
- **For tests**: Can still use old constructor syntax as fallback (legacy API)

## Migration Order (Recommended)

1. **Phase 1 - EASY files** (can be done in parallel):
   - `test/solver_utils.jl`
   - `test/Aqua.jl`
   - `test/Doctest.jl`
   - `test/ExplicitImports.jl`

2. **Phase 2 - MEDIUM files** (requires more care):
   - `test/moment_solver.jl`
   - `test/heisenberg.jl`
   - `test/algebra_constructors.jl`
   - `test/interface.jl`

3. **Phase 3 - HARD files** (comment out for now):
   - `test/state_poly_opt.jl`
   - `test/trace_poly_opt.jl`

4. **Phase 4 - Dualization tests** (most complex):
   - `test/sos_solver.jl` (many tests already wrapped in `if false`)

---

## Detailed File Analysis

### test/solver_utils.jl
**Category:** EASY
**Changes:**
- Line 5: `using NCTSSoS.FastPolynomials: get_basis` → Remove this import (not used in active code)
- Line 8: `@ncpolyvar x y z` → `registry, (x, y, z) = create_noncommutative_variables([("x", 1:1), ("y", 1:1), ("z", 1:1)])`
- Line 24-76: Entire "Simplify" testset is already wrapped in `if false` (basis.jl type mismatch noted)

**State Polynomial Related:** No
**Recommendation:** Enable (uncomment from runtests.jl)
**Notes for Review:**
- Only active test is "VectorConstraint Dim" which uses JuMP, not polynomial API
- The @ncpolyvar on line 8 is unused in the active test
- Could potentially just delete the unused line 8 instead of migrating

---

### test/moment_solver.jl
**Category:** MEDIUM
**Changes:**
- Line 14: `using NCTSSoS.FastPolynomials: get_basis` → Keep but note it's for line 188 usage
- Line 20: `@ncpolyvar x[1:N] y[1:N] z[1:N]` → `registry, (x, y, z) = create_pauli_variables(1:N)`
- Line 28: `cpop = cpolyopt(...)` → **NO CHANGE** (cpolyopt still works as alias)
- Line 36: `variables(mono)` → `variable_indices(mono)`
  - **Context**: Used in `issubset(sort!(variables(mono)), clique)`
  - **New**: `issubset(variable_indices(mono), clique)` (no sort! needed, already Set)
- Line 54: Same pattern as line 20
- Line 68: Same pattern as line 36
- Line 90: `@ncpolyvar x[1:2] y[1:2]` → `registry, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])`
- Line 107: `@ncpolyvar x[1:n]` → `registry, (x,) = create_noncommutative_variables([("x", 1:n)])`
- Line 148: `@ncpolyvar pij[1:length(vec_idx2ij)]` → `registry, (pij,) = create_unipotent_variables([("pij", 1:length(vec_idx2ij))])`
- Line 182: `@ncpolyvar x y z` → `registry, (x, y, z) = create_noncommutative_variables([("x", 1:1), ("y", 1:1), ("z", 1:1)])`
- Line 188: `get_basis([x, y, z], 2)` → `get_ncbasis(registry, 2)` (returns Polynomial[], may need extraction)
  - **Context**: Used in `monomap = Dict(get_basis([x, y, z], 2) .=> jm)`
  - **Issue**: Old get_basis returns Vector{Monomial}, new returns Vector{Polynomial}
  - **Fix**: Use `extract_monomials_from_basis()` or convert polynomials to monomials
- Line 196: `@ncpolyvar x[1:n]` → Same as line 107
- Line 241: `@ncpolyvar x[1:2]` → `registry, (x,) = create_noncommutative_variables([("x", 1:2)])`
- Line 273: `@ncpolyvar x[1:n]` → Same as line 107

**State Polynomial Related:** No
**Recommendation:** Enable
**Notes for Review:**
- Most challenging part: Line 188 `get_basis` usage needs conversion to handle Polynomial[] instead of Monomial[]
- Consider extracting helper function for `get_basis` → `get_ncbasis` + extraction pattern
- All `cpolyopt` calls can stay as-is (backward compatible)
- The `variables(mono)` → `variable_indices(mono)` changes are straightforward

---

### test/heisenberg.jl
**Category:** EASY
**Changes:**
- Line 17: `@ncpolyvar x[1:N] y[1:N] z[1:N]` → `registry, (x, y, z) = create_pauli_variables(1:N)`
- Line 26: `pop = cpolyopt(...)` → **NO CHANGE** (still works)
- Line 46: Same as line 17
- Line 52: Same as line 26
- Line 67: Same as line 17
- Line 73: Same as line 26

**State Polynomial Related:** No
**Recommendation:** Enable
**Notes for Review:**
- Very straightforward migration
- All three testsets follow identical pattern
- `cpolyopt` calls are already correct

---

### test/algebra_constructors.jl
**Category:** MEDIUM
**Changes:**
- Line 82-98: Integration tests use new `pauli_algebra()` helper which already returns proper structure
- Line 82: `sys = pauli_algebra(2)` already creates x, y, z using new API
- Line 89: `pop = cpolyopt(ham, sys)` uses new algebra interface
- Line 109-164: Tests for `cpolyopt` with algebra interface
- Line 168-273: Integration tests all use `pauli_algebra()` helper

**State Polynomial Related:** No
**Recommendation:** Enable
**Notes for Review:**
- **This file is already using the new API!**
- The `pauli_algebra()` function returns a struct with `.variables` containing properly created variables
- All usage is through the new algebra interface
- May need to verify that `pauli_algebra()` function exists and works

---

### test/sos_solver.jl
**Category:** HARD (mostly disabled)
**Changes:**
- Line 5: `using NCTSSoS.FastPolynomials: get_basis` → needed for disabled tests
- Line 16-29: Active test "I_3322 inequality" (LOCAL_TESTING only)
  - Line 16: `@ncpolyvar x[1:3]` → `registry, (x,) = create_projector_variables([("x", 1:3)])`
  - Line 17: `@ncpolyvar y[1:3]` → append to same create call OR create separately and merge registries
  - Line 22: `polyopt(-f; comm_gps=[x, y], is_projective=true)` → Need registry-based call
- Lines 35-68: "CS TS Example" - wrapped in `if false` (basis.jl type mismatch)
- Lines 99-115: "Cαj complex" - wrapped in `if false` (basis.jl type mismatch)
- Lines 118-327: All dualization tests wrapped in `if false` (basis.jl type mismatch)

**State Polynomial Related:** No
**Recommendation:** Comment out active test (line 15-29) for now
**Notes for Review:**
- Only one active test (I_3322) and it's LOCAL_TESTING only
- All other tests already disabled due to basis.jl type mismatch
- The I_3322 test uses projective algebra which may not be fully supported yet
- Recommend wrapping entire file in `if haskey(ENV, "LOCAL_TESTING_DISABLED")`

---

### test/interface.jl
**Category:** MEDIUM
**Changes:**
- Line 15: `@ncpolyvar x[1:N] y[1:N] z[1:N]` → `registry, (x, y, z) = create_pauli_variables(1:N)`
- Line 24: `pop = cpolyopt(...)` → **NO CHANGE**
- Line 35-52: "Naive Example" - wrapped in `if false`
- Line 55-71: "Naive Example 2" - wrapped in `if false`
- Line 76: Same as line 15
- Line 82: Same as line 24
- Line 92: `@ncpolyvar x[1:3] y[1:3]` → `registry, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])`
- Line 97: `polyopt(-f; comm_gps=[x, y], is_projective=true)` → Need registry-based call
- Lines 110-162: "Majumdar Gosh Model" - wrapped in `if false`
- Lines 165-188: "Problem Creation Interface" - wrapped in `if false`
- Lines 191-223: "README Example Unconstrained" - wrapped in `if false`
- Lines 226-257: "README Example Constrained" - wrapped in `if false`

**State Polynomial Related:** No
**Recommendation:** Enable active tests (lines 13-31, 73-107)
**Notes for Review:**
- Two active tests: "1D Transverse Field Ising Model" and "Example"
- Both are LOCAL_TESTING only
- The "Example" test uses projective algebra (line 97) which needs registry-based call
- Most tests already disabled with TODO comments about basis.jl type mismatch

---

### test/state_poly_opt.jl
**Category:** HARD - State Polynomial
**Changes:**
- **This entire file uses state polynomial features**: `ς()`, `NCStateWord`, `expval`, `get_state_basis`
- Line 20: `expval`, `get_state_basis`, `NCStateWord` imported from FastPolynomials
- Line 25-51: Uses `ς(x[1] * y[1])` notation (state polynomial syntax)
- Line 78-132: All tests use state polynomial operators

**State Polynomial Related:** YES - ENTIRE FILE
**Recommendation:** Comment Out (keep wrapped in `if false`)
**Notes for Review:**
- State polynomial feature is not yet fully migrated/tested
- The `ς()` operator and `NCStateWord` are specialized features
- All tests are already wrapped in `if false` with TODO comments
- Should remain disabled until state polynomial integration is complete

---

### test/trace_poly_opt.jl
**Category:** HARD - State Polynomial (Trace variant)
**Changes:**
- **This entire file uses trace polynomial features**: `tr()`, similar to state polynomials
- Line 2: `using NCTSSoS.FastPolynomials:tr, Monomial`
- Line 17: `p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * one(Monomial)`
- All tests use `tr()` operator wrapped around polynomials

**State Polynomial Related:** YES - ENTIRE FILE (trace variant)
**Recommendation:** Comment Out (keep wrapped in `if false`)
**Notes for Review:**
- Trace polynomial is related to state polynomial features
- All tests already wrapped in `if false`
- The `tr()` operator is a specialized feature similar to `ς()`
- Should remain disabled until trace polynomial integration is complete

---

### test/Aqua.jl
**Category:** EASY
**Changes:**
- **NO CHANGES NEEDED**
- This is a code quality test (Aqua.jl testing)
- No polynomial API usage

**State Polynomial Related:** No
**Recommendation:** Enable (uncomment from runtests.jl)
**Notes for Review:**
- Should work immediately after uncommenting

---

### test/Doctest.jl
**Category:** EASY
**Changes:**
- **NO CHANGES NEEDED**
- Tests documentation examples
- Documentation should already be updated

**State Polynomial Related:** No
**Recommendation:** Enable (uncomment from runtests.jl)
**Notes for Review:**
- May reveal documentation that needs updating
- If failures occur, they indicate docs need fixes, not test code

---

### test/ExplicitImports.jl
**Category:** EASY
**Changes:**
- **NO CHANGES NEEDED**
- Tests import hygiene
- No polynomial API usage

**State Polynomial Related:** No
**Recommendation:** Enable (uncomment from runtests.jl)
**Notes for Review:**
- Should work immediately after uncommenting

---

## State Polynomial Tests to Comment Out

1. **test/state_poly_opt.jl** - Uses `ς()` operator, `NCStateWord`, state-specific features
2. **test/trace_poly_opt.jl** - Uses `tr()` operator, trace-specific features

**Explanation**: Both files rely on state/trace polynomial features that are:
- Not fully integrated with new FastPolynomials API
- Already wrapped in `if false` blocks
- Specialized features beyond core polynomial optimization
- Should be migrated in a separate focused effort after core API is stable

---

## Implementation Notes for Software Engineer

### Common Transformation Patterns

#### Pattern 1: Simple Variable Creation
```julia
# OLD
@ncpolyvar x[1:N] y[1:N] z[1:N]

# NEW (Pauli algebra - for physics tests)
registry, (x, y, z) = create_pauli_variables(1:N)

# NEW (Generic non-commutative)
registry, (x, y, z) = create_noncommutative_variables([("x", 1:N), ("y", 1:N), ("z", 1:N)])
```

#### Pattern 2: Variable Indices from Monomial
```julia
# OLD
issubset(sort!(variables(mono)), clique)

# NEW
issubset(variable_indices(mono), clique)  # No sort! needed, returns Set
```

#### Pattern 3: Basis Generation (TRICKY - type mismatch)
```julia
# OLD
basis_monomials = get_basis([x, y, z], order)  # Returns Vector{Monomial}
monomap = Dict(basis_monomials .=> jm_variables)

# NEW
basis_polys = get_ncbasis(registry, order)  # Returns Vector{Polynomial}
# Option 1: Extract monomials if basis is monomial-only
basis_monomials = extract_monomials_from_basis(basis_polys)
monomap = Dict(basis_monomials .=> jm_variables)

# Option 2: If basis has multi-term polynomials, need different approach
# (This is the "basis.jl type mismatch" issue mentioned in TODOs)
```

#### Pattern 4: cpolyopt calls (NO CHANGE NEEDED)
```julia
# OLD
cpop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

# NEW - SAME SYNTAX WORKS (backward compatible)
cpop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

# Or use new registry-based API
pop = polyopt(ham, registry; eq_constraints=eq_cons)
```

### Algebra Type Selection Guide

| Test Feature | Algebra Type | Constructor Function |
|--------------|--------------|---------------------|
| Pauli matrices (σx, σy, σz) | `PauliAlgebra` | `create_pauli_variables(1:N)` |
| Unipotent (U² = I) | `UnipotentAlgebra` | `create_unipotent_variables([("x", 1:N)])` |
| Projector (P² = P) | `ProjectorAlgebra` | `create_projector_variables([("P", 1:N)])` |
| Generic non-commutative | `NonCommutativeAlgebra` | `create_noncommutative_variables([("x", 1:N)])` |

### Test-Specific Notes

1. **moment_solver.jl Line 188 issue**:
   - Old `get_basis` returned `Vector{Monomial{...UInt64...}}`
   - New `get_ncbasis` returns `Vector{Polynomial{...Int64...}}`
   - Type parameter mismatch (UInt64 vs Int64) is the "basis.jl type mismatch" TODO
   - May need to extract monomials AND convert index types

2. **Projective algebra tests**:
   - Tests using `is_projective=true` flag in old API
   - New API should use `create_projector_variables()`
   - Verify projector algebra constraints are auto-applied

3. **Complex polynomial tests**:
   - `cpolyopt` still works as backward-compatible alias
   - `ComplexPolyOpt` type alias still exists
   - No migration needed for complex coefficient handling

### Recommended Implementation Sequence

1. **Start with Aqua/Doctest/ExplicitImports** - instant wins, no changes needed
2. **Do heisenberg.jl** - simple pattern, good test case for Pauli variables
3. **Do solver_utils.jl** - minimal active code, easy verification
4. **Do algebra_constructors.jl** - verify new API is working end-to-end
5. **Tackle moment_solver.jl** - requires solving the get_basis/get_ncbasis conversion
6. **Do interface.jl** - similar to moment_solver but fewer tests active
7. **Leave sos_solver.jl commented** - mostly disabled anyway
8. **Leave state/trace tests commented** - separate feature workstream

### Migration Verification Checklist

For each migrated test file:
- [ ] All `@ncpolyvar` replaced with appropriate `create_*_variables()`
- [ ] Registry variable captured and passed to `polyopt()` if using new API
- [ ] All `variables(mono)` replaced with `variable_indices(mono)`
- [ ] All `get_basis([vars...], d)` replaced with `get_ncbasis(registry, d)` + extraction if needed
- [ ] Tests pass locally with `make test`
- [ ] Tests pass in CI (no LOCAL_TESTING failures if applicable)
- [ ] No deprecation warnings in test output
