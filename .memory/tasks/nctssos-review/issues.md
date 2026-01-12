# Issues Found During Review

## Bugs

### BUG-001: State correlative graph used per-monomial constraint handling [FIXED]
**File**: `src/optimization/sparsity.jl:661-666`
**Impact**: State polynomial constraints with variables in different monomials could become global constraints unnecessarily, reducing sparsity exploitation.
**Fix**: Changed to use `variable_indices(con)` for whole constraint (commit `5dcbae1`)
**Status**: ✅ FIXED

---

## Performance Concerns

### PERF-001: Wasteful Polynomial wrapping in sparsity.jl
**File**: `src/optimization/sparsity.jl` (multiple locations)
**Pattern**: `monomials(Polynomial(simplify(mono)))` wraps Monomial in Polynomial just to extract monomials
**Impact**: Allocates intermediate Polynomial struct unnecessarily (especially for SimpleAlgebra where simplify returns single Monomial)
**Fix**: See MAINT-002 (`simplified_monomials()` helper)

---

## Maintainability Improvements

### MAINT-001: StatePolyOpt/PolyOpt Duplication
**Files**: `problem.jl`, `interface.jl`, `moment.jl`, `sparsity.jl`
**Impact**: ~200-300 lines duplicated per file
**Description**: `StatePolyOpt` and `PolyOpt` have nearly identical structs and parallel `cs_nctssos()` implementations.
**Root cause**: Extra `ST<:StateType` parameter + `NCStatePolynomial` vs `Polynomial`.
**Possible fix**: Unify into `PolyOpt{A, P<:AbstractPolynomial}` with dispatch on `P`.
**Trade-off**: Current = explicit types, easy dispatch. Unified = DRY but complex generics.

### MAINT-002: Add `simplified_monomials()` helper
**Files**: `src/simplification/interface.jl` (new) or `src/util/helpers.jl`
**Description**: Create dispatch-based helper to extract monomials from simplification result without Polynomial wrapping.
```julia
# SimpleAlgebra: simplify → Monomial directly
simplified_monomials(m::Monomial{A,T}) where {A<:SimpleAlgebra,T} = [simplify(m)]

# PauliAlgebra: simplify → PauliMonomial, extract .mono
simplified_monomials(m::Monomial{PauliAlgebra,T}) where {T} = [simplify(m).mono]

# Fermionic/Bosonic: simplify → PhysicsMonomial, extract .monos
function simplified_monomials(m::Monomial{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T}
    pm = simplify(m)
    iszero(pm) && return Monomial{A,T}[]
    pm.monos
end
```
**Impact**: Removes ~6 `Polynomial(simplify(...))` calls in sparsity.jl, cleaner code.

### MAINT-003: Unify Polynomial/State sparsity functions
**Files**: `src/optimization/sparsity.jl`
**Description**: ~400 lines duplicated between Polynomial and NCStatePolynomial versions:
- `CorrelativeSparsity` / `StateCorrelativeSparsity` structs
- `get_correlative_graph` / `get_state_correlative_graph`
- `assign_constraint` / `assign_state_constraint`
- `correlative_sparsity` (two dispatches)
- `init_activated_supp`, `term_sparsities`, etc.

**Approach**:
1. Single `CorrelativeSparsity{A,T,P,M}` struct (no separate State version)
2. Dispatch helpers: `_get_basis(reg, order, ::Type{<:Polynomial})` vs `::Type{<:NCStatePolynomial}`
3. Common `_monomial_var_indices()` dispatch for Monomial vs NCStateWord

**Estimated reduction**: ~400 lines → ~200 lines

---

## Review Notes (Need Investigation)

### REVIEW-001: abs(idx) in get_positions treats a₁† and a₁ as same variable
**File**: `src/optimization/sparsity.jl:115`
**Code**: `abs_idx = T(abs(idx))`
**Question**: For Fermionic/Bosonic, creation (idx=1) and annihilation (idx=-1) are distinct operators. Is it correct to treat them as the "same variable" for sparsity graph construction?
**Intent**: Likely correct - they refer to the same mode/site, so should be in same clique.
**Action**: Verify against NCTSSOS Python implementation or physics test cases.

### REVIEW-002: monosquare asymmetry between Polynomial and State init_activated_supp
**File**: `src/optimization/sparsity.jl`
**Polynomial (line 466-476)**: Includes diagonal entries `b†b` in initial activated support
**State (line 866-874)**: Excludes diagonals (matches NCTSSOS `monosquare=false`)
**Comment in code**: "including all diagonal entries would create spurious edges"
**Question**: Should Polynomial also have `monosquare=false` option? Or is the difference intentional?

---

## Test Coverage Gaps

### TEST-001: compute_sparsity tests are structural only
**File**: `test/relaxations/interface.jl`
**Description**: Current tests only verify field existence and types, not concrete values.
**Action**: When debugging a failing polynomial optimization, add the problem as a concrete test case with expected `initial_activated_supps` and `cliques_term_sparsities` values.

---

## Design Questions

*(To be filled during review)*
