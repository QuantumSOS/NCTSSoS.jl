# Plan: Remove Algebra-Specific Logic from Monomial Multiplication

## Executive Summary

This plan outlines the refactoring to simplify monomial multiplication by removing ALL algebra-specific logic from `*`, `neat_dot`, and `_neat_dot3` functions. After this change, monomial multiplication will simply concatenate `word` fields, delegating all simplification to `simplify!`.

**Key Insight**: Currently, multiplication functions contain algebra-specific simplification logic (site-aware sorting, Pauli products, etc.). This creates tight coupling and makes the code harder to understand. By making multiplication ONLY concatenate words, we achieve clean separation: multiplication creates raw products, `simplify!` applies algebra rules.

## Current State Analysis

### Current Multiplication Implementations

#### 1. NonCommutativeAlgebra (site-aware)
**File**: `src/FastPolynomials/src/simplification/noncommutative.jl`
**Lines**: 163-173

```julia
function Base.:*(m1::Monomial{NonCommutativeAlgebra,T}, m2::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    w1, w2 = m1.word, m2.word

    # Handle empty cases
    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    # Concatenate and simplify using site-aware simplify!
    result = Monomial{NonCommutativeAlgebra,T}(vcat(w1, w2), zero(UInt64))
    simplify!(result)
end
```

**Current behavior**: Concatenates words AND calls `simplify!` which sorts by site.

#### 2. neat_dot and _neat_dot3 Helpers
**File**: `src/FastPolynomials/src/utils.jl`
**Lines**: 116-157

```julia
function neat_dot(a::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # adjoint(a) * b returns a Term, extract monomial
    result = adjoint(a) * b
    result.monomial
end

function _neat_dot3(a::Monomial{A,T}, m::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # (adjoint(a) * m) * b
    temp = adjoint(a) * m
    # temp is a Term, extract monomial for next multiplication
    result = temp.monomial * b
    # Return only the monomial
    result.monomial
end
```

**Current behavior**: Relies on `adjoint` + `*` which may call algebra-specific multiplication.

#### 3. Other Algebra Multiplications

The other algebra types (Pauli, Fermionic, Bosonic, Projector, Unipotent) do NOT have explicit `Base.:*` implementations in their simplification files. This means they fall back to... **NEED TO INVESTIGATE WHERE THE FALLBACK IS**.

Let me search for the default multiplication implementation.

### Usage Sites Analysis

From the grep results, `neat_dot` and `_neat_dot3` are used in:

1. **Moment matrix construction** - `src/moment_solver.jl`, `src/complex_moment_solver.jl`
   - Pattern: `simplify!(expval(_neat_dot3(row_idx, mono, col_idx)), sa)`
   - Creates entries like `⟨row_idx† · mono · col_idx⟩`

2. **Sparsity pattern computation** - `src/sparse.jl`
   - Pattern: `simplify!(neat_dot(b, b), sa)` and `_neat_dot3(bases[i], supp, bases[j])`
   - Used to compute connected supports

3. **GNS construction** - `src/gns.jl`
   - Pattern: `neat_dot(row_mono, col_mono)` for Hankel matrix keys
   - Pattern: `neat_dot(monomial(var), col_mono)` for multiplication matrices

4. **Tests** - Various test files use these helpers

**Critical observation**: ALL usage sites already call `simplify!` or `simplify` AFTER calling `neat_dot`/`_neat_dot3`. This means we can safely remove simplification from multiplication without breaking anything.

## Proposed Changes

### Phase 1: Remove Simplification from NonCommutativeAlgebra Multiplication

**File**: `src/FastPolynomials/src/simplification/noncommutative.jl`
**Lines to modify**: 163-173

**Before**:
```julia
function Base.:*(m1::Monomial{NonCommutativeAlgebra,T}, m2::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    w1, w2 = m1.word, m2.word

    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    result = Monomial{NonCommutativeAlgebra,T}(vcat(w1, w2), zero(UInt64))
    simplify!(result)  # <-- REMOVE THIS CALL
end
```

**After**:
```julia
function Base.:*(m1::Monomial{NonCommutativeAlgebra,T}, m2::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    w1, w2 = m1.word, m2.word

    # Empty cases: return identity coefficient with other monomial
    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    # Simple concatenation - no simplification
    result_word = vcat(w1, w2)
    result = Monomial{NonCommutativeAlgebra,T}(result_word, hash(result_word))
    return Term(1.0, result)
end
```

**Changes**:
- Remove `simplify!(result)` call
- Return `Term(1.0, result)` directly
- Concatenation only - no algebra-aware sorting

### Phase 2: Update All Other Algebra Multiplications

**CONFIRMED**: ALL algebra types have explicit multiplication implementations. They ALL follow the same pattern: `vcat(w1, w2)` + `simplify!`.

#### Files and Current Implementations:

**1. PauliAlgebra**
**File**: `src/FastPolynomials/src/simplification/pauli.jl`

```julia
function Base.:*(m1::Monomial{PauliAlgebra,T}, m2::Monomial{PauliAlgebra,T}) where {T}
    w1, w2 = m1.word, m2.word
    isempty(w1) && return simplify(m2)
    isempty(w2) && return simplify(m1)

    result = vcat(w1, w2)
    m_result = Monomial{PauliAlgebra}(result)
    simplify!(m_result)  # <-- CALLS SIMPLIFY
end
```
**Note**: Has TODO comment about performing simplification during multiplication.

**2. FermionicAlgebra**
**File**: `src/FastPolynomials/src/simplification/fermionic.jl`

```julia
function Base.:*(m1::Monomial{FermionicAlgebra,T}, m2::Monomial{FermionicAlgebra,T}) where {T}
    w1, w2 = m1.word, m2.word
    isempty(w1) && return simplify(m2)
    isempty(w2) && return simplify(m1)

    result = vcat(w1, w2)
    m_result = Monomial{FermionicAlgebra}(result)
    simplify!(m_result)  # <-- CALLS SIMPLIFY (Wick's theorem)
end
```

**3. BosonicAlgebra**
**File**: `src/FastPolynomials/src/simplification/bosonic.jl`

```julia
function Base.:*(m1::Monomial{BosonicAlgebra,T}, m2::Monomial{BosonicAlgebra,T}) where {T}
    w1, w2 = m1.word, m2.word
    isempty(w1) && return simplify(m2)
    isempty(w2) && return simplify(m1)

    result = vcat(w1, w2)
    m_result = Monomial{BosonicAlgebra}(result)
    simplify!(m_result)  # <-- CALLS SIMPLIFY (rook numbers)
end
```

**4. ProjectorAlgebra**
**File**: `src/FastPolynomials/src/simplification/projector.jl`

```julia
function Base.:*(m1::Monomial{ProjectorAlgebra,T}, m2::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}
    w1, w2 = m1.word, m2.word
    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    result = Monomial{ProjectorAlgebra,T}(vcat(w1, w2), zero(UInt64))
    simplify!(result)  # <-- CALLS SIMPLIFY
end
```

**5. UnipotentAlgebra** (has TWO versions!)
**File**: `src/FastPolynomials/src/simplification/unipotent.jl`

```julia
# Unsigned version
function Base.:*(m1::Monomial{UnipotentAlgebra,T}, m2::Monomial{UnipotentAlgebra,T}) where {T<:Unsigned}
    # Similar to Projector
    simplify!(result)  # <-- CALLS SIMPLIFY
end

# Signed version (legacy)
function Base.:*(m1::Monomial{UnipotentAlgebra,T}, m2::Monomial{UnipotentAlgebra,T}) where {T<:Signed}
    # Legacy version
    simplify!(result)  # <-- CALLS SIMPLIFY
end
```

**Pattern**: EVERY algebra type calls `simplify!` during multiplication. This is the coupling we want to remove.

### Phase 3: Simplify neat_dot and _neat_dot3

**File**: `src/FastPolynomials/src/utils.jl`
**Lines**: 116-157

**Current implementation**:
```julia
function neat_dot(a::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    result = adjoint(a) * b
    result.monomial
end

function _neat_dot3(a::Monomial{A,T}, m::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    temp = adjoint(a) * m
    result = temp.monomial * b
    result.monomial
end
```

**After refactoring**:

These functions will remain the SAME in terms of interface, but their behavior will change because:
1. `adjoint(a)` just reverses the word (possibly negating for signed types)
2. `*` will now just concatenate words
3. The result will be an UNSIMPLIFIED monomial

This is CORRECT because all call sites already apply `simplify!` afterwards.

**No code changes needed** - behavior changes due to updated `*` implementation.

### Phase 4: Update Tests

**Files to check**:
- `test/fastpoly_test/monomials.jl` - May expect simplified results from `*`
- `test/fastpoly_test/polynomial.jl` - Polynomial multiplication calls monomial `*`
- `test/fastpoly_test/arithmetic.jl` - Arithmetic tests
- `test/fastpoly_test/simplify.jl` - Simplification tests

**Expected test failures**:
Tests that check `m1 * m2` directly without calling `simplify!` may fail if they expect simplified results.

**Fix strategy**:
Wrap `*` results in `simplify!` calls in tests where simplified results are expected.

## Missing Information / Open Questions

### Q1: Where is the generic monomial multiplication fallback?

The grep results show only ONE explicit `Base.:*` for monomials (in NonCommutativeAlgebra). But multiplication must work for all algebra types. Where is the fallback?

**Action**: Search for generic multiplication in:
- `src/FastPolynomials/src/monomial.jl`
- `src/FastPolynomials/src/term.jl`

### Q2: Do Pauli/Fermionic/Bosonic algebras have multiplication implementations?

Need to check if they have `*` implementations that I missed, or if they rely on a generic fallback.

**Action**: Search more carefully in each simplification file for multiplication logic.

### Q3: What about Polynomial multiplication?

**File**: `src/FastPolynomials/src/polynomial.jl`
**Lines**: 670-692

Polynomial multiplication calls `t1.monomial * t2.monomial` and then handles the result:

```julia
simplified = t1.monomial * t2.monomial
_add_simplified_terms!(result_terms, coef, simplified)
```

This code is designed to handle both `Term` and `Vector{Term}` returns. After our refactoring:
- All monomial `*` will return `Term` with coefficient 1.0
- `simplify!` is called separately when needed
- Polynomial multiplication will continue to work

**Concern**: Polynomial `*` does NOT call `simplify!` on individual monomial products. Is this correct?

**Answer**: Polynomial multiplication creates raw terms from monomial products, then the Polynomial constructor processes terms (sorts, deduplicates, filters zeros). Simplification happens via explicit `simplify(p)` calls in user code, NOT during arithmetic.

**Impact**: NO CHANGE needed for Polynomial multiplication.

## Step-by-Step Implementation Plan

### Step 1: Update NonCommutativeAlgebra Multiplication ✓ CONFIRMED
- [ ] Edit `src/FastPolynomials/src/simplification/noncommutative.jl`
- [ ] Find function at line ~163-173
- [ ] Remove `simplify!(result)` call
- [ ] Return `Term(1.0, result)` with concatenated word only
- [ ] Update docstring to clarify no simplification happens

### Step 2: Update PauliAlgebra Multiplication ✓ CONFIRMED
- [ ] Edit `src/FastPolynomials/src/simplification/pauli.jl`
- [ ] Find `Base.:*(m1::Monomial{PauliAlgebra,T}, m2::Monomial{PauliAlgebra,T})`
- [ ] Remove `simplify!(m_result)` call
- [ ] Change empty cases to return `Term(1.0, m2)` instead of `simplify(m2)`
- [ ] Return `Term(1.0, result)` for normal case
- [ ] Remove TODO comment about simplification during multiplication (resolved!)

### Step 3: Update FermionicAlgebra Multiplication ✓ CONFIRMED
- [ ] Edit `src/FastPolynomials/src/simplification/fermionic.jl`
- [ ] Find `Base.:*(m1::Monomial{FermionicAlgebra,T}, m2::Monomial{FermionicAlgebra,T})`
- [ ] Remove `simplify!(m_result)` call
- [ ] Change empty cases to return `Term(1.0, m2)` instead of `simplify(m2)`
- [ ] Return `Term(1.0, result)` for normal case

### Step 4: Update BosonicAlgebra Multiplication ✓ CONFIRMED
- [ ] Edit `src/FastPolynomials/src/simplification/bosonic.jl`
- [ ] Find `Base.:*(m1::Monomial{BosonicAlgebra,T}, m2::Monomial{BosonicAlgebra,T})`
- [ ] Remove `simplify!(m_result)` call
- [ ] Change empty cases to return `Term(1.0, m2)` instead of `simplify(m2)`
- [ ] Return `Term(1.0, result)` for normal case

### Step 5: Update ProjectorAlgebra Multiplication ✓ CONFIRMED
- [ ] Edit `src/FastPolynomials/src/simplification/projector.jl`
- [ ] Find `Base.:*(m1::Monomial{ProjectorAlgebra,T}, m2::Monomial{ProjectorAlgebra,T})`
- [ ] Remove `simplify!(result)` call
- [ ] Return `Term(1.0, result)` with concatenated word only
- [ ] Already has correct empty case handling (returns Term, not simplify!)

### Step 6: Update UnipotentAlgebra Multiplication (2 versions!) ✓ CONFIRMED
- [ ] Edit `src/FastPolynomials/src/simplification/unipotent.jl`
- [ ] Find BOTH versions (Unsigned and Signed)
- [ ] Remove `simplify!(result)` calls from both
- [ ] Return `Term(1.0, result)` for both versions

### Step 7: Run Tests and Fix Failures
- [ ] Run `julia --project=. test/fastpoly_test/runtests.jl`
- [ ] Identify tests expecting simplified multiplication results
- [ ] Update tests to call `simplify!` explicitly when needed
- [ ] Document any unexpected failures

### Step 8: Verify neat_dot and _neat_dot3 Behavior
- [ ] Write test showing `neat_dot(a, b)` returns UNSIMPLIFIED result
- [ ] Verify all usage sites call `simplify!` afterwards
- [ ] Document that `neat_dot` is now "raw adjoint multiplication"

### Step 9: Update Documentation
- [ ] Update multiplication docstrings in all 6 algebra files to clarify "concatenation only"
- [ ] Update `simplify!` docstrings to clarify it's the ONLY simplification point
- [ ] Add architectural note in `utils.jl` explaining neat_dot behavior
- [ ] Update README or architecture docs if they mention multiplication

## Testing Strategy

### Unit Tests
For each algebra type, verify:
```julia
m1 = Monomial{AlgebraType}([...])
m2 = Monomial{AlgebraType}([...])
term = m1 * m2
@test term.coefficient == 1.0
@test term.monomial.word == vcat(m1.word, m2.word)  # Raw concatenation

# Simplification happens separately
simplified_term = simplify!(term.monomial)
# ... check simplified result
```

### Integration Tests
Check that moment solver and sparse code still work:
```julia
# From moment_solver.jl pattern
row_idx = ...
mono = ...
col_idx = ...
result = simplify!(expval(_neat_dot3(row_idx, mono, col_idx)), sa)
# Should produce same results as before
```

## Risks and Mitigation

### Risk 1: Performance Regression
**Concern**: Calling `simplify!` separately might be slower than simplifying during multiplication.

**Mitigation**:
- Benchmark before/after
- Most code already calls `simplify!` explicitly, so no change
- Potential speedup from simpler multiplication

### Risk 2: Test Failures
**Concern**: Many tests may expect simplified multiplication results.

**Mitigation**:
- Run tests incrementally
- Fix tests to call `simplify!` when expecting simplified results
- Document which tests changed and why

### Risk 3: Breaking Higher-Level Code
**Concern**: Some code might assume `*` returns simplified results.

**Mitigation**:
- Grep for all uses of monomial `*`
- Verify they either (a) don't care about simplification or (b) call `simplify!`
- The pattern `simplify!(..._neat_dot3(...))` is already universal

## Success Criteria

- [ ] All monomial `*` implementations are pure concatenation
- [ ] All tests pass
- [ ] No performance regressions
- [ ] Code is simpler and easier to understand
- [ ] Clear separation: `*` creates products, `simplify!` applies algebra rules

## References

### Files to Modify
1. `src/FastPolynomials/src/simplification/noncommutative.jl` - Remove simplify from `*`
2. `src/FastPolynomials/src/utils.jl` - Document new behavior of `neat_dot`/_`neat_dot3`
3. `test/fastpoly_test/*.jl` - Update tests expecting simplified multiplication

### Files to Review (No Changes Expected)
1. `src/FastPolynomials/src/polynomial.jl` - Polynomial `*` already handles raw terms
2. `src/moment_solver.jl` - Already calls `simplify!` on `_neat_dot3` results
3. `src/sparse.jl` - Already calls `simplify!` on `neat_dot` results

### Search Patterns Used
```bash
rg "function Base\.:\*" src/FastPolynomials/src/simplification/
rg "neat_dot" --type julia
rg "_neat_dot3" --type julia
rg "Base\.\*\(" --type julia -A 5
```

## Summary: Complete File List for Refactoring

### Files to MODIFY (6 algebra simplification files)
1. **`src/FastPolynomials/src/simplification/noncommutative.jl`**
   - Line ~163: `Base.:*(::Monomial{NonCommutativeAlgebra,T<:Unsigned}, ...)`

2. **`src/FastPolynomials/src/simplification/pauli.jl`**
   - Find: `Base.:*(::Monomial{PauliAlgebra,T}, ...)`
   - Extra: Remove TODO comment about simplification during multiplication

3. **`src/FastPolynomials/src/simplification/fermionic.jl`**
   - Find: `Base.:*(::Monomial{FermionicAlgebra,T}, ...)`

4. **`src/FastPolynomials/src/simplification/bosonic.jl`**
   - Find: `Base.:*(::Monomial{BosonicAlgebra,T}, ...)`

5. **`src/FastPolynomials/src/simplification/projector.jl`**
   - Find: `Base.:*(::Monomial{ProjectorAlgebra,T<:Unsigned}, ...)`

6. **`src/FastPolynomials/src/simplification/unipotent.jl`**
   - Find TWO functions:
     - `Base.:*(::Monomial{UnipotentAlgebra,T<:Unsigned}, ...)`
     - `Base.:*(::Monomial{UnipotentAlgebra,T<:Signed}, ...)`

### Files to UPDATE DOCS ONLY
7. **`src/FastPolynomials/src/utils.jl`**
   - Update docstrings for `neat_dot` and `_neat_dot3`
   - Add architectural note explaining concatenation-only behavior

### Test Files (may need updates)
- `test/fastpoly_test/monomials.jl`
- `test/fastpoly_test/polynomial.jl`
- `test/fastpoly_test/arithmetic.jl`
- `test/fastpoly_test/simplify.jl`

## Change Pattern for All 6 Algebra Files

**Before** (example from any algebra):
```julia
function Base.:*(m1::Monomial{SomeAlgebra,T}, m2::Monomial{SomeAlgebra,T}) where {T}
    w1, w2 = m1.word, m2.word

    # Empty cases call simplify
    isempty(w1) && return simplify(m2)
    isempty(w2) && return simplify(m1)

    # Concatenate and simplify
    result = vcat(w1, w2)
    m_result = Monomial{SomeAlgebra}(result)
    simplify!(m_result)  # <-- REMOVE THIS
end
```

**After** (uniform across all algebras):
```julia
function Base.:*(m1::Monomial{SomeAlgebra,T}, m2::Monomial{SomeAlgebra,T}) where {T}
    w1, w2 = m1.word, m2.word

    # Empty cases return identity coefficient
    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    # Pure concatenation - no simplification
    result_word = vcat(w1, w2)
    result = Monomial{SomeAlgebra}(result_word)
    return Term(1.0, result)  # <-- SIMPLE RETURN
end
```

**Key changes**:
1. Remove `simplify!(...)` call
2. Change empty case returns from `simplify(m)` to `Term(1.0, m)`
3. Return `Term(1.0, result)` instead of simplified term
4. Update docstring to clarify "concatenation only, no simplification"
