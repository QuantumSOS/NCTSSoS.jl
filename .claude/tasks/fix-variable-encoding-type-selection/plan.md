# Implementation Plan: Fix Variable Encoding Type Selection Bug

## 1. Analysis: Understanding the Encoding Scheme

### 1.1 The Encoding Mathematics

The codebase uses a **bit-packed encoding scheme** for variable indices in algebras that support site-based commutation (ProjectorAlgebra, UnipotentAlgebra, NonCommutativeAlgebra).

**Encoding Formula** (from `algebra_types.jl` line 243):
```julia
encoded_index = (operator_id << site_bits) | site
```

**Bit Allocation** (from `algebra_types.jl` line 186):
```julia
site_bits(T) = sizeof(T) * 2  # n/4 where n is bit width
```

This means:
- **UInt8** (8 bits total):
  - Site bits: `sizeof(UInt8) * 2 = 1 * 2 = 2 bits`
  - Operator bits: `8 - 2 = 6 bits`
  - Max sites: `2^2 - 1 = 3` (values 1, 2, 3)
  - Max operators: `2^6 - 1 = 63` (values 1-63)

- **UInt16** (16 bits total):
  - Site bits: `sizeof(UInt16) * 2 = 2 * 2 = 4 bits`
  - Operator bits: `16 - 4 = 12 bits`
  - Max sites: `2^4 - 1 = 15` (values 1-15)
  - Max operators: `2^12 - 1 = 4095` (values 1-4095)

- **UInt32** (32 bits total):
  - Site bits: `sizeof(UInt32) * 2 = 4 * 2 = 8 bits`
  - Operator bits: `32 - 8 = 24 bits`
  - Max sites: `2^8 - 1 = 255` (values 1-255)
  - Max operators: `2^24 - 1 = 16777215` (values 1-16777215)

- **UInt64** (64 bits total):
  - Site bits: `sizeof(UInt64) * 2 = 8 * 2 = 16 bits`
  - Operator bits: `64 - 16 = 48 bits`
  - Max sites: `2^16 - 1 = 65535` (values 1-65535)
  - Max operators: `2^48 - 1 = 281474976710655` (values 1-281474976710655)

### 1.2 The Bug: Type Selection Ignores Encoding Constraints

**Current buggy implementation** (`variable_registry.jl` lines 656-661):
```julia
function _select_contiguous_type(n::Int)
    n <= typemax(UInt8) && return UInt8    # 255 - WRONG!
    n <= typemax(UInt16) && return UInt16  # 65535 - WRONG!
    n <= typemax(UInt32) && return UInt32  # 4294967295 - WRONG!
    return UInt64
end
```

This function incorrectly uses `typemax(T)` which represents the maximum VALUE that can be stored in the type, NOT the maximum number of operators that can be encoded.

**Why it fails:**
- When creating 64 variables in a single group:
  - `n = 64`
  - `_select_contiguous_type(64)` returns `UInt8` (since 64 ≤ 255)
  - But `_encode_index(UInt8, 64, 1)` tries to compute: `64 << 2 | 1 = 256 | 1 = 257`
  - `UInt8(257)` throws `InexactError` because 257 > 255 (max value for UInt8)

**Root cause:** The function should check against `max_operators(T)` (63 for UInt8), not `typemax(T)` (255 for UInt8).

### 1.3 Additional Constraint: Multiple Prefix Groups

The `_create_noncommutative_variables` function (line 687-727) supports **multiple prefix groups**, where each group is assigned a different `physical_site` value.

**Example:**
```julia
create_projector_variables([("P", 1:100), ("Q", 1:100)])
# Creates P₁...P₁₀₀ (site 1) and Q₁...Q₁₀₀ (site 2)
# n_total = 200 variables
# n_groups = 2 physical sites
```

**This means the type selection must satisfy TWO constraints:**
1. `n_total ≤ max_operators(T)` — enough bits for operator indices
2. `n_groups ≤ max_sites(T)` — enough bits for site indices

**Current implementation only checks constraint #1 indirectly (and incorrectly).**

### 1.4 Relationship to `_select_unsigned_type`

The codebase ALREADY has a correct implementation in `algebra_types.jl` (lines 263-269) and duplicated in `variable_registry.jl` (lines 563-570):

```julia
function _select_unsigned_type(n_operators::Int, n_sites::Int)
    for T in (UInt8, UInt16, UInt32, UInt64)
        if n_sites <= max_sites(T) && n_operators <= max_operators(T)
            return T
        end
    end
    error("Cannot fit $n_operators operators × $n_sites sites in any UInt type")
end
```

This function correctly uses `max_operators(T)` and `max_sites(T)`.

**The bug is that `_select_contiguous_type` does NOT use this approach.**

## 2. Proposed Solution

### 2.1 Fix Strategy

Replace the buggy `_select_contiguous_type` function with logic that:
1. Accepts BOTH `n_operators` (total variables) AND `n_groups` (number of prefix groups)
2. Checks both `max_operators(T)` and `max_sites(T)` constraints
3. Returns the smallest `T` that satisfies both constraints
4. Throws a clear error if no type can accommodate the request

### 2.2 Pseudocode for Corrected Function

**Option A: Create a new dual-constraint function**
```julia
function _select_encoded_type(n_operators::Int, n_sites::Int)
    for T in (UInt8, UInt16, UInt32, UInt64)
        if n_operators <= max_operators(T) && n_sites <= max_sites(T)
            return T
        end
    end
    error("Cannot fit $n_operators operators across $n_sites sites in any UInt type")
end
```

**Option B: Unify with existing `_select_unsigned_type`**

Since `_select_unsigned_type` already implements the correct logic and appears twice in the codebase:
- Use `_select_unsigned_type` directly in `_create_noncommutative_variables`
- Deprecate `_select_contiguous_type` entirely

**Recommended: Option B** (reduces duplication, uses proven implementation)

### 2.3 Call Site Modification

**Current code** (`_create_noncommutative_variables`, line 692):
```julia
n = sum(x -> length(x[2]), prefix_subscripts)
IndexT = _select_contiguous_type(n)
```

**Fixed code:**
```julia
n_operators = sum(x -> length(x[2]), prefix_subscripts)
n_sites = length(prefix_subscripts)
IndexT = _select_unsigned_type(n_operators, n_sites)
```

**Key changes:**
1. Rename `n` to `n_operators` for clarity
2. Extract `n_sites = length(prefix_subscripts)` to count physical sites
3. Call `_select_unsigned_type(n_operators, n_sites)` instead of `_select_contiguous_type(n)`

### 2.4 Handling Edge Cases

**Edge Case 1: Single prefix group with many variables**
```julia
create_projector_variables([("P", 1:100)])
# n_operators = 100, n_sites = 1
# UInt8: max_operators = 63 → FAIL
# UInt16: max_operators = 4095, max_sites = 15 → PASS
# Returns UInt16 ✓
```

**Edge Case 2: Many prefix groups with few variables each**
```julia
create_projector_variables([("P", 1:2), ("Q", 1:2), ..., ("Z", 1:2)])  # 10 groups
# n_operators = 20, n_sites = 10
# UInt8: max_sites = 3 → FAIL
# UInt16: max_operators = 4095, max_sites = 15 → PASS
# Returns UInt16 ✓
```

**Edge Case 3: Massive number of operators**
```julia
create_projector_variables([("P", 1:100000)])
# n_operators = 100000, n_sites = 1
# UInt8: max_operators = 63 → FAIL
# UInt16: max_operators = 4095 → FAIL
# UInt32: max_operators = 16777215 → PASS
# Returns UInt32 ✓
```

**Edge Case 4: Too many sites for any type**
```julia
create_projector_variables([(string(i), 1:1) for i in 1:100000])
# n_operators = 100000, n_sites = 100000
# UInt64: max_sites = 65535 → FAIL
# Throws error: "Cannot fit 100000 operators across 100000 sites in any UInt type"
```

## 3. Implementation Steps

### Step 1: Remove `_select_contiguous_type` function
**File:** `src/FastPolynomials/src/variable_registry.jl`
**Lines:** 651-661 (documentation + implementation)
**Action:** Delete this function entirely

### Step 2: Update `_create_noncommutative_variables` to use `_select_unsigned_type`
**File:** `src/FastPolynomials/src/variable_registry.jl`
**Lines:** 687-727
**Changes:**
- Line 691: Rename `n` to `n_operators` (variable name)
- Line 691: Add new line: `n_sites = length(prefix_subscripts)`
- Line 692: Replace `_select_contiguous_type(n)` with `_select_unsigned_type(n_operators, n_sites)`
- Update all subsequent references to `n` → `n_operators` (lines 694, 695)

### Step 3: Remove test for deprecated function
**File:** `test/fastpoly_test/variables.jl`
**Lines:** 112-114
**Action:** Remove these test lines:
```julia
# _select_contiguous_type
@test _select_contiguous_type(255) == UInt8
@test _select_contiguous_type(256) == UInt16
```

**Rationale:** The function is being removed, so its tests should also be removed.

### Step 4: Remove import of deprecated function
**File:** `test/fastpoly_test/variables.jl`
**Line:** 6
**Action:** Remove `_select_contiguous_type` from the import statement:
```julia
# Before:
using NCTSSoS.FastPolynomials:
    _subscript_string,
    _select_pauli_type,
    _select_signed_index_type,
    _select_contiguous_type

# After:
using NCTSSoS.FastPolynomials:
    _subscript_string,
    _select_pauli_type,
    _select_signed_index_type
```

## 4. Testing Strategy

### 4.1 Verification of Fix

**Test Case 1: Reproduce original bug (should now pass)**
```julia
# Previously failed at n=64 for UInt8 (max_operators = 63)
reg, (x,) = create_noncommutative_variables([("x", 1:64)])
@test index_type(reg) == UInt16  # Should select UInt16
@test length(reg) == 64
```

**Test Case 2: Boundary at max_operators for UInt8**
```julia
# n=63 should use UInt8 (exactly at limit)
reg63, (x,) = create_noncommutative_variables([("x", 1:63)])
@test index_type(reg63) == UInt8
@test length(reg63) == 63

# n=64 should use UInt16 (exceeds UInt8 limit)
reg64, (x,) = create_noncommutative_variables([("x", 1:64)])
@test index_type(reg64) == UInt16
@test length(reg64) == 64
```

**Test Case 3: Multiple prefix groups hitting site limit**
```julia
# 3 groups should use UInt8 (max_sites = 3)
reg3, (P, Q, R) = create_projector_variables([("P", 1:2), ("Q", 1:2), ("R", 1:2)])
@test index_type(reg3) == UInt8

# 4 groups should use UInt16 (exceeds UInt8 max_sites = 3)
reg4, (P, Q, R, S) = create_projector_variables([("P", 1:2), ("Q", 1:2), ("R", 1:2), ("S", 1:2)])
@test index_type(reg4) == UInt16
```

**Test Case 4: Boundary at max_operators for UInt16**
```julia
# n=4095 should use UInt16 (exactly at limit)
reg4095, (x,) = create_noncommutative_variables([("x", 1:4095)])
@test index_type(reg4095) == UInt16

# n=4096 should use UInt32 (exceeds UInt16 limit)
reg4096, (x,) = create_noncommutative_variables([("x", 1:4096)])
@test index_type(reg4096) == UInt32
```

**Test Case 5: Combined constraints**
```julia
# 10 groups × 50 variables = 500 total
# UInt8: max_sites = 3 → FAIL (need 10 sites)
# UInt16: max_sites = 15 → PASS
reg, groups = create_projector_variables([(string('A' + i - 1), 1:50) for i in 1:10])
@test index_type(reg) == UInt16
@test length(reg) == 500
```

### 4.2 Regression Testing

Run the full test suite to ensure no existing functionality is broken:
```bash
LOCAL_TESTING=true make test
```

Specifically check:
- `test/fastpoly_test/variables.jl` — Variable creation tests
- `test/fastpoly_test/algebra_types.jl` — Encoding/decoding tests
- All algebra-specific variable creation (Pauli, Fermionic, Bosonic, Projector, Unipotent, NonCommutative)

### 4.3 Documentation Verification

Ensure that removing `_select_contiguous_type` doesn't break:
1. Any docstrings that reference it (unlikely, as it's internal)
2. Any examples in documentation (unlikely, as it's internal)

Search for references:
```bash
rg "_select_contiguous_type" --type md
```

## 5. Edge Cases and Error Handling

### 5.1 Error Message Quality

The existing `_select_unsigned_type` function (line 569) provides a clear error:
```julia
error("Cannot fit $n_operators operators × $n_sites sites in any UInt type")
```

This is superior to a silent overflow or cryptic `InexactError`.

### 5.2 Performance Considerations

**Current approach:**
- Iterates through (UInt8, UInt16, UInt32, UInt64) checking constraints
- Early exit on first match
- O(1) in practice (max 4 iterations)

**No performance regression expected** — the corrected logic is equally efficient.

### 5.3 Type Stability

The function returns a **Type** (UInt8, UInt16, UInt32, or UInt64), not a value. This is used at compile-time for type parameterization, so type stability is maintained.

## 6. Summary of Changes

| File | Lines | Action |
|------|-------|--------|
| `src/FastPolynomials/src/variable_registry.jl` | 651-661 | **DELETE** `_select_contiguous_type` function |
| `src/FastPolynomials/src/variable_registry.jl` | 691-692 | **MODIFY** `_create_noncommutative_variables` to use `_select_unsigned_type` |
| `test/fastpoly_test/variables.jl` | 6 | **MODIFY** Remove `_select_contiguous_type` from imports |
| `test/fastpoly_test/variables.jl` | 112-114 | **DELETE** Tests for deprecated function |
| **(NEW)** `test/fastpoly_test/variables.jl` | N/A | **ADD** New tests for encoding boundary conditions |

## 7. References

**Internal Functions:**
- `max_operators(::Type{T})` — `src/FastPolynomials/src/algebra_types.jl:215`
- `max_sites(::Type{T})` — `src/FastPolynomials/src/algebra_types.jl:201`
- `_select_unsigned_type(n_operators, n_sites)` — `src/FastPolynomials/src/variable_registry.jl:563`
- `_encode_index(T, global_idx, physical_site)` — `src/FastPolynomials/src/variable_registry.jl:663`
- `encode_index(T, operator_id, site)` — `src/FastPolynomials/src/algebra_types.jl:236` (public version with assertions)

**Test Evidence:**
- `test/fastpoly_test/algebra_types.jl` — Tests for `max_operators` and `max_sites` (lines 181-212)
- `test_bug_verification.jl` — User-created test demonstrating the bug

## 8. Post-Implementation Checklist

After implementing the fix, verify:

- [ ] All tests pass: `LOCAL_TESTING=true make test`
- [ ] Bug verification test passes (64 variables in single group)
- [ ] Boundary tests for UInt8 (63 vs 64 variables) pass
- [ ] Boundary tests for UInt8 (3 vs 4 prefix groups) pass
- [ ] No regressions in existing variable creation tests
- [ ] Code compiles without warnings
- [ ] Git commit message follows conventional commit format: `fix(FastPolynomials): use max_operators limit in type selection`
