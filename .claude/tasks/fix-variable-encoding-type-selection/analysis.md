# Detailed Analysis: Variable Encoding Type Selection Bug

## Executive Summary

The `_select_contiguous_type` function incorrectly selects UInt8 for 64 variables because it checks `64 <= typemax(UInt8)` (255), but the encoding scheme `(operator_id << 2) | site` can only fit operator_id values up to 63 in UInt8's 6 available bits.

**Fix:** Use the existing `_select_unsigned_type(n_operators, n_sites)` function which correctly validates both encoding constraints.

---

## Visual Breakdown: Why UInt8 Fails at n=64

### Bit Layout for UInt8 Encoding

```
UInt8 (8 bits total)
┌──────────┬──┐
│ operator │site│
│  (6 bits)│(2) │
└──────────┴──┘
  bits 7-2  1-0

Max values:
- operator_id: 0-63 (2^6 - 1)
- site: 0-3 (2^2 - 1, but 1-indexed so 1-3)
```

### Encoding Process for operator_id=64, site=1

```julia
_encode_index(UInt8, 64, 1)
= UInt8(64 << 2 | 1)
= UInt8(256 | 1)
= UInt8(257)  # ❌ OVERFLOW! 257 > typemax(UInt8) = 255
```

**Visual:**
```
operator_id = 64 (binary: 01000000)
Shift left by 2:  0100000000  # Now 10 bits!
OR with site=1:   0100000001  # 257 in decimal
                  ^^^^^^^^^^
                  10 bits needed, but UInt8 only has 8!
```

### Why operator_id=63 Works

```julia
_encode_index(UInt8, 63, 1)
= UInt8(63 << 2 | 1)
= UInt8(252 | 1)
= UInt8(253)  # ✓ OK! 253 <= 255
```

**Visual:**
```
operator_id = 63 (binary: 00111111)
Shift left by 2:  11111100  # Exactly 8 bits
OR with site=1:   11111101  # 253 in decimal ✓
```

---

## Encoding Limits Table

| Type   | sizeof(T) | site_bits | operator_bits | typemax(T) | max_operators(T) | max_sites(T) |
|--------|-----------|-----------|---------------|------------|------------------|--------------|
| UInt8  | 1 byte    | 2         | 6             | 255        | **63**           | 3            |
| UInt16 | 2 bytes   | 4         | 12            | 65535      | **4095**         | 15           |
| UInt32 | 4 bytes   | 8         | 24            | 4294967295 | **16777215**     | 255          |
| UInt64 | 8 bytes   | 16        | 48            | 2^64-1     | **281474976710655** | 65535     |

**Key Insight:** The "usable" limit is `max_operators(T)`, NOT `typemax(T)`.

---

## Code Comparison: Buggy vs Correct

### Buggy: `_select_contiguous_type` (current)

```julia
function _select_contiguous_type(n::Int)
    n <= typemax(UInt8) && return UInt8    # 255 ❌ WRONG
    n <= typemax(UInt16) && return UInt16  # 65535 ❌ WRONG
    n <= typemax(UInt32) && return UInt32  # 4294967295 ❌ WRONG
    return UInt64
end
```

**Problems:**
1. Checks wrong limit (typemax vs max_operators)
2. Ignores site constraint entirely
3. No error handling for edge cases

### Correct: `_select_unsigned_type` (already exists!)

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

**Benefits:**
1. Checks correct limits (max_operators and max_sites)
2. Validates BOTH encoding constraints
3. Clear error message if impossible

---

## Mathematical Derivation: Why sizeof(T)*2 for site_bits?

The design allocates **1/4 of total bits** to site encoding:

```
site_bits(T) = sizeof(T) * 2
             = (byte_count) * 2
             = (bit_count / 8) * 2
             = bit_count / 4
```

**Examples:**
- UInt8: 8 bits / 4 = 2 bits for site
- UInt16: 16 bits / 4 = 4 bits for site
- UInt32: 32 bits / 4 = 8 bits for site
- UInt64: 64 bits / 4 = 16 bits for site

This leaves **3/4 of total bits** for operator_id:

```
operator_bits(T) = sizeof(T) * 8 - site_bits(T)
                 = bit_count - (bit_count / 4)
                 = 3 * bit_count / 4
```

**Design Rationale:** Most applications have many operators per site, but few total sites, so the 3:1 ratio is a reasonable trade-off.

---

## Call Graph: How the Bug is Triggered

```
User calls:
  create_projector_variables([("P", 1:64)])
    ↓
  _create_noncommutative_variables(ProjectorAlgebra, [("P", 1:64)])
    ↓
  n = sum(length(x[2]) for x in prefix_subscripts) = 64
    ↓
  IndexT = _select_contiguous_type(64)
    ↓
  Returns UInt8 (because 64 <= 255) ❌ BUG!
    ↓
  for global_idx in 1:64, physical_site = 1
    encoded_idx = _encode_index(UInt8, global_idx, physical_site)
    ↓ (when global_idx = 64)
    = UInt8(64 << 2 | 1)
    = UInt8(257)  ❌ InexactError!
```

**Fixed call graph:**

```
User calls:
  create_projector_variables([("P", 1:64)])
    ↓
  _create_noncommutative_variables(ProjectorAlgebra, [("P", 1:64)])
    ↓
  n_operators = 64, n_sites = 1
    ↓
  IndexT = _select_unsigned_type(64, 1)
    ↓
  Check UInt8: 64 <= max_operators(UInt8)? 64 <= 63? NO
  Check UInt16: 64 <= max_operators(UInt16)? 64 <= 4095? YES
  Returns UInt16 ✓ CORRECT!
    ↓
  for global_idx in 1:64, physical_site = 1
    encoded_idx = _encode_index(UInt16, global_idx, physical_site)
    ↓ (when global_idx = 64)
    = UInt16(64 << 4 | 1)
    = UInt16(1024 | 1)
    = UInt16(1025)  ✓ OK! 1025 <= 65535
```

---

## Test Cases: Boundary Conditions

### Test 1: Operator Limit for UInt8 (63 vs 64)

```julia
# Exactly at limit: should use UInt8
reg63, (x,) = create_noncommutative_variables([("x", 1:63)])
@test index_type(reg63) == UInt8

# One over limit: should use UInt16
reg64, (x,) = create_noncommutative_variables([("x", 1:64)])
@test index_type(reg64) == UInt16
```

### Test 2: Site Limit for UInt8 (3 vs 4 groups)

```julia
# Exactly at limit: should use UInt8
reg3, (P, Q, R) = create_projector_variables([
    ("P", 1:10), ("Q", 1:10), ("R", 1:10)
])
@test index_type(reg3) == UInt8

# One over limit: should use UInt16
reg4, (P, Q, R, S) = create_projector_variables([
    ("P", 1:10), ("Q", 1:10), ("R", 1:10), ("S", 1:10)
])
@test index_type(reg4) == UInt16
```

### Test 3: Combined Constraints

```julia
# 10 groups × 50 vars = 500 total
# UInt8: max_sites=3 ❌ (need 10)
# UInt16: max_sites=15 ✓, max_operators=4095 ✓
reg, _ = create_projector_variables([
    (string('A' + i - 1), 1:50) for i in 1:10
])
@test index_type(reg) == UInt16
@test length(reg) == 500
```

### Test 4: Operator Limit for UInt16 (4095 vs 4096)

```julia
# Exactly at limit: should use UInt16
reg4095, (x,) = create_noncommutative_variables([("x", 1:4095)])
@test index_type(reg4095) == UInt16

# One over limit: should use UInt32
reg4096, (x,) = create_noncommutative_variables([("x", 1:4096)])
@test index_type(reg4096) == UInt32
```

---

## Implementation Diff Preview

### File: `src/FastPolynomials/src/variable_registry.jl`

**DELETE lines 651-661:**
```julia
-"""
-    _select_contiguous_type(n::Int) -> Type{<:Unsigned}
-
-Select smallest unsigned type for contiguous 1-indexed variables.
-"""
-function _select_contiguous_type(n::Int)
-    n <= typemax(UInt8) && return UInt8
-    n <= typemax(UInt16) && return UInt16
-    n <= typemax(UInt32) && return UInt32
-    return UInt64
-end
```

**MODIFY lines 691-692:**
```diff
  function _create_noncommutative_variables(
      ::Type{A},
      prefix_subscripts::Vector{Tuple{String, VT}}
  ) where {A<:AlgebraType, T<:Integer, VT<:AbstractVector{T}}
-     n = sum(x -> length(x[2]), prefix_subscripts)
-     IndexT = _select_contiguous_type(n)
+     n_operators = sum(x -> length(x[2]), prefix_subscripts)
+     n_sites = length(prefix_subscripts)
+     IndexT = _select_unsigned_type(n_operators, n_sites)

      all_symbols = Vector{Symbol}(undef, n_operators)
      all_indices = Vector{IndexT}(undef, n_operators)
```

### File: `test/fastpoly_test/variables.jl`

**MODIFY line 6 (remove import):**
```diff
  using NCTSSoS.FastPolynomials:
      _subscript_string,
      _select_pauli_type,
-     _select_signed_index_type,
-     _select_contiguous_type
+     _select_signed_index_type
```

**DELETE lines 112-114:**
```diff
-     # _select_contiguous_type
-     @test _select_contiguous_type(255) == UInt8
-     @test _select_contiguous_type(256) == UInt16
```

---

## Verification Checklist

After implementing the fix:

**Functional Tests:**
- [ ] `create_noncommutative_variables([("x", 1:64)])` succeeds (was failing)
- [ ] Registry uses UInt16 for 64 variables (not UInt8)
- [ ] Boundary test: 63 vars → UInt8, 64 vars → UInt16
- [ ] Boundary test: 3 groups → UInt8, 4 groups → UInt16
- [ ] All existing tests pass: `make test`

**Code Quality:**
- [ ] No references to `_select_contiguous_type` remain
- [ ] No compiler warnings
- [ ] Git commit message follows convention: `fix(FastPolynomials): correct type selection for encoded indices`

**Documentation:**
- [ ] No broken docstring references (run `rg "_select_contiguous_type" --type md`)
- [ ] Plan document archived in `.claude/tasks/fix-variable-encoding-type-selection/`

---

## References

**Source Code:**
- [Source: `/src/FastPolynomials/src/variable_registry.jl:663`] — `_encode_index` definition
- [Source: `/src/FastPolynomials/src/algebra_types.jl:186`] — `site_bits` definition
- [Source: `/src/FastPolynomials/src/algebra_types.jl:201`] — `max_sites` definition
- [Source: `/src/FastPolynomials/src/algebra_types.jl:215`] — `max_operators` definition
- [Source: `/src/FastPolynomials/src/variable_registry.jl:563`] — Correct `_select_unsigned_type` implementation

**Test Evidence:**
- [Source: `/test/fastpoly_test/algebra_types.jl:181-212`] — Tests validating `max_operators` and `max_sites` formulas
- [Source: `/test_bug_verification.jl`] — User-created test reproducing the bug

**Design Documentation:**
- [Source: `/src/FastPolynomials/src/algebra_types.jl:162-174`] — Comments explaining site-based index encoding scheme
