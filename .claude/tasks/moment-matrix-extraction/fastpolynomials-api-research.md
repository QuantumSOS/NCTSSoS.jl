# FastPolynomials API Research: Monomial Operations

**Date**: 2025-11-16
**Purpose**: Understand FastPolynomials API for implementing monomial adjoint products in moment matrix extraction

---

## Executive Summary

FastPolynomials provides all necessary operations for computing adjoint products of monomials:
- **`star(mono)`**: Computes adjoint/conjugate by reversing variable order and exponents
- **`neat_dot(mono1, mono2)`**: Computes `star(mono1) * mono2` efficiently
- **`_neat_dot3(m1, m2, m3)`**: Computes `star(m1) * m2 * m3` efficiently
- **Direct multiplication**: `mono1 * mono2` concatenates monomials

The **recommended approach** is to use the existing `neat_dot` function directly, as it already implements the exact operation needed: `mono1^† * mono2`.

---

## 1. Core Operations

### 1.1 Monomial Adjoint: `star()`

**Location**: `/Users/yushengzhao/projects/NCTSSoS.jl/src/FastPolynomials/src/monomials.jl` (lines 126-140)

**Function Signature**:
```julia
star(m::Monomial) -> Monomial
```

**Implementation**:
```julia
"""
    star(m::Monomial)

Computes the adjoint (star) of a monomial by reversing variable order and exponents.

# Arguments
- `m::Monomial`: The monomial to compute the adjoint of

# Returns
- `Monomial`: Adjoint monomial with reversed variables and exponents
"""
star(m::Monomial) = star!(copy(m))

function star!(m::Monomial)
    (length(m.vars) <= 1 && return m; reverse!(m.vars); reverse!(m.z); return m)
end
```

**Usage Example** (from test file):
```julia
@ncpolyvar x y z
mono = monomial([x, y, z], [1, 1, 1])   # xyz
mono_star = star(mono)
# Result: mono_star.vars == [z, y, x]
#         mono_star.z == [1, 1, 1]
```

**Important Notes**:
- Reverses **both** the variable list and exponent list
- Assumes all variables are Hermitian (self-adjoint)
- In-place version `star!()` mutates the monomial
- For single variable or empty monomial, returns unchanged

---

### 1.2 Monomial Multiplication: `*`

**Location**: `/Users/yushengzhao/projects/NCTSSoS.jl/src/FastPolynomials/src/arithmetic.jl` (lines 80-84)

**Function Signature**:
```julia
Base.:(*)(x::Monomial, y::Monomial) -> Monomial
```

**Implementation**:
```julia
function Base.:(*)(x::Monomial, y::Monomial)
    isempty(x.z) && return y
    isempty(y.z) && return x
    return Monomial(_concat_var_expos(x.vars, x.z, y.vars, y.z)...)
end
```

**Usage Example**:
```julia
@ncpolyvar x y z
mono1 = monomial([x, y], [1, 2])  # xy²
mono2 = monomial([x, z], [3, 4])  # x³z⁴
result = mono1 * mono2
# Result: monomial([x, y, x, z], [1, 2, 3, 4])
#         which simplifies to xy²x³z⁴
```

**Key Feature**: `_concat_var_expos()` handles automatic consolidation when last variable of first monomial equals first variable of second monomial.

---

### 1.3 Adjoint Product: `neat_dot()`

**Location**: `/Users/yushengzhao/projects/NCTSSoS.jl/src/FastPolynomials/src/arithmetic.jl` (lines 87-99)

**Function Signature**:
```julia
neat_dot(x::Monomial, y::Monomial) -> Monomial
```

**Implementation**:
```julia
"""
    neat_dot(x::Monomial, y::Monomial)

Computes the "neat dot" product of two monomials as star(x) * y.

# Arguments
- `x::Monomial`: First monomial
- `y::Monomial`: Second monomial

# Returns
- `Monomial`: Product of star(x) and y
"""
neat_dot(x::Monomial, y::Monomial) =
    Monomial(_concat_var_expos(reverse(x.vars), reverse(x.z), y.vars, y.z)...)
```

**Usage Example** (from test file):
```julia
@ncpolyvar x y z
mono1 = monomial([x, y], [1, 0])  # x
mono2 = monomial([x, y], [1, 1])  # xy

result = neat_dot(mono1, mono2)
# Result: monomial([x, y], [2, 1])
#         This is star(x) * xy = x * xy = x²y

result2 = neat_dot(mono2, mono2)
# Result: monomial([y, x, y], [1, 2, 1])
#         This is star(xy) * xy = yx * xy = yxxy = yx²y
```

**Performance Note**: More efficient than `star(x) * y` because it avoids intermediate allocation.

---

### 1.4 Three-Way Adjoint Product: `_neat_dot3()`

**Location**: `/Users/yushengzhao/projects/NCTSSoS.jl/src/FastPolynomials/src/arithmetic.jl` (lines 211-215)

**Function Signature**:
```julia
_neat_dot3(x::Monomial, y::Monomial, z::Monomial) -> Monomial
```

**Implementation**:
```julia
function _neat_dot3(x::Monomial, y::Monomial, z::Monomial)
    return Monomial(
        _concat_var_expos3(reverse(x.vars), reverse(x.z), y.vars, y.z, z.vars, z.z)...
    )
end
```

**Semantic Meaning**:
```julia
_neat_dot3(x, y, z) == star(x) * y * z
```

**Usage Example** (from test file):
```julia
@ncpolyvar x y z
mono1 = monomial([x, y], [1, 0])  # x
mono2 = monomial([x, y], [1, 1])  # xy
mono3 = monomial([z, x], [2, 1])  # z²x

result = _neat_dot3(mono1, mono2, mono3)
# Equivalent to: neat_dot(mono1, mono2 * mono3)
```

**Performance Note**: Optimized for three-way products, avoiding intermediate monomial allocations.

---

## 2. Usage in NCTSSoS.jl Codebase

### 2.1 Moment Matrix Construction

**File**: `/Users/yushengzhao/projects/NCTSSoS.jl/src/moment_solver.jl` (lines 100-101)

```julia
moment_mtx = [
    sum([T_prom(coef) * monomap[simplify!(expval(_neat_dot3(row_idx, mono, col_idx)), sa)]
        for (coef, mono) in zip(coefficients(poly), monomials(poly))])
    for row_idx in local_basis, col_idx in local_basis
]
```

**Pattern**: Computes `row_idx^† * mono * col_idx` for moment matrix entries.

---

### 2.2 GNS Reconstruction

**File**: `/Users/yushengzhao/projects/NCTSSoS.jl/src/gns.jl` (lines 180, 237)

```julia
# Hankel dictionary construction
key = neat_dot(row_mono, col_mono)
dict[key] = hankel[i, j]

# Localizing matrix construction
key = neat_dot(row_mono, neat_dot(monomial(var), col_mono))
K[i, j] = get(hankel_dict, key, zero(T))
```

**Pattern**: Uses `neat_dot` to compute inner product keys for dictionary lookups.

---

### 2.3 Term Sparsity Graph

**File**: `/Users/yushengzhao/projects/NCTSSoS.jl/src/sparse.jl` (lines 167, 214-215)

```julia
# Initial support includes adjoint products of basis elements
[simplify!(neat_dot(b, b), sa) for b in mom_mtx_bases]

# Term sparsity graph connections
connected_mono_lr = _neat_dot3(bases[i], supp, bases[j])
connected_mono_rl = _neat_dot3(bases[j], supp, bases[i])
```

**Pattern**: Both left-right and right-left products needed for graph connectivity.

---

## 3. Implementation Recommendations

### 3.1 For `adjoint_product(mono1, mono2)` Function

**Recommendation**: Use `neat_dot` directly as it already implements the required operation.

```julia
"""
    adjoint_product(mono1::M, mono2::M) where {M<:Monomial}

Compute mono1^† * mono2 (adjoint of mono1 times mono2).

This is equivalent to reversing mono1's variable order and exponents,
then multiplying by mono2. Used extensively in moment matrix construction
where entries are indexed by basis products.

# Arguments
- `mono1::M`: First monomial (will be adjointed)
- `mono2::M`: Second monomial

# Returns
- `M`: The product star(mono1) * mono2

# Examples
```jldoctest
julia> @ncpolyvar x y
(x, y)

julia> mono1 = monomial([x, y], [1, 1])  # xy
xy¹

julia> mono2 = monomial([x], [2])  # x²
x²

julia> adjoint_product(mono1, mono2)
yx²

julia> # This is star(xy) * x² = yx * x² = yx²
```

# Implementation Note
This function is a semantic wrapper around FastPolynomials.neat_dot,
which efficiently computes the adjoint product without intermediate allocations.
"""
adjoint_product(mono1::M, mono2::M) where {M<:Monomial} = neat_dot(mono1, mono2)
```

**Rationale**:
1. **Already implemented**: `neat_dot` does exactly what's needed
2. **Performance optimized**: Avoids intermediate allocation of `star(mono1)`
3. **Well-tested**: Used throughout the codebase
4. **Semantic clarity**: The wrapper makes intent explicit in moment matrix code

---

### 3.2 Alternative: Direct Implementation

If you prefer not to depend on `neat_dot`:

```julia
function adjoint_product(mono1::Monomial, mono2::Monomial)
    return Monomial(
        _concat_var_expos(
            reverse(mono1.vars),
            reverse(mono1.z),
            mono2.vars,
            mono2.z
        )...
    )
end
```

**Note**: This is literally the implementation of `neat_dot`, so using the existing function is preferred.

---

## 4. Important Caveats and Notes

### 4.1 Hermitian Variables

All variables in NCTSSoS.jl are assumed to be **Hermitian** (self-adjoint):
```julia
star(x) == x  # for all variables x
```

This means:
- `star(x²y) = yx²` (reverse order, but variables themselves unchanged)
- Complex conjugation only affects coefficients in polynomials
- The `star` operation only reverses the **order** of variables

### 4.2 Monomial Consolidation

The multiplication functions automatically consolidate consecutive identical variables:

```julia
@ncpolyvar x y
mono1 = monomial([x, y], [1, 0])  # Becomes just x (y⁰ removed)
mono2 = monomial([x], [2])        # x²

mono1 * mono2  # Results in x³, not x¹x²
```

**Key Point**: `_concat_var_expos` handles this automatically when the last variable of `mono1` equals the first variable of `mono2`.

### 4.3 Simplification Required

After computing adjoint products, you typically need to simplify the result:

```julia
result = neat_dot(mono1, mono2)
simplified = simplify!(result, sa)  # sa is SimplifyAlgorithm
canonical = canonicalize(simplified, sa)
```

**From codebase usage**:
```julia
simplify!(expval(_neat_dot3(row_idx, mono, col_idx)), sa)
```

### 4.4 Type Parameters

All monomial operations preserve the monomial type:
- `Monomial` → `Monomial`
- No type conversions or promotions
- Works with any monomial-like type that implements the interface

### 4.5 Performance Considerations

**Optimal performance ordering** (fastest to slowest):
1. `_neat_dot3(x, y, z)` - single allocation for three-way product
2. `neat_dot(x, y)` - single allocation for two-way product
3. `star(x) * y` - two allocations (star creates copy, then multiply)

**For moment matrix construction**:
```julia
# GOOD: Direct use of _neat_dot3
_neat_dot3(basis[i], mono, basis[j])

# SUBOPTIMAL: Multiple intermediate allocations
star(basis[i]) * mono * basis[j]
```

---

## 5. Integration with SimplifyAlgorithm

Many monomial operations need to be followed by simplification:

### 5.1 SimplifyAlgorithm Structure

**Location**: Used throughout codebase, definition needs investigation

**Key operations**:
```julia
simplify(mono::Monomial, sa::SimplifyAlgorithm) -> Monomial
simplify!(mono::Monomial, sa::SimplifyAlgorithm) -> Monomial  # in-place
canonicalize(mono::Monomial, sa::SimplifyAlgorithm) -> Monomial
```

**Purpose**:
- Apply commutation relations
- Handle projective/unipotent constraints
- Ensure canonical ordering

### 5.2 Usage Pattern

From `src/moment_solver.jl`:
```julia
simplify!(expval(_neat_dot3(row_idx, mono, col_idx)), sa)
```

From `src/sparse.jl`:
```julia
simplify!(neat_dot(b, b), sa)
```

**Recommendation**: Always wrap adjoint products in `simplify!()` when building support structures.

---

## 6. Recommended Implementation Pattern

For moment matrix extraction, follow this pattern:

```julia
# 1. Define wrapper for semantic clarity (optional but recommended)
adjoint_product(m1::M, m2::M) where {M<:Monomial} = neat_dot(m1, m2)

# 2. Use in support construction
function build_moment_support(...)
    for i in 1:length(basis), j in i:length(basis)
        # Compute basis product: basis[i]^† * basis[j]
        product = adjoint_product(basis[i], basis[j])

        # Simplify to canonical form
        canonical_product = simplify!(product, sa)

        # Look up or store in global support
        # ...
    end
end

# 3. For moment matrix entries with polynomial
function compute_moment_entry(row_basis, poly_mono, col_basis, sa)
    # This computes row_basis^† * poly_mono * col_basis
    product = _neat_dot3(row_basis, poly_mono, col_basis)
    return simplify!(expval(product), sa)
end
```

---

## 7. Code Examples from Codebase

### 7.1 Complete Example: GNS Localizing Matrix

**File**: `src/gns.jl` (lines 227-242)

```julia
function construct_localizing_matrix(
    hankel_dict::Dict{Monomial,T},
    var::Variable,
    basis::Vector{Monomial},
) where {T<:Number}
    n = length(basis)
    K = zeros(T, n, n)

    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            # Compute row_mono^† * var * col_mono
            key = neat_dot(row_mono, neat_dot(monomial(var), col_mono))
            K[i, j] = get(hankel_dict, key, zero(T))
        end
    end

    return K
end
```

**Pattern**: Nested `neat_dot` calls compute `row^† * (var * col)`.

---

### 7.2 Complete Example: Moment Matrix Constraint

**File**: `src/moment_solver.jl` (lines 91-105)

```julia
function constrain_moment_matrix!(
    model::GenericModel{T1},
    poly::P,
    local_basis::Vector{M1},
    monomap::Dict{M2,JS},
    cone,
    sa::SimplifyAlgorithm
) where {T,T1,P<:AbstractPolynomial{T},M1,M2,JS<:AbstractJuMPScalar}
    T_prom = promote_type(T, T1)
    moment_mtx = [
        sum([T_prom(coef) * monomap[simplify!(expval(_neat_dot3(row_idx, mono, col_idx)), sa)]
            for (coef, mono) in zip(coefficients(poly), monomials(poly))])
        for row_idx in local_basis, col_idx in local_basis
    ]
    return @constraint(model, moment_mtx in cone)
end
```

**Pattern**: For each (row, col) pair and each monomial in polynomial, compute `row^† * mono * col`, simplify, and look up in monomial map.

---

## 8. Summary Table

| Operation | Function | Semantics | Performance | Use Case |
|-----------|----------|-----------|-------------|----------|
| Adjoint | `star(m)` | Reverse vars/expos | Medium (copy) | Standalone adjoint |
| Multiply | `m1 * m2` | Concatenate | Fast | Direct products |
| Adjoint product | `neat_dot(m1, m2)` | `star(m1) * m2` | Fast | Moment matrix indexing |
| 3-way adjoint product | `_neat_dot3(m1,m2,m3)` | `star(m1)*m2*m3` | Fastest | Localizing matrices |

---

## 9. Final Recommendations

### For Moment Matrix Extraction Implementation

1. **Use `neat_dot` directly** for adjoint products
   - Optionally wrap in `adjoint_product()` for semantic clarity
   - Already exported from FastPolynomials module

2. **Always simplify results** when building support structures
   - Use `simplify!(product, sa)` after computing products
   - Ensures canonical form for dictionary lookups

3. **Use `_neat_dot3` for efficiency** when computing moment matrix entries
   - Pattern: `_neat_dot3(row_basis, poly_mono, col_basis)`
   - Avoids intermediate allocations

4. **Follow existing patterns** from `src/gns.jl` and `src/moment_solver.jl`
   - These files show correct usage throughout

5. **Import from FastPolynomials**:
   ```julia
   using ..FastPolynomials: neat_dot, _neat_dot3, simplify!, canonicalize
   ```

---

## 10. Open Questions for Implementation

### Resolved:
- ✅ How to compute monomial adjoint: `star(mono)`
- ✅ How to multiply monomials: `mono1 * mono2`
- ✅ How to compute adjoint product: `neat_dot(mono1, mono2)`
- ✅ Existing helper functions: `neat_dot`, `_neat_dot3` are perfect

### Remaining:
- ⚠️  Where is `SimplifyAlgorithm` stored in the problem flow?
  - Appears in `MomentProblem.sa` field
  - Appears as parameter to `moment_relax()` and `sos_dualize()`
  - Need to verify exact structure during implementation

---

## Appendices

### A. Related Files

**Core FastPolynomials Files**:
- `src/FastPolynomials/src/monomials.jl` - Monomial type and `star()`
- `src/FastPolynomials/src/arithmetic.jl` - Multiplication and `neat_dot()`
- `src/FastPolynomials/src/compare.jl` - Ordering for sorting/searching
- `src/FastPolynomials/src/variables.jl` - Variable type and basis generation

**NCTSSoS.jl Usage**:
- `src/gns.jl` - GNS reconstruction using `neat_dot`
- `src/moment_solver.jl` - Moment matrix construction using `_neat_dot3`
- `src/sparse.jl` - Sparsity graph construction using adjoint products

**Tests**:
- `test/fastpoly_test/monomials.jl` - Tests for `star` and `neat_dot`
- `test/fastpoly_test/arithmetic.jl` - Tests for monomial multiplication

### B. Key Exports

From `src/NCTSSoS.jl`:
```julia
using .FastPolynomials:
    sorted_union, monomials, sorted_unique, maxdegree, get_basis,
    neat_dot, _neat_dot3, monomials, coefficients, terms, expval
```

These functions are already available for use in the codebase.

---

**End of Research Report**
