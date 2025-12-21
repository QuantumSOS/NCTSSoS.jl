# Migration Plan: State/Trace Polynomial Tests

## Executive Summary

- **Test Files**: `test/state_poly_opt.jl`, `test/trace_poly_opt.jl`
- **Status**: BLOCKED - Tests migrated syntactically but solver pipeline doesn't support StatePolyOpt
- **Root Cause**: `correlative_sparsity()` and `moment_relax()` only accept `Polynomial`, not `NCStatePolynomial`
- **Estimated Complexity**: HIGH (requires significant solver infrastructure extension)
- **Blockers**: Multiple solver methods need StatePolyOpt support (see below)

## Current Status (Updated 2024)

### What's Done
1. ✅ `get_state_basis()` implemented in `src/FastPolynomials/src/state_word.jl`
2. ✅ `StatePolyOpt` type added to `src/pop.jl`
3. ✅ `polyopt(::NCStatePolynomial, ...)` method added
4. ✅ Test files migrated to new variable API (`create_unipotent_variables`, etc.)
5. ✅ Test files use `one(typeof(x[1]))` for correct algebra-typed identity

### What's Blocking

**The solver pipeline (`cs_nctssos`) calls functions that only accept `Polynomial`, not `NCStatePolynomial`:**

```
cs_nctssos(spop::StatePolyOpt, ...)
    └── correlative_sparsity(pop, order, elim_algo)  ❌ NO METHOD
            └── get_correlative_graph(registry, obj, cons)  ❌ uses variable_indices(::Polynomial)
            └── assign_constraint(cliques, cons, registry)  ❌ uses variable_indices(::Polynomial)
    └── moment_relax(pop, corr_sparsity, term_sparsities)  ❌ NO METHOD
            └── _build_constraint_matrix(poly, basis, cone)  ❌ expects Polynomial
```

### Missing Infrastructure

1. **`variable_indices(::NCStatePolynomial)`** - Extract variable indices from state polynomials
2. **`correlative_sparsity(::StatePolyOpt, ...)`** - Build sparsity graph from state poly
3. **`moment_relax(::StatePolyOpt, ...)`** - Build moment matrices for state words
4. **`maxdegree(::NCStatePolynomial)`** - May need implementation
5. **State-aware term sparsity** - `term_sparsities()` needs adaptation

### Error When Running

```julia
ERROR: MethodError: no method matching correlative_sparsity(
  ::StatePolyOpt{UnipotentAlgebra, Arbitrary, NCStatePolynomial{...}}, 
  ::Int64, 
  ::NoElimination
)
```

## Decision: Defer Full Implementation

Given the scope (multiple solver methods need rewrite), this migration is deferred.

**Committed WIP includes:**
- `StatePolyOpt` type (problem definition works)
- `get_state_basis()` function  
- Test files with new API syntax (disabled)

**Future work required:**
- Extend solver pipeline for NCStatePolynomial support
- OR: Design alternative approach (e.g., convert NCStatePolynomial → Polynomial internally)

## Background

### What These Tests Do

1. **State Polynomial Tests (`state_poly_opt.jl`)**
   - Test `ς()` operator (creates `StateWord{Arbitrary}` from monomials)
   - Optimize expressions like `-1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2])...`
   - Use `UnipotentAlgebra` (`is_unipotent=true` in old API)
   - Test sparse optimization with `MMD()` and `MaximalElimination()`

2. **Trace Polynomial Tests (`trace_poly_opt.jl`)**
   - Test `tr()` operator (creates `StateWord{MaxEntangled}` from monomials)
   - Optimize trace expressions like `tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])`
   - Use `ProjectorAlgebra` (`is_projective=true`) and `UnipotentAlgebra`

### Old API Patterns (to migrate FROM)

```julia
# Variable creation (OLD)
@ncpolyvar x[1:2] y[1:2]

# Problem setup (OLD)
spop = polyopt(sp * one(Monomial); is_unipotent=true, comm_gps=[x, y])

# State basis (OLD - MISSING FUNCTION)
sa = SimplifyAlgorithm(comm_gps=[x], is_unipotent=false, is_projective=false)
basis = get_state_basis(Arbitrary, x, 1, sa)
```

### New API Patterns (to migrate TO)

```julia
# Variable creation (NEW)
reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])

# Problem setup (NEW)
spop = polyopt(sp * one(Monomial), reg)

# State basis (NEW - NEEDS IMPLEMENTATION)
basis = get_state_basis(reg, order)  # Returns Vector{NCStateWord}
```

## Detailed Analysis

### Test File: `state_poly_opt.jl`

#### Test 1: "State Polynomial Opt 7.2.0" (Lines 24-51)
**Purpose**: Basic CHSH-like Bell inequality optimization

```julia
# OLD CODE
@ncpolyvar x[1:2] y[1:2]
sp = -1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) + 1.0 * ς(x[2] * y[2])
spop = polyopt(sp * one(Monomial); is_unipotent=true, comm_gps=[x, y])
```

**Migration**:
```julia
# NEW CODE
reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
sp = -1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) + 1.0 * ς(x[2] * y[2])
spop = polyopt(sp * one(Monomial), reg)
```

**Expected Result**: `-2.8284271321623202` (Tsirelson bound)

#### Test 2: "State Polynomial Opt 7.2.1" (Lines 56-76)
**Purpose**: Higher-order state polynomial optimization

Similar pattern - multiply state polynomials and optimize.

#### Test 3: "State Polynomial Opt 7.2.2" (Lines 85-100)
**Purpose**: Covariance-based optimization

Uses `cov(a, b) = ς(x[a] * y[b]) - ς(x[a]) * ς(y[b])` pattern.

#### Test 4: "Constrain Moment matrix" (Lines 106-149)
**Purpose**: Test internal moment matrix construction

**CRITICAL**: Uses `get_state_basis()` which doesn't exist:
```julia
sa = SimplifyAlgorithm(comm_gps=[x], is_unipotent=false, is_projective=false)
basis = get_state_basis(Arbitrary, x, 1, sa)  # MISSING FUNCTION
```

### Test File: `trace_poly_opt.jl`

#### Test 1: "Example 6.1" (Lines 14-34)
**Purpose**: Basic trace polynomial with projectors

```julia
# OLD CODE
@ncpolyvar x[1:3]
p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * one(Monomial)
spop = polyopt(p; is_projective=true, comm_gps=[x])
```

**Migration**:
```julia
# NEW CODE
reg, (x,) = create_projector_variables([("x", 1:3)])
p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * one(Monomial)
tpop = polyopt(p, reg)
```

#### Test 2: "Example 6.2.0" (Lines 39-52)
**Purpose**: Trace with unipotent variables

#### Test 3: "Example 6.2.1" (Lines 56-69)
**Purpose**: Squared trace expressions

#### Test 4: "Example 6.2.2" (Lines 75-94)
**Purpose**: Trace covariance patterns

## API Changes Summary

| Old API | New API | Notes |
|---------|---------|-------|
| `@ncpolyvar x[1:2]` | `reg, (x,) = create_unipotent_variables([("x", 1:2)])` | Registry-based |
| `polyopt(p; is_unipotent=true)` | `polyopt(p, reg)` | Algebra type from registry |
| `polyopt(p; is_projective=true)` | Use `create_projector_variables` | Algebra type inference |
| `polyopt(p; comm_gps=[x, y])` | Removed | Commutation handled by algebra |
| `SimplifyAlgorithm(...)` | Removed | Simplification via algebra type |
| `get_state_basis(ST, vars, d, sa)` | `get_state_basis(reg, d)` | **NEEDS IMPLEMENTATION** |

## Missing Components

### 1. `get_state_basis()` Function (BLOCKER)

**Location**: Should be in `src/FastPolynomials/src/basis.jl` or new file

**Required Signature**:
```julia
function get_state_basis(registry::VariableRegistry{A,T}, d::Int; 
                         state_type::Type{ST}=Arbitrary) where {A<:AlgebraType, T<:Integer, ST<:StateType}
    # Generate basis of NCStateWord elements up to degree d
    # Returns Vector{NCStateWord{ST,A,T}}
end
```

**Implementation Approach**:
1. Get monomial basis using existing `get_ncbasis(registry, d)`
2. Convert each monomial to `NCStateWord` with identity StateWord
3. Return sorted, unique NCStateWord vector

**Draft Implementation**:
```julia
function get_state_basis(registry::VariableRegistry{A,T}, d::Int;
                         state_type::Type{ST}=Arbitrary) where {A<:AlgebraType, T<:Integer, ST<:StateType}
    # Get polynomial basis (handles simplification)
    poly_basis = get_ncbasis(registry, d)
    
    # Convert to NCStateWord basis
    result = NCStateWord{ST,A,T}[]
    identity_sw = one(StateWord{ST,A,T})
    
    for poly in poly_basis
        for term in poly.terms
            # Each term.monomial becomes an NCStateWord
            ncsw = NCStateWord(identity_sw, term.monomial)
            push!(result, ncsw)
        end
    end
    
    # Sort and deduplicate
    unique!(sort!(result))
    return result
end
```

### 2. Type Compatibility

The existing `ς()` and `tr()` functions in `state_word.jl` (lines 668 and 695) work correctly:
```julia
ς(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer} = StateWord{Arbitrary}(m)
tr(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer} = StateWord{MaxEntangled}(m)
```

The issue is that the polynomial-level `ς(p::Polynomial)` conversion in `utils.jl` (line 177) returns `StatePolynomial`, which then needs to multiply with `one(Monomial)` to create `NCStatePolynomial` for optimization.

## Migration Steps

### Step 1: Implement `get_state_basis()` (Required)

1. Add function to `src/FastPolynomials/src/basis.jl`
2. Export from `FastPolynomials.jl`
3. Import in `NCTSSoS.jl`

### Step 2: Migrate `state_poly_opt.jl`

| Line | Change |
|------|--------|
| 1 | Keep imports |
| 20 | Remove `get_state_basis` from import (will be different signature) |
| 25 | Replace `@ncpolyvar x[1:2] y[1:2]` → `reg, (x, y) = create_unipotent_variables(...)` |
| 29 | Replace `polyopt(...; is_unipotent=true, comm_gps=[x, y])` → `polyopt(..., reg)` |
| 57 | Same pattern as line 25 |
| 62 | Same pattern as line 29 |
| 93-98 | Uncomment and apply same patterns |
| 107-149 | Update to use new `get_state_basis(reg, d)` signature |

### Step 3: Migrate `trace_poly_opt.jl`

| Line | Change |
|------|--------|
| 15 | Replace `@ncpolyvar x[1:3]` → `reg, (x,) = create_projector_variables([("x", 1:3)])` |
| 17 | Keep `tr()` usage (works with new Monomial types) |
| 19 | Replace `polyopt(...; is_projective=true, comm_gps=[x])` → `polyopt(..., reg)` |
| 40 | Use `create_unipotent_variables` for unipotent tests |
| 44 | Remove `is_unipotent=true` |
| 57-68 | Same pattern |
| 87-93 | Same pattern |

### Step 4: Verify Test Results

Expected values from original tests:
- State Poly Opt 7.2.0: `-2.8284271321623202` (atol=1e-5)
- State Poly Opt 7.2.1: `-4.0` (atol=1e-4)
- State Poly Opt 7.2.2: `-5.0` (atol=1e-2)
- Trace Example 6.1 (order 2): `-0.046717378455438933` (atol=1e-6)
- Trace Example 6.1 (order 3): `-0.03124998978001017` (atol=1e-6)
- Trace Example 6.2.0: `-2.8284271157283083` (atol=1e-5)
- Trace Example 6.2.1: `-4.000000007460838` (atol=1e-5)
- Trace Example 6.2.2: `-5.0` (atol=1e-5)

## Risk Assessment

### High Risk
- **get_state_basis() implementation**: May have subtle differences from old behavior
- **StatePolynomial × Monomial interaction**: May have type inference issues

### Medium Risk
- **Different sparsity patterns**: New algebra types may produce different clique structures
- **Numerical precision**: Solver behavior may differ slightly

### Low Risk
- **Variable creation**: Well-tested pattern from other migrations
- **polyopt() signature**: Straightforward change

## Dependencies

### Requires
- FastPolynomials basis.jl (for get_state_basis implementation)
- StateWord, NCStateWord types (already exist)
- VariableRegistry (already exists)

### No Changes Needed
- StatePolynomial arithmetic (already implemented)
- NCStatePolynomial arithmetic (already implemented)
- ς() and tr() functions (already work with new Monomial types)

## Testing Strategy

1. **Unit test `get_state_basis()`** first
2. **Enable tests incrementally** (remove `if false` one test at a time)
3. **Compare numerical results** against expected values
4. **Check for type stability** using `@code_warntype`

## Timeline Estimate

| Task | Estimate |
|------|----------|
| Implement `get_state_basis()` | 2-4 hours |
| Migrate `state_poly_opt.jl` | 1-2 hours |
| Migrate `trace_poly_opt.jl` | 1 hour |
| Debug and verify | 2-4 hours |
| **Total** | **6-11 hours** |

## Reference: Existing State/Trace Types

### StateWord (`src/FastPolynomials/src/state_word.jl`)
```julia
struct StateWord{ST<:StateType,A<:AlgebraType,T<:Integer} <: AbstractMonomial
    state_monos::Vector{Monomial{A,T}}  # Sorted expectations <M1><M2>...
    hash::UInt64
end
```

### NCStateWord (`src/FastPolynomials/src/state_word.jl`)
```julia
struct NCStateWord{ST<:StateType,A<:AlgebraType,T<:Integer}
    sw::StateWord{ST,A,T}      # Commutative state word part
    nc_word::Monomial{A,T}     # Non-commutative operator part
    hash::UInt64
end
```

### StatePolynomial (`src/FastPolynomials/src/state_polynomial.jl`)
```julia
struct StatePolynomial{C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    coeffs::Vector{C}
    state_words::Vector{StateWord{ST,A,T}}
end
```

### NCStatePolynomial (`src/FastPolynomials/src/state_polynomial.jl`)
```julia
struct NCStatePolynomial{C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    coeffs::Vector{C}
    nc_state_words::Vector{NCStateWord{ST,A,T}}
end
```

## Appendix: Full Test File Migration Templates

### state_poly_opt.jl (Migrated Template)

```julia
using Test, NCTSSoS, NCTSSoS.FastPolynomials

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

using COSMO
const QUICK_SOLVER = COSMO.Optimizer
using JuMP
using NCTSSoS.FastPolynomials: expval, terms, Arbitrary, NCStateWord, ς, Monomial

@testset "State Polynomial Opt 7.2.0" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp =
        -1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) +
        1.0 * ς(x[2] * y[2])
    spop = polyopt(sp * one(Monomial), reg)

    d = 1
    solver_config = SolverConfig(; optimizer=SOLVER, order=d)

    if haskey(ENV, "LOCAL_TESTING")
        result_mom = cs_nctssos(spop, solver_config; dualize=false)
        @test isapprox(result_mom.objective, -2.8284271321623202, atol=1e-5)
    end

    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-5)

    @testset "Sparse" begin
        solver_config = SolverConfig(; optimizer=SOLVER, order=d, cs_algo=NoElimination(), ts_algo=MMD())
        result = cs_nctssos(spop, solver_config)
        @test result.objective ≈ -2.8284271321623202 atol = 1e-5
    end
end

# ... additional tests follow same pattern
```

### trace_poly_opt.jl (Migrated Template)

```julia
using Test, NCTSSoS
using NCTSSoS.FastPolynomials: tr, Monomial

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

@testset "Example 6.1" begin
    reg, (x,) = create_projector_variables([("x", 1:3)])

    p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * one(Monomial)

    tpop = polyopt(p, reg)

    solver_config = SolverConfig(; optimizer=SOLVER, order=2)

    result = cs_nctssos(tpop, solver_config)

    @test result.objective ≈ -0.046717378455438933 atol = 1e-6

    if haskey(ENV, "LOCAL_TESTING")
        solver_config = SolverConfig(; optimizer=SOLVER, order=3)
        result = cs_nctssos(tpop, solver_config)
        @test result.objective ≈ -0.03124998978001017 atol = 1e-6
    end
end

# ... additional tests follow same pattern
```
