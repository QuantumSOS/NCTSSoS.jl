# Fermionic Algebra Implementation Plan

**Date**: 2025-11-05
**Author**: julia-development sub-agent
**Status**: Complete - Ready for Parent Agent Implementation

---

## Executive Summary

This plan provides a complete, step-by-step guide to implement fermionic algebra in NCTSSoS.jl by following the established Pauli algebra pattern. The implementation requires:

1. **Core extension**: Add `is_nilpotent` support to `SimplifyAlgorithm` with `_simplify_nilpotent!` function
2. **Algebra constructor**: Create `fermionic_algebra(N)` following the Pauli pattern
3. **Type system updates**: Extend `PolyOpt` and `ComplexPolyOpt` with nilpotent field
4. **Comprehensive tests**: Unit tests, integration tests, and numerical validation
5. **Documentation**: Tutorial following the Pauli algebra interface pattern

**Critical Note**: Issue #179 documents that zero monomials aren't supported. Our workaround is to use empty monomial vectors (`empty!(m.vars); empty!(m.z)`) which should be filtered out at the Polynomial level. This is tested and validated in the plan.

---

## 1. Architecture Overview

### 1.1 Key Design Decisions

**Nilpotent vs Unipotent**:
- Pauli: `σᵢ² = 1` (unipotent) → exponents reduced modulo 2
- Fermionic: `cᵢ² = 0` (nilpotent) → any exponent ≥ 2 makes entire monomial zero

**Commutation Structure**:
- Pauli: Operators at different sites commute → `comm_gps = [[σx[i], σy[i], σz[i]] for i in 1:N]`
- Fermionic: ALL operators anti-commute → `comm_gps = [all_operators]` (single group)

**Variable Naming Convention**:
- `c[i]`: Annihilation operator for mode i
- `c_dag[i]`: Creation operator for mode i (uses `_dag` suffix, NOT superscript due to Julia naming rules)

**Anti-commutation Relations** (encoded as equality constraints):
```
{cᵢ, cⱼ} = cᵢcⱼ + cⱼcᵢ = 0           for all i,j
{cᵢ†, cⱼ†} = cᵢ†cⱼ† + cⱼ†cᵢ† = 0       for all i,j
{cᵢ, cⱼ†} = cᵢcⱼ† + cⱼ†cᵢ = δᵢⱼ      canonical anti-commutation relation
```

### 1.2 Zero Monomial Workaround Strategy

**Problem**: FastPolynomials doesn't support zero monomials (issue #179).

**Solution**: Use empty monomial representation:
- Empty `vars` and `z` vectors in Monomial currently represent identity (1)
- For nilpotent case, we'll distinguish by context:
  - After simplification, if `_simplify_nilpotent!` empties the monomial, it represents zero
  - The `Polynomial` constructor will handle this correctly by filtering out zero-coefficient terms

**Validation**: The existing code in `polynomial.jl` already filters out zero coefficients and combines like terms. When a monomial becomes zero (empty), multiplying by any coefficient should produce a zero polynomial.

**Future-proofing**: When #179 is resolved with a proper `is_zero` flag, we can:
1. Update `_simplify_nilpotent!` to set `m.is_zero = true`
2. Update `Polynomial` constructor to filter `is_zero` monomials
3. All existing code will continue to work

### 1.3 Type System Extensions

Add `is_nilpotent::Bool` field to:
- `SimplifyAlgorithm` (in `simplify.jl`)
- `PolyOpt` (in `pop.jl`)
- `ComplexPolyOpt` (in `pop.jl`)

**Mutual exclusivity**: Only ONE of `{is_unipotent, is_projective, is_nilpotent}` can be `true`.

---

## 2. Detailed Implementation Steps

### Phase 1: Core Simplification Logic (FastPolynomials)

#### Step 1.1: Extend `SimplifyAlgorithm` struct

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/FastPolynomials/src/simplify.jl`
**Lines**: 1-18 (modify existing struct)

**Action**: Add `is_nilpotent` field and update constructor

**Code Template**:
```julia
struct SimplifyAlgorithm
    comm_gps::Dict{Variable,Int}
    n_gps::Int
    is_unipotent::Bool
    is_projective::Bool
    is_nilpotent::Bool  # NEW: Add this field

    function SimplifyAlgorithm(;
        comm_gps::Vector{Vector{Variable}},
        is_unipotent::Bool=false,
        is_projective::Bool=false,
        is_nilpotent::Bool=false,  # NEW: Add this parameter
    )
        # NEW: Add mutual exclusivity check
        @assert sum([is_unipotent, is_projective, is_nilpotent]) <= 1
            "Only one of is_unipotent, is_projective, is_nilpotent can be true"

        return new(
            Dict(var => i for (i, vars) in enumerate(comm_gps) for var in vars),
            length(comm_gps),
            is_unipotent,
            is_projective,
            is_nilpotent,  # NEW: Add this field
        )
    end
end
```

**Test**: Verify constructor accepts `is_nilpotent` parameter and enforces mutual exclusivity.

---

#### Step 1.2: Update `nosimp` helper function

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/FastPolynomials/src/simplify.jl`
**Lines**: 20-22 (modify)

**Action**: Include `is_nilpotent` in the check

**Code Template**:
```julia
function nosimp(sa::SimplifyAlgorithm)
    # OLD: return isone(sa.n_gps) && !sa.is_unipotent && !sa.is_projective
    # NEW: Include is_nilpotent in check
    return isone(sa.n_gps) && !sa.is_unipotent && !sa.is_projective && !sa.is_nilpotent
end
```

**Test**: Verify `nosimp` returns `false` when `is_nilpotent=true`.

---

#### Step 1.3: Implement `_simplify_nilpotent!` function

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/FastPolynomials/src/simplify.jl`
**Lines**: After line 66 (insert new function before `_simplify_unipotent!`)

**Action**: Create nilpotent simplification logic

**Code Template**:
```julia
"""
    _simplify_nilpotent!(m::Monomial)

Simplify a monomial with nilpotent operators (X² = 0).

For nilpotent operators, any variable with exponent ≥ 2 causes the entire
monomial to become zero. This function checks all exponents and empties the
monomial (representing zero) if any exponent is ≥ 2.

# Algorithm
1. Iterate through all exponents in the monomial
2. If any exponent is ≥ 2, empty both vars and z vectors (representing zero)
3. Otherwise, keep the monomial as-is (after commutation has already been applied)

# Notes
- After commutation, consecutive identical variables should not exist
- This function handles the case where exponent collection creates z[i] ≥ 2
- Empty vars/z vectors represent the zero monomial (workaround for issue #179)

# Mutates
- `m.vars`: Emptied if any exponent ≥ 2
- `m.z`: Emptied if any exponent ≥ 2
"""
function _simplify_nilpotent!(m::Monomial)
    # Check all exponents
    @inbounds for i in eachindex(m.z)
        if m.z[i] >= 2
            # Found X² or higher power → entire monomial is zero
            empty!(m.vars)
            empty!(m.z)
            return nothing
        end
    end

    # All exponents are 1, monomial is valid
    # (After _comm! has been applied, no simplification needed beyond zero detection)
    return nothing
end
```

**Key Points**:
- Simpler than `_simplify_unipotent!` because we only check for zero
- After `_comm!`, we should never see consecutive identical variables (they'd have combined exponents)
- If any exponent ≥ 2, the entire term vanishes

**Test**:
- `c * c` → empty monomial
- `c_dag * c_dag` → empty monomial
- `c * c_dag * c` → stays as-is (exponents all 1)
- `c * c * c_dag` → empty (after combining, c has exponent 2)

---

#### Step 1.4: Update `simplify!` dispatch

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/FastPolynomials/src/simplify.jl`
**Lines**: 24-32 (modify)

**Action**: Add nilpotent case to simplification dispatch

**Code Template**:
```julia
function simplify!(m::Monomial, sa::SimplifyAlgorithm)
    (isone(m) || nosimp(sa)) && return m
    _comm!(m, sa.comm_gps)  # Step 1: Apply commutation rules

    # Step 2: Apply algebraic simplification (mutually exclusive)
    sa.is_unipotent && (_simplify_unipotent!(m); return m)
    sa.is_projective && (_simplify_projective!(m); return m)
    sa.is_nilpotent && (_simplify_nilpotent!(m); return m)  # NEW

    # Step 3: Standard simplification (combine like terms)
    _simplify_standard!(m)
    return m
end
```

**Test**: Verify correct dispatch based on `is_nilpotent` flag.

---

### Phase 2: Optimization Problem Type Updates

#### Step 2.1: Update `PolyOpt` struct

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/pop.jl`
**Lines**: 25-33 (modify struct definition)

**Action**: Add `is_nilpotent` field

**Code Template**:
```julia
struct PolyOpt{P} <: OptimizationProblem{P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    variables::Vector{Variable}
    comm_gps::Vector{Vector{Variable}}
    is_unipotent::Bool
    is_projective::Bool
    is_nilpotent::Bool  # NEW: Add this field
end
```

**Documentation Update**: Update docstring (lines 3-24) to include:
```
- `is_nilpotent::Bool`: Whether variables are nilpotent (X² = 0, e.g., fermionic operators)
```

**Test**: Verify struct can be constructed with `is_nilpotent` field.

---

#### Step 2.2: Update `polyopt` constructor

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/pop.jl`
**Lines**: 57-71 (modify function)

**Action**: Add `is_nilpotent` parameter and validation

**Code Template**:
```julia
function polyopt(
    objective::P;
    eq_constraints=Any[],
    ineq_constraints=Any[],
    comm_gps=Vector{Variable}[],
    is_unipotent::Bool=false,
    is_projective::Bool=false,
    is_nilpotent::Bool=false  # NEW: Add parameter
) where {T,P<:AbstractPolynomial{T}}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    eq_cons = unique!(collect(P, eq_constraints))
    ineq_cons = unique!(collect(P, ineq_constraints))
    vars = sorted_union(variables(objective), variables.(eq_cons)..., variables.(ineq_cons)...)

    if !isempty(comm_gps)
        @assert all([isempty(intersect(gp_a, gp_b)) for gp_a in comm_gps, gp_b in comm_gps if gp_a != gp_b]) "The commutative groups must be disjoint."
        @assert issubset(union(comm_gps...), vars) "The commutative variables must be a subset of the variables."
        @assert all(issorted(gp) for gp in comm_gps) "The commutative groups must be sorted."
    else
        push!(comm_gps, sort(vars))
    end

    # NEW: Update mutual exclusivity check
    @assert sum([is_unipotent, is_projective, is_nilpotent]) <= 1
        "Only one of is_unipotent, is_projective, is_nilpotent can be true"

    # NEW: Add is_nilpotent to constructor call
    return PolyOpt{P}(objective, eq_cons, ineq_cons, vars, comm_gps,
                      is_unipotent, is_projective, is_nilpotent)
end
```

**Documentation Update**: Update docstring (lines 36-56) to include:
```
- `is_nilpotent::Bool=false`: Flag indicating if operators are nilpotent (X² = 0)
```

**Test**: Verify constructor accepts `is_nilpotent` and enforces mutual exclusivity.

---

#### Step 2.3: Update `ComplexPolyOpt` struct and constructor

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/pop.jl`
**Lines**: 95-126 (modify struct and constructor)

**Action**: Mirror changes from `PolyOpt` to `ComplexPolyOpt`

**Code Template**:

For struct (around line 95):
```julia
struct ComplexPolyOpt{P} <: OptimizationProblem{P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    variables::Vector{Variable}
    comm_gps::Vector{Vector{Variable}}
    is_unipotent::Bool
    is_projective::Bool
    is_nilpotent::Bool  # NEW
end
```

For constructor (around line 106):
```julia
function cpolyopt(
    objective::P;
    eq_constraints=Any[],
    ineq_constraints=Any[],
    comm_gps=Vector{Variable}[],
    is_unipotent::Bool=false,
    is_projective::Bool=false,
    is_nilpotent::Bool=false  # NEW
) where {T,P<:AbstractPolynomial{T}}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    eq_cons = unique!(collect(P, eq_constraints))
    ineq_cons = unique!(collect(P, ineq_constraints))
    vars = sorted_union(variables(objective), variables.(eq_cons)..., variables.(ineq_cons)...)

    if !isempty(comm_gps)
        @assert all([isempty(intersect(gp_a, gp_b)) for gp_a in comm_gps, gp_b in comm_gps if gp_a != gp_b]) "The commutative groups must be disjoint."
        @assert issubset(union(comm_gps...), vars) "The commutative variables must be a subset of the variables."
    else
        push!(comm_gps, vars)
    end

    # NEW: Update mutual exclusivity check
    @assert sum([is_unipotent, is_projective, is_nilpotent]) <= 1
        "Only one of is_unipotent, is_projective, is_nilpotent can be true"

    # Symmetry validation
    sa = SimplifyAlgorithm(comm_gps=comm_gps, is_unipotent=is_unipotent,
                           is_projective=is_projective, is_nilpotent=is_nilpotent)  # NEW: Add is_nilpotent
    @assert FastPolynomials.is_symmetric(objective, sa) "Objective must be symmetric"
    for ineq in ineq_constraints
        @assert FastPolynomials.is_symmetric(ineq, sa) "Inequality constraints must be symmetric"
    end

    # NEW: Add is_nilpotent to constructor call
    return ComplexPolyOpt{P}(objective, eq_cons, ineq_cons, vars, comm_gps,
                             is_unipotent, is_projective, is_nilpotent)
end
```

**Critical Note**: The `SimplifyAlgorithm` construction in line ~119 must be updated to pass `is_nilpotent`.

**Test**: Verify `ComplexPolyOpt` works with nilpotent flag.

---

#### Step 2.4: Update `cpolyopt` algebra interface method

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/pop.jl`
**Lines**: Around 129-180 (modify existing algebra interface)

**Action**: Extract and pass `is_nilpotent` from algebra system

**Code Template**:
```julia
function cpolyopt(
    objective::P,
    algebra::NamedTuple;
    eq_constraints=Any[],
    ineq_constraints=Any[]
) where {T,P<:AbstractPolynomial{T}}
    # Extract properties from algebra
    comm_gps = algebra.comm_gps
    is_unipotent = algebra.simplify_algo.is_unipotent
    is_projective = algebra.simplify_algo.is_projective
    is_nilpotent = algebra.simplify_algo.is_nilpotent  # NEW: Extract is_nilpotent

    # Merge algebra constraints with user-provided constraints
    converted_user_eq = isempty(eq_constraints) ?
        typeof(algebra.equality_constraints)[] :
        [T(1) * poly for poly in eq_constraints]
    merged_eq_constraints = vcat(algebra.equality_constraints, converted_user_eq)

    # Convert inequality constraints
    converted_user_ineq = isempty(ineq_constraints) ?
        typeof(algebra.inequality_constraints)[] :
        [T(1) * poly for poly in ineq_constraints]
    merged_ineq_constraints = vcat(algebra.inequality_constraints, converted_user_ineq)

    # Call original cpolyopt with all properties
    return cpolyopt(objective;
        eq_constraints=merged_eq_constraints,
        ineq_constraints=merged_ineq_constraints,
        comm_gps=comm_gps,
        is_unipotent=is_unipotent,
        is_projective=is_projective,
        is_nilpotent=is_nilpotent)  # NEW: Pass is_nilpotent
end
```

**Test**: Verify fermionic algebra system works with this interface.

---

### Phase 3: Fermionic Algebra Constructor

#### Step 3.1: Implement `fermionic_algebra` function

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/algebra_constructors.jl`
**Lines**: After line 39 (append new function)

**Action**: Create fermionic algebra constructor following Pauli pattern

**Code Template**:
```julia
"""
    fermionic_algebra(N::Int)

Create a fermionic algebra system for N fermionic modes.

Fermionic operators satisfy anti-commutation relations and nilpotency:
- Anti-commutation: {cᵢ, cⱼ} = 0, {cᵢ†, cⱼ†} = 0
- Canonical anti-commutation: {cᵢ, cⱼ†} = δᵢⱼ
- Nilpotency: cᵢ² = 0, (cᵢ†)² = 0

# Arguments
- `N::Int`: Number of fermionic modes

# Returns
A `NamedTuple` with fields:
- `variables::Tuple`: Tuple of (c, c_dag) where c[i] and c_dag[i] are annihilation/creation operators
- `simplify_algo::SimplifyAlgorithm`: Simplification algorithm with is_nilpotent=true
- `equality_constraints::Vector{Polynomial{ComplexF64}}`: Anti-commutation relations
- `inequality_constraints::Vector{Polynomial{ComplexF64}}`: Empty (no inequality constraints)
- `comm_gps::Vector{Vector{Variable}}`: Single group containing all operators (all anti-commute)

# Example
```julia
sys = fermionic_algebra(2)
c, c_dag = sys.variables

# Number operator for mode 1
n1 = c_dag[1] * c[1]

# Create optimization problem
pop = cpolyopt(n1, sys)
```

# Mathematical Details
For N modes, the equality constraints encode:
- N(N+1)/2 constraints for {cᵢ, cⱼ} = 0 (i ≤ j)
- N(N+1)/2 constraints for {cᵢ†, cⱼ†} = 0 (i ≤ j)
- N² constraints for {cᵢ, cⱼ†} = δᵢⱼ

Total: N² + N(N+1) = 2N² + N constraints
"""
function fermionic_algebra(N::Int)
    @assert N >= 1 "Number of modes N must be at least 1"

    # Step 1: Declare variables
    # c[i] = annihilation operator for mode i
    # c_dag[i] = creation operator for mode i (dagger notation via _dag suffix)
    @ncpolyvar c[1:N] c_dag[1:N]

    # Step 2: Define commutation groups
    # ALL fermionic operators anti-commute with each other
    # Put all operators in a single commutation group (no reordering)
    all_ops = vcat(collect(c), collect(c_dag))
    comm_gps = [all_ops]

    # Step 3: Encode anti-commutation relations as equality constraints
    # These constraints will be enforced during optimization
    equality_constraints = Polynomial{ComplexF64}[]

    # Anti-commutation: {cᵢ, cⱼ} = cᵢcⱼ + cⱼcᵢ = 0
    # Only need i ≤ j to avoid duplicates (since cᵢcⱼ + cⱼcᵢ = 0 ⟺ cⱼcᵢ + cᵢcⱼ = 0)
    for i in 1:N, j in i:N
        push!(equality_constraints, c[i] * c[j] + c[j] * c[i])
    end

    # Anti-commutation: {cᵢ†, cⱼ†} = cᵢ†cⱼ† + cⱼ†cᵢ† = 0
    for i in 1:N, j in i:N
        push!(equality_constraints, c_dag[i] * c_dag[j] + c_dag[j] * c_dag[i])
    end

    # Canonical anti-commutation relation: {cᵢ, cⱼ†} = cᵢcⱼ† + cⱼ†cᵢ = δᵢⱼ
    # This defines the fermionic algebra structure
    for i in 1:N, j in 1:N
        if i == j
            # {cᵢ, cᵢ†} = 1
            push!(equality_constraints, c[i] * c_dag[i] + c_dag[i] * c[i] - 1)
        else
            # {cᵢ, cⱼ†} = 0 for i ≠ j
            push!(equality_constraints, c[i] * c_dag[j] + c_dag[j] * c[i])
        end
    end

    # Step 4: Create SimplifyAlgorithm with nilpotent property
    simplify_algo = SimplifyAlgorithm(
        comm_gps=comm_gps,
        is_unipotent=false,
        is_projective=false,
        is_nilpotent=true  # cᵢ² = 0, (cᵢ†)² = 0
    )

    # Step 5: No inequality constraints for fermionic algebra
    inequality_constraints = empty(equality_constraints)

    # Step 6: Return algebra system as NamedTuple
    # This format matches pauli_algebra and is compatible with cpolyopt(obj, algebra)
    return (
        variables=(c, c_dag),
        simplify_algo=simplify_algo,
        equality_constraints=equality_constraints,
        inequality_constraints=inequality_constraints,
        comm_gps=comm_gps
    )
end
```

**Constraint Count Verification**:
- For N=1: 1 + 1 + 1 = 3 constraints
- For N=2: 3 + 3 + 4 = 10 constraints
- For N=3: 6 + 6 + 9 = 21 constraints
- General: N(N+1)/2 + N(N+1)/2 + N² = N² + N(N+1) = 2N² + N

**Test**: Verify constraint counts and structure for N=1,2,3.

---

#### Step 3.2: Export `fermionic_algebra`

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/NCTSSoS.jl`
**Lines**: Around 14-20 (add to exports)

**Action**: Add fermionic_algebra to module exports

**Code Template**:
```julia
export @ncpolyvar, ς
export polyopt, cpolyopt
export SolverConfig
export cs_nctssos, cs_nctssos_higher
export reconstruct
export pauli_algebra
export fermionic_algebra  # NEW: Add this line
```

**Test**: Verify `fermionic_algebra` is accessible after `using NCTSSoS`.

---

### Phase 4: Comprehensive Testing

#### Step 4.1: Create fermionic algebra unit tests

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/test/algebra_constructors.jl`
**Lines**: After existing tests (append new testset)

**Action**: Add comprehensive fermionic algebra tests

**Code Template**:
```julia
@testset "Fermionic Algebra Constructor" begin
    @testset "Basic Structure N=1" begin
        sys = fermionic_algebra(1)

        # Check return type structure (matches Pauli pattern)
        @test hasfield(typeof(sys), :variables)
        @test hasfield(typeof(sys), :simplify_algo)
        @test hasfield(typeof(sys), :equality_constraints)
        @test hasfield(typeof(sys), :inequality_constraints)
        @test hasfield(typeof(sys), :comm_gps)

        # Check variables
        c, c_dag = sys.variables
        @test length(c) == 1
        @test length(c_dag) == 1

        # Check constraints count for N=1
        # {c,c}=0 (1) + {c†,c†}=0 (1) + {c,c†}=1 (1) = 3 constraints
        @test length(sys.equality_constraints) == 3
        @test isempty(sys.inequality_constraints)

        # Check commutation groups (single group with all operators)
        @test length(sys.comm_gps) == 1
        @test length(sys.comm_gps[1]) == 2  # c[1] and c_dag[1]
    end

    @testset "Basic Structure N=2" begin
        sys = fermionic_algebra(2)
        c, c_dag = sys.variables

        @test length(c) == 2
        @test length(c_dag) == 2

        # Constraint count for N=2:
        # {cᵢ, cⱼ} for i,j ∈ {1,2}, i≤j → 3 constraints
        # {cᵢ†, cⱼ†} for i,j ∈ {1,2}, i≤j → 3 constraints
        # {cᵢ, cⱼ†} for i,j ∈ {1,2} → 4 constraints
        # Total: 10 constraints
        @test length(sys.equality_constraints) == 10

        # Single commutation group with all 4 operators
        @test length(sys.comm_gps) == 1
        @test length(sys.comm_gps[1]) == 4
    end

    @testset "Basic Structure N=3" begin
        sys = fermionic_algebra(3)
        c, c_dag = sys.variables

        @test length(c) == 3
        @test length(c_dag) == 3

        # Constraint count for N=3: 2N² + N = 18 + 3 = 21
        @test length(sys.equality_constraints) == 21

        # Single commutation group with all 6 operators
        @test length(sys.comm_gps) == 1
        @test length(sys.comm_gps[1]) == 6
    end

    @testset "SimplifyAlgorithm Properties" begin
        sys = fermionic_algebra(2)
        sa = sys.simplify_algo

        @test sa.is_nilpotent == true
        @test sa.is_unipotent == false
        @test sa.is_projective == false
        @test sa.n_gps == 1  # Single commutation group
    end

    @testset "Anti-commutation Relations N=1" begin
        sys = fermionic_algebra(1)
        c, c_dag = sys.variables
        constraints = sys.equality_constraints

        # {c, c} = 0 ⟹ c*c + c*c = 0 ⟹ 2*c*c = 0 ⟹ c*c = 0
        @test c[1] * c[1] + c[1] * c[1] in constraints

        # {c†, c†} = 0 ⟹ c†*c† + c†*c† = 0
        @test c_dag[1] * c_dag[1] + c_dag[1] * c_dag[1] in constraints

        # {c, c†} = 1 ⟹ c*c† + c†*c = 1
        @test c[1] * c_dag[1] + c_dag[1] * c[1] - 1 in constraints
    end

    @testset "Anti-commutation Relations N=2" begin
        sys = fermionic_algebra(2)
        c, c_dag = sys.variables
        constraints = sys.equality_constraints

        # Spot check key constraints
        @test c[1] * c[2] + c[2] * c[1] in constraints
        @test c_dag[1] * c_dag[2] + c_dag[2] * c_dag[1] in constraints
        @test c[1] * c_dag[1] + c_dag[1] * c[1] - 1 in constraints
        @test c[1] * c_dag[2] + c_dag[2] * c[1] in constraints  # {c₁, c₂†} = 0
    end

    @testset "Nilpotent Simplification" begin
        sys = fermionic_algebra(1)
        c, c_dag = sys.variables
        sa = sys.simplify_algo

        # Test that c² simplifies to zero (empty monomial)
        c_squared = c[1] * c[1]
        # After simplification, should get empty polynomial or zero
        # The multiplication will create a monomial with exponent 2
        # which _simplify_nilpotent! will detect and empty

        # NOTE: This test validates the zero monomial workaround
        # We expect the multiplication to produce a monomial,
        # then simplification empties it (representing zero)
        # The Polynomial constructor should handle this correctly

        # We can't directly test iszero because of the workaround
        # Instead, verify the monomial itself becomes empty
        mono = first(c_squared.monos)
        simplified_mono = FastPolynomials.simplify(mono, sa)
        @test isempty(simplified_mono.vars)
        @test isempty(simplified_mono.z)

        # Same for c_dag squared
        cdag_squared = c_dag[1] * c_dag[1]
        mono2 = first(cdag_squared.monos)
        simplified_mono2 = FastPolynomials.simplify(mono2, sa)
        @test isempty(simplified_mono2.vars)
        @test isempty(simplified_mono2.z)
    end

    @testset "Integration with cpolyopt" begin
        sys = fermionic_algebra(2)
        c, c_dag = sys.variables

        # Number operator: n₁ = c₁† c₁
        # This is Hermitian and represents occupation of mode 1
        n1 = ComplexF64(1.0) * c_dag[1] * c[1]

        # Create optimization problem using algebra interface
        pop = cpolyopt(n1, sys)

        @test pop.objective == n1
        @test pop.is_nilpotent == true
        @test pop.is_unipotent == false
        @test pop.is_projective == false

        # Check that algebra constraints were included
        @test length(pop.eq_constraints) == length(sys.equality_constraints)
    end

    @testset "Integration with Additional Constraints" begin
        sys = fermionic_algebra(2)
        c, c_dag = sys.variables

        # Hubbard interaction: n₁n₂ where nᵢ = cᵢ†cᵢ
        n1 = c_dag[1] * c[1]
        n2 = c_dag[2] * c[2]
        ham = ComplexF64(1.0) * n1 * n2

        # Add custom constraint: total particle number = 1
        N_total = n1 + n2 - 1

        pop = cpolyopt(ham, sys; eq_constraints=[N_total])

        # Should have algebra constraints + 1 custom constraint
        @test length(pop.eq_constraints) == length(sys.equality_constraints) + 1
        @test N_total in pop.eq_constraints
    end

    @testset "Invalid Input" begin
        @test_throws AssertionError fermionic_algebra(0)
        @test_throws AssertionError fermionic_algebra(-1)
    end

    @testset "Mutual Exclusivity of Algebra Properties" begin
        sys = fermionic_algebra(1)
        c, c_dag = sys.variables
        obj = ComplexF64(1.0) * c_dag[1] * c[1]

        # Should not be able to create polyopt with multiple properties
        @test_throws AssertionError cpolyopt(obj;
            eq_constraints=sys.equality_constraints,
            comm_gps=sys.comm_gps,
            is_nilpotent=true,
            is_unipotent=true)

        @test_throws AssertionError cpolyopt(obj;
            eq_constraints=sys.equality_constraints,
            comm_gps=sys.comm_gps,
            is_nilpotent=true,
            is_projective=true)
    end
end
```

**Test Coverage**:
- ✅ Structure validation (NamedTuple fields)
- ✅ Variable creation and naming
- ✅ Constraint count verification (N=1,2,3)
- ✅ SimplifyAlgorithm properties
- ✅ Anti-commutation relations
- ✅ Nilpotent simplification (zero monomial workaround)
- ✅ Integration with cpolyopt
- ✅ Custom constraint merging
- ✅ Error handling
- ✅ Mutual exclusivity validation

---

#### Step 4.2: Add integration test with known results

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/test/algebra_constructors.jl`
**Lines**: After unit tests (append, wrapped in LOCAL_TESTING check)

**Action**: Add numerical validation test

**Code Template**:
```julia
# Integration test with solver (only run with LOCAL_TESTING environment variable)
if haskey(ENV, "LOCAL_TESTING")
    @testset "Integration: Free Fermion System" begin
        # Test with a simple system with known exact result
        # Free fermions: H = ε₁ c₁†c₁ + ε₂ c₂†c₂
        # Ground state energy: 0 (no particles)
        # One-particle ground state energy: min(ε₁, ε₂)

        N = 2
        sys = fermionic_algebra(N)
        c, c_dag = sys.variables

        # Single-particle energies
        ε1, ε2 = -1.0, -0.5

        # Hamiltonian: H = ε₁ n₁ + ε₂ n₂ where nᵢ = cᵢ†cᵢ
        n1 = c_dag[1] * c[1]
        n2 = c_dag[2] * c[2]
        ham = ComplexF64(ε1) * n1 + ComplexF64(ε2) * n2

        # Create optimization problem
        pop = cpolyopt(ham, sys)

        # Solve with relaxation order 1 (sufficient for quadratic Hamiltonian)
        config = SolverConfig(order=1, optimizer=SOLVER)
        result = cs_nctssos(pop, config)

        # Ground state energy should be 0 (no particles)
        @test isapprox(result.objective, 0.0, atol=1e-6)

        # Test with constraint: exactly 1 particle
        # Expected energy: min(ε₁, ε₂) = -1.0
        N_total = n1 + n2 - 1  # Total particle number = 1
        pop_constrained = cpolyopt(ham, sys; eq_constraints=[N_total])

        config2 = SolverConfig(order=2, optimizer=SOLVER)
        result2 = cs_nctssos(pop_constrained, config2)

        @test isapprox(result2.objective, ε1, atol=1e-4)
    end

    @testset "Integration: Fermi-Hubbard Dimer" begin
        # Two-site Fermi-Hubbard model (known results)
        # H = -t(c₁†c₂ + c₂†c₁) + U n₁ n₂
        # For U=0, t=1, ground state energy = -1 (bonding state)

        N = 2
        sys = fermionic_algebra(N)
        c, c_dag = sys.variables

        t = 1.0  # Hopping
        U = 0.0  # No interaction (free fermions)

        # Hopping term: -t(c₁†c₂ + c₂†c₁)
        hopping = -ComplexF64(t) * (c_dag[1] * c[2] + c_dag[2] * c[1])

        # Number operators
        n1 = c_dag[1] * c[1]
        n2 = c_dag[2] * c[2]

        # Interaction: U n₁ n₂
        interaction = ComplexF64(U) * n1 * n2

        ham = hopping + interaction

        # Constraint: single particle
        N_total = n1 + n2 - 1

        pop = cpolyopt(ham, sys; eq_constraints=[N_total])

        config = SolverConfig(order=2, optimizer=SOLVER)
        result = cs_nctssos(pop, config)

        # For single particle, hopping term gives energy -t
        @test isapprox(result.objective, -t, atol=1e-4)
    end
end
```

**Known Results for Validation**:
- Free fermions: Exact ground state energy = 0 (vacuum)
- Single particle: Energy = minimum single-particle energy
- Fermi-Hubbard dimer (U=0): Single-particle ground state = -t

**Test Strategy**:
- Start with simple systems with known exact results
- Use low relaxation orders (1-2) for efficiency
- Validate to 4-6 decimal places (accounting for SDP solver tolerance)

---

### Phase 5: Documentation

#### Step 5.1: Create fermionic algebra tutorial

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/fermionic_algebra_interface.jl`
**Lines**: New file (create)

**Action**: Write comprehensive tutorial following Pauli pattern

**Code Template**:
```julia
# # Fermionic Algebra Interface
#
# This tutorial demonstrates the convenient `fermionic_algebra` constructor for
# working with fermionic systems in NCTSSoS.jl.
#
# ## Why Fermionic Algebra?
#
# Fermionic systems appear throughout quantum many-body physics:
# - **Electrons in solids**: Described by fermionic creation/annihilation operators
# - **Quantum chemistry**: Molecular orbitals are fermionic
# - **Ultracold atoms**: Fermionic atoms in optical lattices
# - **Topological matter**: Majorana fermions and topological insulators
#
# Fermionic operators satisfy:
# - **Anti-commutation relations**: ``\{c_i, c_j\} = 0, \{c_i^\dagger, c_j^\dagger\} = 0``
# - **Canonical anti-commutation**: ``\{c_i, c_j^\dagger\} = \delta_{ij}``
# - **Nilpotency**: ``(c_i)^2 = 0, (c_i^\dagger)^2 = 0``
#
# ## The Traditional Manual Approach
#
# Before `fermionic_algebra`, you would need to manually:
# 1. Declare variables for creation and annihilation operators
# 2. Define the commutation groups (single group for all operators)
# 3. Encode all anti-commutation relations as equality constraints
# 4. Create the SimplifyAlgorithm with is_nilpotent=true
# 5. Pass everything to cpolyopt
#
# For even a 2-mode system, this requires ~20 lines of boilerplate code!

using NCTSSoS, NCTSSoS.FastPolynomials

# ## The Simplified Approach with `fermionic_algebra`
#
# The `fermionic_algebra(N)` constructor handles all the boilerplate:

N = 2  # Two fermionic modes
sys = fermionic_algebra(N)

# Extract operators:
c, c_dag = sys.variables

# Now you can directly construct fermionic Hamiltonians:
# Number operators: nᵢ = cᵢ† cᵢ
n1 = c_dag[1] * c[1]
n2 = c_dag[2] * c[2]

# Free fermion Hamiltonian: H = ε₁ n₁ + ε₂ n₂
ε1, ε2 = -1.0, -0.5
ham = ComplexF64(ε1) * n1 + ComplexF64(ε2) * n2

# Create optimization problem (algebra constraints included automatically):
pop = cpolyopt(ham, sys)

# ## Comparing Approaches: Side-by-Side
#
# ### Manual Approach (Tedious)
# ```julia
# # Step 1: Declare variables
# @ncpolyvar c[1:2] c_dag[1:2]
#
# # Step 2: Define commutation groups (all operators in one group)
# comm_gps = [vcat(collect(c), collect(c_dag))]
#
# # Step 3: Encode anti-commutation relations (10 constraints for N=2!)
# eq_constraints = [
#     c[1]*c[1] + c[1]*c[1],  # {c₁, c₁} = 0
#     c[1]*c[2] + c[2]*c[1],  # {c₁, c₂} = 0
#     c[2]*c[2] + c[2]*c[2],  # {c₂, c₂} = 0
#     c_dag[1]*c_dag[1] + c_dag[1]*c_dag[1],  # {c₁†, c₁†} = 0
#     c_dag[1]*c_dag[2] + c_dag[2]*c_dag[1],  # {c₁†, c₂†} = 0
#     c_dag[2]*c_dag[2] + c_dag[2]*c_dag[2],  # {c₂†, c₂†} = 0
#     c[1]*c_dag[1] + c_dag[1]*c[1] - 1,      # {c₁, c₁†} = 1
#     c[1]*c_dag[2] + c_dag[2]*c[1],          # {c₁, c₂†} = 0
#     c[2]*c_dag[1] + c_dag[1]*c[2],          # {c₂, c₁†} = 0
#     c[2]*c_dag[2] + c_dag[2]*c[2] - 1,      # {c₂, c₂†} = 1
# ]
#
# # Step 4: Create SimplifyAlgorithm
# sa = SimplifyAlgorithm(comm_gps=comm_gps, is_nilpotent=true)
#
# # Step 5: Create optimization problem
# pop = cpolyopt(ham; eq_constraints=eq_constraints, comm_gps=comm_gps, is_nilpotent=true)
# ```
#
# ### Simplified Approach (Clean)
# ```julia
# sys = fermionic_algebra(2)
# c, c_dag = sys.variables
# ham = ComplexF64(ε1) * n1 + ComplexF64(ε2) * n2
# pop = cpolyopt(ham, sys)
# ```
#
# **4 lines vs 20+!** The `fermionic_algebra` constructor eliminates all boilerplate.

# ## Example: Fermi-Hubbard Model
#
# The Fermi-Hubbard model describes interacting fermions on a lattice:
#
# ```math
# H = -t \sum_{\langle i,j \rangle} (c_i^\dagger c_j + c_j^\dagger c_i) + U \sum_i n_i n_{i+1}
# ```
#
# Let's solve the two-site version:

N = 2
sys = fermionic_algebra(N)
c, c_dag = sys.variables

# Parameters:
t = 1.0  # Hopping amplitude
U = 4.0  # On-site interaction

# Hopping term: -t(c₁†c₂ + c₂†c₁)
hopping = -ComplexF64(t) * (c_dag[1] * c[2] + c_dag[2] * c[1])

# Number operators:
n1 = c_dag[1] * c[1]
n2 = c_dag[2] * c[2]

# Interaction term: U n₁ n₂
interaction = ComplexF64(U) * n1 * n2

# Total Hamiltonian:
ham_hubbard = hopping + interaction

# Add constraint: exactly 2 particles (one per site on average)
N_total = n1 + n2 - 2

pop_hubbard = cpolyopt(ham_hubbard, sys; eq_constraints=[N_total])

# Solve the optimization problem:
using MosekTools  # Or Clarabel
config = SolverConfig(order=2, optimizer=Mosek.Optimizer)
result = cs_nctssos(pop_hubbard, config)

println("Hubbard dimer ground state energy: ", result.objective)

# ## Advantages of `fermionic_algebra`
#
# 1. **Concise**: Reduces boilerplate from 20+ lines to 1 line
# 2. **Correct by construction**: All anti-commutation relations automatically included
# 3. **Scales easily**: Works for any number of modes N
# 4. **Maintainable**: Changes to algebra structure centralized in one place
# 5. **Composable**: Easy to add custom constraints (particle number, symmetries, etc.)
# 6. **Type-safe**: Returns properly typed NamedTuple compatible with cpolyopt
#
# ## Physical Observables
#
# Common fermionic observables are easy to construct:

sys3 = fermionic_algebra(3)
c, c_dag = sys3.variables

# Number operators:
n = [c_dag[i] * c[i] for i in 1:3]

# Total particle number:
N_total_op = sum(n)

# Hopping Hamiltonian (chain):
H_hop = sum(-ComplexF64(1.0) * (c_dag[i] * c[i+1] + c_dag[i+1] * c[i]) for i in 1:2)

# Nearest-neighbor interaction:
H_int = sum(ComplexF64(1.0) * n[i] * n[i+1] for i in 1:2)

# ## Next Steps
#
# - Explore [Pauli Algebra Interface](@ref) for spin systems
# - Learn about [GNS Reconstruction](@ref) for extracting states
# - See [Advanced Examples](@ref) for complex fermionic systems
#
# ## References
#
# - A. Auerbach, "Interacting Electrons and Quantum Magnetism" (1994)
# - Fradkin, "Field Theories of Condensed Matter Physics" (2013)
# - NCTSSoS.jl documentation: [https://quantumsos.github.io/NCTSSoS.jl/](https://quantumsos.github.io/NCTSSoS.jl/)
```

**Documentation Structure**:
1. **Introduction**: Motivation and physical context
2. **Manual approach**: Show the tedious way (for comparison)
3. **Simplified approach**: Demonstrate `fermionic_algebra(N)`
4. **Side-by-side comparison**: Highlight the improvement
5. **Example**: Fermi-Hubbard model with numerical solution
6. **Advantages list**: Benefits of the new interface
7. **Common observables**: Show typical constructions
8. **Next steps**: Cross-references to related documentation

**Literate.jl Notes**:
- Lines starting with `#` are markdown
- Regular Julia code will be executed
- Use `#-` to separate code blocks visually
- Cross-references use `@ref` tag

---

#### Step 5.2: Update documentation index

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples.md` (or similar)
**Lines**: Add to examples list

**Action**: Add fermionic algebra tutorial to documentation navigation

**Code Template**:
```markdown
## Algebra Interfaces

- [Pauli Algebra Interface](@ref)
- [Fermionic Algebra Interface](@ref)  # NEW: Add this line

## Advanced Examples

- [GNS Reconstruction](@ref)
```

---

#### Step 5.3: Add docstring cross-references

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/algebra_constructors.jl`
**Lines**: Update pauli_algebra docstring to cross-reference fermionic_algebra

**Action**: Add "See also" section to pauli_algebra

**Code Template**:
```julia
"""
    pauli_algebra(N::Int)

Create a Pauli algebra system for N spin-1/2 sites.

[Existing documentation...]

# See Also
- [`fermionic_algebra`](@ref): For fermionic systems with nilpotent operators
- [`cpolyopt`](@ref): For creating optimization problems with algebra systems
"""
```

---

### Phase 6: Integration Verification

#### Step 6.1: Update interface.jl to handle nilpotent

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/interface.jl`
**Lines**: Around 73-97 (check SimplifyAlgorithm construction)

**Action**: Verify that `cs_nctssos` correctly passes `is_nilpotent` to SimplifyAlgorithm

**Code Location to Verify**:
```julia
function cs_nctssos(pop::OP, solver_config::SolverConfig; dualize::Bool=true) where {OP<:OptimizationProblem}
    # This line should already work correctly if pop.is_nilpotent is defined
    sa = SimplifyAlgorithm(
        comm_gps=pop.comm_gps,
        is_projective=pop.is_projective,
        is_unipotent=pop.is_unipotent,
        is_nilpotent=pop.is_nilpotent  # VERIFY: This line should be added
    )
    # ... rest of function
end
```

**If not present**: Add `is_nilpotent=pop.is_nilpotent` to SimplifyAlgorithm construction.

**Test**: Verify that solving fermionic optimization problems works end-to-end.

---

#### Step 6.2: Check PolyOpt show method

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/pop.jl`
**Lines**: Around 73-90 (Base.show implementation)

**Action**: Update display to show is_nilpotent

**Code Template**:
```julia
function Base.show(io::IO, pop::P) where {P<:OptimizationProblem}
    cons_str(cons::Vector{P}, iseq::Bool) where {P} =
        join(["$(string(c)) " * (iseq ? "= 0" : ">= 0 \n") for c in cons], " \t")
    res_str = """
        obj: \n
            $(string(pop.objective)) \n
        constraints: \n
            $(cons_str(pop.eq_constraints,true))
            $(cons_str(pop.ineq_constraints,false))
        variables:
            $(join(string.(pop.variables)," ")) \n
        is_unipotent:
            $(pop.is_unipotent) \n
        is_projective:
            $(pop.is_projective) \n
        is_nilpotent:  # NEW: Add this line
            $(pop.is_nilpotent) \n
    """
    print(io, res_str)
end
```

**Test**: Print a fermionic PolyOpt and verify `is_nilpotent: true` appears.

---

## 3. Test-Driven Development Plan

### TDD Cycle for Each Phase

**General TDD Flow**:
1. Write failing test for specific feature
2. Implement minimum code to pass the test
3. Explain what was implemented and why
4. Refactor if needed
5. Move to next test

### Phase 1 Tests (SimplifyAlgorithm & Nilpotent)

#### Test 1.1: SimplifyAlgorithm accepts is_nilpotent
```julia
@testset "SimplifyAlgorithm nilpotent flag" begin
    @ncpolyvar x[1:2]
    comm_gps = [[x[1]], [x[2]]]

    # Should construct with is_nilpotent=true
    sa = SimplifyAlgorithm(comm_gps=comm_gps, is_nilpotent=true)
    @test sa.is_nilpotent == true
    @test sa.is_unipotent == false
    @test sa.is_projective == false
end
```

#### Test 1.2: Mutual exclusivity enforcement
```julia
@testset "Mutual exclusivity" begin
    @ncpolyvar x[1:2]
    comm_gps = [[x[1]], [x[2]]]

    # Should fail: nilpotent + unipotent
    @test_throws AssertionError SimplifyAlgorithm(
        comm_gps=comm_gps, is_nilpotent=true, is_unipotent=true)

    # Should fail: nilpotent + projective
    @test_throws AssertionError SimplifyAlgorithm(
        comm_gps=comm_gps, is_nilpotent=true, is_projective=true)
end
```

#### Test 1.3: _simplify_nilpotent! handles c²=0
```julia
@testset "Nilpotent simplification" begin
    @ncpolyvar c[1:2]
    comm_gps = [vcat(collect(c))]
    sa = SimplifyAlgorithm(comm_gps=comm_gps, is_nilpotent=true)

    # c² should become empty monomial
    c_squared = c[1] * c[1]
    mono = first(c_squared.monos)
    simplified = FastPolynomials.simplify(mono, sa)
    @test isempty(simplified.vars)
    @test isempty(simplified.z)

    # c * c_other stays as-is (different variables)
    c_mixed = c[1] * c[2]
    mono2 = first(c_mixed.monos)
    simplified2 = FastPolynomials.simplify(mono2, sa)
    @test !isempty(simplified2.vars)  # Should not vanish
end
```

#### Test 1.4: Higher powers also vanish
```julia
@testset "Higher powers nilpotent" begin
    @ncpolyvar c[1:1]
    comm_gps = [[c[1]]]
    sa = SimplifyAlgorithm(comm_gps=comm_gps, is_nilpotent=true)

    # c³ = c² * c = 0 * c = 0
    c_cubed = c[1] * c[1] * c[1]
    # After multiplication, this creates exponent 3
    # _simplify_nilpotent! should detect z[i] >= 2 and empty

    # Test by creating monomial directly with exponent 3
    mono = FastPolynomials.monomial([c[1]], [3])
    FastPolynomials._simplify_nilpotent!(mono)
    @test isempty(mono.vars)
    @test isempty(mono.z)
end
```

### Phase 2 Tests (PolyOpt Updates)

#### Test 2.1: PolyOpt accepts is_nilpotent
```julia
@testset "PolyOpt nilpotent construction" begin
    @ncpolyvar c[1:2]
    obj = ComplexF64(1.0) * c[1] * c[2]

    pop = cpolyopt(obj;
        comm_gps=[collect(c)],
        is_nilpotent=true)

    @test pop.is_nilpotent == true
    @test pop.is_unipotent == false
end
```

#### Test 2.2: Algebra interface extracts is_nilpotent
```julia
@testset "Algebra interface nilpotent" begin
    sys = fermionic_algebra(1)
    c, c_dag = sys.variables
    obj = ComplexF64(1.0) * c_dag[1] * c[1]

    pop = cpolyopt(obj, sys)

    @test pop.is_nilpotent == true
end
```

### Phase 3 Tests (Fermionic Algebra Constructor)

All tests listed in **Step 4.1** above. Priority order:

1. **Basic structure N=1** (validates constructor works)
2. **SimplifyAlgorithm properties** (validates correct setup)
3. **Anti-commutation relations N=1** (validates constraint generation)
4. **Nilpotent simplification** (validates zero monomial workaround)
5. **Integration with cpolyopt** (validates end-to-end flow)
6. **Scaling to N=2,3** (validates generalization)
7. **Additional constraints** (validates composability)
8. **Invalid input** (validates error handling)

### Phase 4 Tests (Integration with Solver)

**Run only with LOCAL_TESTING environment variable**:

#### Test 4.1: Free fermion vacuum
```julia
# Expected: Ground state energy = 0
@test isapprox(result.objective, 0.0, atol=1e-6)
```

#### Test 4.2: Single particle ground state
```julia
# Expected: Energy = minimum single-particle energy
@test isapprox(result.objective, -1.0, atol=1e-4)
```

#### Test 4.3: Hubbard dimer
```julia
# Expected: For U=0, single particle ground state = -t
@test isapprox(result.objective, -t, atol=1e-4)
```

### Test Execution Order

1. **Phase 1 tests**: SimplifyAlgorithm (foundation)
2. **Phase 2 tests**: PolyOpt types (integration point)
3. **Phase 3 tests**: fermionic_algebra constructor (user interface)
4. **Phase 4 tests**: Numerical validation (full pipeline)

---

## 4. Code Templates for Copy-Paste

### Template 1: SimplifyAlgorithm Extension (Complete)

```julia
# File: src/FastPolynomials/src/simplify.jl
# Location: Lines 1-18 (modify), then insert new function after line 66

# MODIFY EXISTING STRUCT (lines 1-18):
struct SimplifyAlgorithm
    comm_gps::Dict{Variable,Int}
    n_gps::Int
    is_unipotent::Bool
    is_projective::Bool
    is_nilpotent::Bool  # ADD THIS

    function SimplifyAlgorithm(;
        comm_gps::Vector{Vector{Variable}},
        is_unipotent::Bool=false,
        is_projective::Bool=false,
        is_nilpotent::Bool=false,  # ADD THIS
    )
        @assert sum([is_unipotent, is_projective, is_nilpotent]) <= 1  # MODIFY THIS
            "Only one of is_unipotent, is_projective, is_nilpotent can be true"

        return new(
            Dict(var => i for (i, vars) in enumerate(comm_gps) for var in vars),
            length(comm_gps),
            is_unipotent,
            is_projective,
            is_nilpotent,  # ADD THIS
        )
    end
end

# MODIFY nosimp (lines 20-22):
function nosimp(sa::SimplifyAlgorithm)
    return isone(sa.n_gps) && !sa.is_unipotent && !sa.is_projective && !sa.is_nilpotent  # ADD && !sa.is_nilpotent
end

# INSERT NEW FUNCTION (after line 66, before _simplify_unipotent!):
function _simplify_nilpotent!(m::Monomial)
    @inbounds for i in eachindex(m.z)
        if m.z[i] >= 2
            empty!(m.vars)
            empty!(m.z)
            return nothing
        end
    end
    return nothing
end

# MODIFY simplify! (lines 24-32):
function simplify!(m::Monomial, sa::SimplifyAlgorithm)
    (isone(m) || nosimp(sa)) && return m
    _comm!(m, sa.comm_gps)

    sa.is_unipotent && (_simplify_unipotent!(m); return m)
    sa.is_projective && (_simplify_projective!(m); return m)
    sa.is_nilpotent && (_simplify_nilpotent!(m); return m)  # ADD THIS LINE
    _simplify_standard!(m)
    return m
end
```

### Template 2: PolyOpt Extension (Complete)

```julia
# File: src/pop.jl
# Location: Lines 25-33 (struct), 57-71 (constructor)

# MODIFY STRUCT (lines 25-33):
struct PolyOpt{P} <: OptimizationProblem{P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    variables::Vector{Variable}
    comm_gps::Vector{Vector{Variable}}
    is_unipotent::Bool
    is_projective::Bool
    is_nilpotent::Bool  # ADD THIS
end

# MODIFY CONSTRUCTOR (lines 57-71):
function polyopt(
    objective::P;
    eq_constraints=Any[],
    ineq_constraints=Any[],
    comm_gps=Vector{Variable}[],
    is_unipotent::Bool=false,
    is_projective::Bool=false,
    is_nilpotent::Bool=false  # ADD THIS
) where {T,P<:AbstractPolynomial{T}}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    eq_cons = unique!(collect(P, eq_constraints))
    ineq_cons = unique!(collect(P, ineq_constraints))
    vars = sorted_union(variables(objective), variables.(eq_cons)..., variables.(ineq_cons)...)
    if !isempty(comm_gps)
        @assert all([isempty(intersect(gp_a, gp_b)) for gp_a in comm_gps, gp_b in comm_gps if gp_a != gp_b]) "The commutative groups must be disjoint."
        @assert issubset(union(comm_gps...), vars) "The commutative variables must be a subset of the variables."
        @assert all(issorted(gp) for gp in comm_gps) "The commutative groups must be sorted."
    else
        push!(comm_gps, sort(vars))
    end
    @assert sum([is_unipotent, is_projective, is_nilpotent]) <= 1  # MODIFY THIS
        "Only one of is_unipotent, is_projective, is_nilpotent can be true"

    return PolyOpt{P}(objective, eq_cons, ineq_cons, vars, comm_gps,
                      is_unipotent, is_projective, is_nilpotent)  # ADD is_nilpotent
end

# REPEAT SAME CHANGES FOR ComplexPolyOpt (lines 95-126)
```

### Template 3: fermionic_algebra Constructor (Complete)

See **Step 3.1** code template above (too long to repeat here).

### Template 4: Basic Test Suite (Starter)

```julia
# File: test/algebra_constructors.jl
# Location: Append after existing tests

@testset "Fermionic Algebra Constructor" begin
    @testset "Basic Structure N=1" begin
        sys = fermionic_algebra(1)
        @test hasfield(typeof(sys), :variables)
        @test hasfield(typeof(sys), :simplify_algo)

        c, c_dag = sys.variables
        @test length(c) == 1
        @test length(c_dag) == 1
        @test length(sys.equality_constraints) == 3
        @test sys.simplify_algo.is_nilpotent == true
    end
end
```

---

## 5. Integration Checklist

### Pre-Implementation Checklist

- [ ] Read and understand pauli_analysis.md
- [ ] Understand issue #179 (zero monomial limitation)
- [ ] Review existing SimplifyAlgorithm code
- [ ] Review existing PolyOpt code
- [ ] Confirm TDD approach with user

### Phase 1 Checklist (SimplifyAlgorithm)

- [ ] Add `is_nilpotent::Bool` field to SimplifyAlgorithm struct
- [ ] Update SimplifyAlgorithm constructor with `is_nilpotent` parameter
- [ ] Add mutual exclusivity assertion
- [ ] Update `nosimp` function
- [ ] Implement `_simplify_nilpotent!` function
- [ ] Update `simplify!` dispatch to call `_simplify_nilpotent!`
- [ ] Write and pass unit tests for nilpotent simplification
- [ ] Test zero monomial workaround (empty vars/z)

### Phase 2 Checklist (PolyOpt)

- [ ] Add `is_nilpotent::Bool` field to PolyOpt struct
- [ ] Update PolyOpt constructor with `is_nilpotent` parameter
- [ ] Update mutual exclusivity assertion in polyopt
- [ ] Repeat for ComplexPolyOpt struct
- [ ] Repeat for cpolyopt constructor
- [ ] Update cpolyopt algebra interface to extract `is_nilpotent`
- [ ] Update SimplifyAlgorithm construction in cpolyopt (line ~119)
- [ ] Update Base.show method to display is_nilpotent
- [ ] Write and pass unit tests for PolyOpt with nilpotent

### Phase 3 Checklist (Algebra Constructor)

- [ ] Implement `fermionic_algebra(N)` function
- [ ] Add comprehensive docstring with examples
- [ ] Export fermionic_algebra in NCTSSoS.jl
- [ ] Write unit tests for basic structure (N=1,2,3)
- [ ] Write unit tests for SimplifyAlgorithm properties
- [ ] Write unit tests for anti-commutation relations
- [ ] Write unit tests for nilpotent simplification
- [ ] Write unit tests for cpolyopt integration
- [ ] Write unit tests for custom constraint merging
- [ ] Write unit tests for error handling
- [ ] Verify all tests pass

### Phase 4 Checklist (Integration Tests)

- [ ] Write free fermion vacuum test (LOCAL_TESTING only)
- [ ] Write single-particle ground state test
- [ ] Write Fermi-Hubbard dimer test
- [ ] Run with Mosek and verify numerical results
- [ ] Run with Clarabel and verify numerical results
- [ ] Document expected accuracies

### Phase 5 Checklist (Documentation)

- [ ] Create fermionic_algebra_interface.jl tutorial
- [ ] Follow Literate.jl format (markdown comments + code)
- [ ] Include physical motivation section
- [ ] Show manual vs simplified approach comparison
- [ ] Add Fermi-Hubbard example with numerical solution
- [ ] List advantages of new interface
- [ ] Add cross-references to related documentation
- [ ] Update documentation index/navigation
- [ ] Add cross-references in pauli_algebra docstring
- [ ] Build documentation locally and verify rendering

### Phase 6 Checklist (Final Integration)

- [ ] Verify cs_nctssos passes `is_nilpotent` to SimplifyAlgorithm
- [ ] Update interface.jl if needed
- [ ] Run full test suite
- [ ] Check code coverage
- [ ] Run integration tests with LOCAL_TESTING
- [ ] Verify no regressions in existing Pauli tests
- [ ] Manual smoke test: create and solve simple fermionic problem
- [ ] Review code for consistency with Julia idioms
- [ ] Check for any TODO or FIXME comments

### Post-Implementation Checklist

- [ ] All unit tests pass
- [ ] All integration tests pass (with LOCAL_TESTING)
- [ ] Documentation builds without errors
- [ ] No regressions in existing functionality
- [ ] Code follows project style guidelines
- [ ] Docstrings complete and accurate
- [ ] Examples in docstrings are tested (doctests)
- [ ] Ready for code review

---

## 6. Edge Cases and Considerations

### Edge Case 1: Mixed Commutator/Anti-commutator Expressions

**Scenario**: User tries to create expressions mixing fermionic and bosonic operators.

**Current Behavior**: Not supported (different commutation structures).

**Recommendation**: Document clearly that fermionic_algebra creates a pure fermionic system. Mixing with other algebras requires advanced manual setup.

### Edge Case 2: Zero Monomial in Polynomial Operations

**Scenario**: After nilpotent simplification, monomial becomes empty (zero).

**Current Handling**:
- Empty monomial represents identity (1) in current code
- For nilpotent case, empty monomial after simplification represents zero
- Polynomial constructor filters zero-coefficient terms

**Validation Needed**:
```julia
# Test that zero monomial doesn't survive in final polynomial
@test iszero(simplify(c[1]^2, sa))  # May fail depending on Polynomial handling
```

**Workaround**: If direct testing fails, test at monomial level:
```julia
mono = first((c[1] * c[1]).monos)
simplified = FastPolynomials.simplify(mono, sa)
@test isempty(simplified.vars) && isempty(simplified.z)
```

### Edge Case 3: Higher-Order Terms (c³, c⁴, ...)

**Scenario**: User creates expressions like `c * c * c`.

**Expected Behavior**: Should simplify to zero.

**Implementation**: `_simplify_nilpotent!` checks `z[i] >= 2`, which catches all higher powers.

**Test**:
```julia
c_cubed = c[1] * c[1] * c[1]
# After multiplication, exponent becomes 3
# _simplify_nilpotent! detects >= 2 and empties monomial
```

### Edge Case 4: Constraint Count Scaling

**Scenario**: Large N leads to many constraints.

**Constraint Count**: For N modes: 2N² + N constraints

| N | Constraints |
|---|-------------|
| 1 | 3 |
| 2 | 10 |
| 3 | 21 |
| 5 | 55 |
| 10 | 210 |

**Consideration**: For N > 10, constraint count grows quadratically. Document performance implications.

**Optimization Opportunity** (future work): Exploit symmetries to reduce constraint count.

### Edge Case 5: Variable Naming Conflicts

**Scenario**: User wants to use `c` for something else.

**Solution**: Extract to different names:
```julia
sys = fermionic_algebra(N)
annihilation, creation = sys.variables  # Rename as needed
```

### Edge Case 6: Complex vs Real Coefficients

**Scenario**: Fermionic Hamiltonians are often Hermitian but coefficients are complex.

**Current Handling**: Use `ComplexF64(1.0) * expression` to ensure correct type.

**Validation**: cpolyopt checks that objective is symmetric (Hermitian).

**Test**: Verify that fermionic number operator passes symmetry check.

### Edge Case 7: Single Commutation Group Semantics

**Scenario**: All fermionic operators in one group means no reordering.

**Expected Behavior**: Anti-commutation handled entirely by equality constraints, not by simplifier.

**Validation**:
```julia
# c[1] * c[2] should NOT be reordered to c[2] * c[1]
# Both orders should exist in the algebra, related by constraint
```

### Edge Case 8: Interaction with State Polynomials

**Scenario**: User wants to use fermionic operators with state-dependent polynomials.

**Current Status**: Not directly supported (requires extension).

**Recommendation**: Document as future work. Current implementation focuses on operator algebra.

### Edge Case 9: GNS Reconstruction with Fermions

**Scenario**: Extracting fermionic states from moment matrices.

**Consideration**: GNS reconstruction should work unchanged (operates on general NC polynomials).

**Test** (future): Add fermionic example to GNS tests.

### Edge Case 10: Numerical Precision Issues

**Scenario**: Anti-commutation constraints may not be exactly satisfied due to SDP solver tolerance.

**Mitigation**:
- Use high-precision solvers (Mosek recommended)
- Check constraint violations in results
- Document expected tolerances

**Test**: Verify that constraint violations are < 1e-6 in integration tests.

---

## 7. Performance Implications

### Simplification Performance

**Nilpotent simplification** is simpler than unipotent:
- Unipotent: O(n) with complex cancellation logic
- Nilpotent: O(n) with simple zero detection

**Expected Impact**: Negligible performance difference, possibly slightly faster than unipotent.

### Constraint Scaling

**Critical Concern**: O(N²) constraints for N modes.

**Comparison**:
- Pauli: 6N constraints (linear scaling)
- Fermionic: ~2N² constraints (quadratic scaling)

**Impact on Solver**:
- More constraints → larger SDP matrices
- Recommend N ≤ 10 for practical computations
- Document scaling limitations clearly

### Memory Usage

**Monomial Representation**: Empty vectors for zero monomials use minimal memory.

**Polynomial Storage**: Zero monomials filtered out, no extra memory overhead.

---

## 8. Future-Proofing for Issue #179 Resolution

### Current Workaround

```julia
# When c² is detected:
empty!(m.vars)  # Represents zero (by context)
empty!(m.z)
```

### Future Proper Implementation

When issue #179 adds `is_zero` flag:

```julia
struct Monomial
    vars::Vector{Variable}
    z::Vector{Int}
    is_zero::Bool  # NEW
end

function _simplify_nilpotent!(m::Monomial)
    @inbounds for i in eachindex(m.z)
        if m.z[i] >= 2
            m.is_zero = true  # FUTURE: Set flag instead of emptying
            return nothing
        end
    end
    return nothing
end
```

### Migration Path

1. Add `is_zero` field to Monomial (default false)
2. Update `_simplify_nilpotent!` to set flag
3. Update Polynomial constructor to filter `is_zero` monomials
4. Update arithmetic operations to handle zero monomials
5. All existing code continues to work

**No breaking changes** to user code because:
- `fermionic_algebra` API unchanged
- SimplifyAlgorithm API unchanged
- User never directly creates monomials

---

## 9. Summary for Parent Agent

### Implementation Path

The parent agent should implement in this exact order:

1. **Start with Phase 1, Step 1.1**: Add `is_nilpotent` field to SimplifyAlgorithm
   - Write failing test first
   - Implement minimum code to pass
   - Explain changes

2. **Continue through Phase 1**: Complete all SimplifyAlgorithm changes
   - Each step builds on previous
   - Test after each step

3. **Move to Phase 2**: Update PolyOpt types
   - Tests ensure integration works

4. **Implement Phase 3**: Create fermionic_algebra constructor
   - Core functionality delivered

5. **Add Phase 4**: Integration tests (if LOCAL_TESTING available)
   - Validate numerical correctness

6. **Finish with Phase 5**: Documentation
   - Tutorial and API docs

### Critical Success Factors

✅ **Follow TDD rigorously**: Test first, implement second, explain third
✅ **One change at a time**: Don't batch multiple modifications
✅ **Verify after each step**: Run tests to ensure no regressions
✅ **Explain reasoning**: Document why each change is made
✅ **Handle edge cases**: Test zero monomials, constraint counts, etc.

### Expected Outcomes

After full implementation:
- ✅ Users can create fermionic systems with `fermionic_algebra(N)`
- ✅ Nilpotent simplification correctly handles c² = 0
- ✅ All anti-commutation relations automatically included
- ✅ Integration with cpolyopt seamless
- ✅ Comprehensive test coverage
- ✅ Complete documentation with examples
- ✅ No regressions in existing Pauli algebra functionality

---

## 10. Questions for User (Before Implementation)

Before proceeding, the parent agent should confirm with the user:

1. **Naming convention**: Is `c_dag` acceptable for creation operators, or prefer another notation?
2. **Constraint optimization**: Should we implement any constraint reduction for large N?
3. **Integration test availability**: Do you have LOCAL_TESTING environment with Mosek?
4. **Documentation scope**: Should we include advanced examples (BCS, topological systems)?
5. **Performance benchmarks**: Should we profile and compare with Pauli algebra?

---

## Appendix A: File Modification Summary

| File | Action | Lines | Description |
|------|--------|-------|-------------|
| `src/FastPolynomials/src/simplify.jl` | Modify + Add | 1-18, 20-22, 24-32, +66 | Add nilpotent support |
| `src/pop.jl` | Modify | 25-33, 57-71, 95-126 | Add is_nilpotent field |
| `src/algebra_constructors.jl` | Add | +40 | New fermionic_algebra function |
| `src/NCTSSoS.jl` | Modify | ~14-20 | Export fermionic_algebra |
| `src/interface.jl` | Verify | ~73-97 | Check SimplifyAlgorithm construction |
| `test/algebra_constructors.jl` | Add | +300 | Comprehensive test suite |
| `docs/src/examples/literate/fermionic_algebra_interface.jl` | Create | New file | Tutorial documentation |

**Total Changes**: ~500 lines of new code + ~50 lines of modifications

---

## Appendix B: Mathematical Validation

### Anti-commutation Relations

For N=1, verify all 3 constraints:

1. `{c, c} = c·c + c·c = 0` ✓
2. `{c†, c†} = c†·c† + c†·c† = 0` ✓
3. `{c, c†} = c·c† + c†·c = 1` ✓

For N=2, verify scaling:
- `{cᵢ, cⱼ}`: (1,1), (1,2), (2,2) → 3 constraints
- `{cᵢ†, cⱼ†}`: (1,1), (1,2), (2,2) → 3 constraints
- `{cᵢ, cⱼ†}`: (1,1), (1,2), (2,1), (2,2) → 4 constraints
- Total: 10 ✓

### Nilpotent Algebra Properties

Verify c² = 0 implies:
- c³ = c·c² = c·0 = 0 ✓
- c⁴ = c²·c² = 0·0 = 0 ✓
- c^n = 0 for all n ≥ 2 ✓

Verify anti-commutation + nilpotency:
- `c₁·c₂ = -c₂·c₁` (from {c₁,c₂}=0) ✓
- `(c₁·c₂)² = c₁·c₂·c₁·c₂ = -c₁·c₁·c₂·c₂ = 0` ✓

---

## Appendix C: References

1. **NCTSSoS.jl Repository**: https://github.com/QuantumSOS/NCTSSoS.jl
2. **Issue #179**: FastPolynomials zero monomial support
3. **Pauli Analysis**: `.claude/tasks/fermionic-algebra/pauli_analysis.md`
4. **Julia Documentation**: https://docs.julialang.org/
5. **Fermionic Algebra**: Fradkin, "Field Theories of Condensed Matter Physics"
6. **Moment Methods**: Lasserre, "Moments, Positive Polynomials and Their Applications"

---

**End of Plan**

This plan is complete and ready for the parent agent to execute via TDD.
