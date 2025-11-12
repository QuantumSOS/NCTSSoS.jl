# Fermionic Algebra Implementation Plan (Revised - Simplified)

**Date**: 2025-11-12
**Author**: julia-development sub-agent
**Status**: REVISED - Constraint-Based Approach (No Core Modifications)
**Previous Version**: 2025-11-05 (archived in git history)

---

## Executive Summary

This revised plan implements fermionic algebra in NCTSSoS.jl using a **simplified, constraint-based approach** that requires **NO modifications** to FastPolynomials, PolyOpt, or ComplexPolyOpt. The implementation follows the Pauli algebra pattern exactly:

1. **Define AbstractAlgebra type hierarchy** with struct definitions
2. **Create FermionicAlgebra struct** to encapsulate algebra data
3. **Create fermionic_algebra(N) constructor** returning FermionicAlgebra instance
4. **Encode all fermionic properties as polynomial constraints**:
   - Anti-commutation relations: {c_i, c_j} = 0, {c†_i, c†_j} = 0
   - Canonical anti-commutation: {c_i, c†_j} = δ_ij
   - Nilpotent constraints: c_i² = 0, (c†_i)² = 0
5. **Use existing SimplifyAlgorithm** without modifications
6. **Follow pauli_algebra interface** for consistency

**Key Difference from Previous Plan**: Nilpotency is enforced through polynomial equality constraints (c² = 0), NOT through monomial simplification rules. This eliminates the need for zero monomial support or any core type modifications.

---

## 1. Architecture Overview

### 1.1 Core Design Philosophy

**Constraint-Based Approach**:
- ALL fermionic algebra properties are encoded as polynomial equality constraints
- NO special simplification algorithms needed
- NO modifications to existing types
- SimplifyAlgorithm handles commutation groups ONLY (not nilpotency)

**Pattern to Follow**: `pauli_algebra` implementation
- Return concrete algebra type (FermionicAlgebra) with standard fields
- Use SimplifyAlgorithm with existing flags (no new flags)
- Encode algebra structure entirely through constraints
- Integration with cpolyopt(objective, algebra) interface

### 1.2 Fermionic Algebra Structure

**Variables**:
- `c[i]`: Annihilation operators (destroy a fermion in mode i)
- `c_dag[i]`: Creation operators (create a fermion in mode i)

**Commutation Groups**:
- Single group containing ALL operators: `[c[1], c_dag[1], c[2], c_dag[2], ...]`
- No automatic reordering (operators in same group don't commute/anti-commute automatically)
- All algebraic relations encoded explicitly as constraints

**Constraint Types**:
1. **Anti-commutation (c, c)**: c_i·c_j + c_j·c_i = 0 for all i,j
2. **Anti-commutation (c†, c†)**: c†_i·c†_j + c†_j·c†_i = 0 for all i,j
3. **Canonical anti-commutation**: c_i·c†_j + c†_j·c_i = δ_ij
4. **Nilpotent constraints**: c_i² = 0 and (c†_i)² = 0 for all i

### 1.3 Why This Approach Works

**Nilpotency as Constraint**:
- Instead of `_simplify_nilpotent!` detecting c² and making it zero
- We add explicit constraint: c² = 0
- SDP solver enforces this during optimization
- No need for zero monomial representation

**Comparison with Previous Plan**:

| Aspect | Previous Plan | Revised Plan |
|--------|---------------|---------------|
| Nilpotency | Simplification rule | Equality constraint |
| FastPolynomials | Modify (_simplify_nilpotent!) | No changes |
| PolyOpt types | Add is_nilpotent field | No changes |
| SimplifyAlgorithm | Add is_nilpotent flag | No changes |
| Zero monomials | Workaround needed | Not needed |
| Implementation complexity | High (100+ line changes) | Low (40 lines total) |

**Trade-offs**:
- ✅ Much simpler implementation
- ✅ No core type modifications
- ✅ No workarounds for zero monomials
- ✅ Easier to maintain
- ⚠️ More constraints for solver (includes c_i² = 0 explicitly)
- ⚠️ No automatic simplification of c² terms in expressions

---

## 2. Detailed Implementation Plan

### Phase 0: Define AbstractAlgebra Type Hierarchy and Structs

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/algebra_constructors.jl` (add at top)

**Action**: Define abstract type and concrete struct for algebras

**Structure**:
```julia
"""
    AbstractAlgebra

Abstract base type for all algebra systems in NCTSSoS.
Concrete subtypes include PauliAlgebra and FermionicAlgebra.
"""
abstract type AbstractAlgebra end

"""
    PauliAlgebra <: AbstractAlgebra

Pauli algebra system for spin-1/2 sites.

# Fields
- `N::Int`: Number of spin-1/2 sites
- `variables::Tuple`: Tuple of (x, y, z) variable arrays
- `simplify_algo::SimplifyAlgorithm`: Simplification algorithm configuration
- `equality_constraints::Vector`: Polynomial equality constraints
- `inequality_constraints::Vector`: Polynomial inequality constraints
- `comm_gps::Vector`: Commutation groups
"""
struct PauliAlgebra <: AbstractAlgebra
    N::Int
    variables::Tuple
    simplify_algo::SimplifyAlgorithm
    equality_constraints::Vector
    inequality_constraints::Vector
    comm_gps::Vector
end

"""
    FermionicAlgebra <: AbstractAlgebra

Fermionic algebra system for fermionic modes.

# Fields
- `N::Int`: Number of fermionic modes
- `variables::Tuple`: Tuple of (c, c_dag) variable arrays
- `simplify_algo::SimplifyAlgorithm`: Simplification algorithm configuration
- `equality_constraints::Vector`: Polynomial equality constraints
- `inequality_constraints::Vector`: Polynomial inequality constraints
- `comm_gps::Vector`: Commutation groups
"""
struct FermionicAlgebra <: AbstractAlgebra
    N::Int
    variables::Tuple
    simplify_algo::SimplifyAlgorithm
    equality_constraints::Vector
    inequality_constraints::Vector
    comm_gps::Vector
end
```

**Implementation Notes**:
1. Add this code at the top of `algebra_constructors.jl` before function definitions
2. Update `pauli_algebra(N)` to return `PauliAlgebra(N, ...)` instead of NamedTuple
3. `fermionic_algebra(N)` will return `FermionicAlgebra(N, ...)`
4. This provides type safety and better IDE support
5. Export both concrete types from `NCTSSoS.jl`

**Testing Approach**:
- Test that types are defined and can be constructed
- Test that field access works correctly
- Test type hierarchy with `<:` operator

---

### Phase 1: Implement fermionic_algebra Constructor

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/algebra_constructors.jl`

**Action**: Create `fermionic_algebra(N)` function following `pauli_algebra` pattern

**Structure**:
```julia
function fermionic_algebra(N::Int)
    # 1. Input validation
    # 2. Declare variables (@ncpolyvar)
    # 3. Define commutation groups
    # 4. Build equality constraints (anti-commutation + nilpotent)
    # 5. Create SimplifyAlgorithm
    # 6. Return FermionicAlgebra instance
end
```

**Detailed Steps**:

#### Step 1.1: Input Validation
```julia
@assert N >= 1 "Number of modes N must be at least 1"
```

#### Step 1.2: Declare Variables
```julia
@ncpolyvar c[1:N] c_dag[1:N]
```

**Naming Convention**:
- `c[i]`: Annihilation operator for mode i
- `c_dag[i]`: Creation operator for mode i (using _dag suffix, not superscript †)

#### Step 1.3: Define Commutation Groups
```julia
# All fermionic operators in a single group (no automatic commutation)
all_ops = vcat(collect(c), collect(c_dag))
comm_gps = [all_ops]
```

**Rationale**:
- Fermionic operators neither commute nor anti-commute automatically
- Putting all in one group prevents automatic reordering
- All relations handled explicitly via constraints

#### Step 1.4: Build Equality Constraints

**Constraint Categories**:

1. **Anti-commutation {c_i, c_j} = 0**:
```julia
# For all i,j with i ≤ j (to avoid duplicates)
for i in 1:N, j in i:N
    push!(equality_constraints, c[i] * c[j] + c[j] * c[i])
end
```
Count: N(N+1)/2 constraints

2. **Anti-commutation {c†_i, c†_j} = 0**:
```julia
for i in 1:N, j in i:N
    push!(equality_constraints, c_dag[i] * c_dag[j] + c_dag[j] * c_dag[i])
end
```
Count: N(N+1)/2 constraints

3. **Canonical anti-commutation {c_i, c†_j} = δ_ij**:
```julia
for i in 1:N, j in 1:N
    if i == j
        # {c_i, c†_i} = 1
        push!(equality_constraints, c[i] * c_dag[i] + c_dag[i] * c[i] - 1)
    else
        # {c_i, c†_j} = 0 for i ≠ j
        push!(equality_constraints, c[i] * c_dag[j] + c_dag[j] * c[i])
    end
end
```
Count: N² constraints

**Total so far**: N(N+1)/2 + N(N+1)/2 + N² = N² + N(N+1) = 2N² + N

4. **Nilpotent constraints** (NEW - Key simplification):
```julia
# c_i² = 0 for all i
for i in 1:N
    push!(equality_constraints, c[i] * c[i])
end

# (c†_i)² = 0 for all i
for i in 1:N
    push!(equality_constraints, c_dag[i] * c_dag[i])
end
```
Count: 2N additional constraints

**IMPORTANT NOTE**: These nilpotent constraints (c² = 0) are REDUNDANT with the anti-commutation constraints ({c,c} = 2c² = 0 ⟹ c² = 0). However, we include them for two reasons:
1. **Clarity**: Makes nilpotency explicit in the constraint set
2. **Solver efficiency**: Some solvers may benefit from having this constraint explicitly

**Final constraint count**: 2N² + N + 2N = 2N² + 3N

**Constraint Count Verification**:
- N=1: 2(1) + 3(1) = 5 constraints
- N=2: 2(4) + 3(2) = 14 constraints
- N=3: 2(9) + 3(3) = 27 constraints

#### Step 1.5: Create SimplifyAlgorithm
```julia
simplify_algo = SimplifyAlgorithm(
    comm_gps=comm_gps,
    is_unipotent=false,
    is_projective=false
)
```

**Note**: We use the EXISTING SimplifyAlgorithm without any modifications. No `is_nilpotent` flag needed.

#### Step 1.6: Create Inequality Constraints
```julia
inequality_constraints = empty(equality_constraints)
```

No inequality constraints for basic fermionic algebra.

#### Step 1.7: Return FermionicAlgebra Instance
```julia
return FermionicAlgebra(
    N,
    (c, c_dag),
    simplify_algo,
    equality_constraints,
    inequality_constraints,
    comm_gps
)
```

**Complete Implementation Template**:
```julia
"""
    fermionic_algebra(N::Int)

Create a fermionic algebra system for N fermionic modes.

Fermionic operators satisfy anti-commutation relations and nilpotency:
- Anti-commutation: {cᵢ, cⱼ} = 0, {cᵢ†, cⱼ†} = 0
- Canonical anti-commutation: {cᵢ, cⱼ†} = δᵢⱼ
- Nilpotency: cᵢ² = 0, (cᵢ†)² = 0

All properties are enforced through polynomial equality constraints.

# Arguments
- `N::Int`: Number of fermionic modes

# Returns
A `NamedTuple` with fields:
- `variables::Tuple`: (c, c_dag) where c[i] and c_dag[i] are annihilation/creation operators
- `simplify_algo::SimplifyAlgorithm`: Simplification algorithm for commutation groups
- `equality_constraints::Vector`: All fermionic algebra constraints
- `inequality_constraints::Vector`: Empty (no inequality constraints)
- `comm_gps::Vector{Vector{Variable}}`: Commutation group structure

# Example
```julia
using NCTSSoS

# Create 2-mode fermionic system
sys = fermionic_algebra(2)
c, c_dag = sys.variables

# Number operator for mode 1
n1 = c_dag[1] * c[1]

# Create optimization problem
pop = cpolyopt(n1, sys)
```

# Mathematical Details
For N modes, the system includes:
- N(N+1)/2 constraints: {cᵢ, cⱼ} = 0
- N(N+1)/2 constraints: {cᵢ†, cⱼ†} = 0
- N² constraints: {cᵢ, cⱼ†} = δᵢⱼ
- 2N constraints: cᵢ² = 0, (cᵢ†)² = 0

Total: 2N² + 3N equality constraints
"""
function fermionic_algebra(N::Int)
    @assert N >= 1 "Number of modes N must be at least 1"

    # Step 1: Declare variables
    @ncpolyvar c[1:N] c_dag[1:N]

    # Step 2: Define commutation groups
    # All operators in single group (no automatic commutation)
    all_ops = vcat(collect(c), collect(c_dag))
    comm_gps = [all_ops]

    # Step 3: Build equality constraints
    equality_constraints = Polynomial{ComplexF64}[]

    # Anti-commutation: {cᵢ, cⱼ} = 0
    for i in 1:N, j in i:N
        push!(equality_constraints, c[i] * c[j] + c[j] * c[i])
    end

    # Anti-commutation: {cᵢ†, cⱼ†} = 0
    for i in 1:N, j in i:N
        push!(equality_constraints, c_dag[i] * c_dag[j] + c_dag[j] * c_dag[i])
    end

    # Canonical anti-commutation: {cᵢ, cⱼ†} = δᵢⱼ
    for i in 1:N, j in 1:N
        if i == j
            push!(equality_constraints, c[i] * c_dag[i] + c_dag[i] * c[i] - 1)
        else
            push!(equality_constraints, c[i] * c_dag[j] + c_dag[j] * c[i])
        end
    end

    # Nilpotent constraints: cᵢ² = 0, (cᵢ†)² = 0
    for i in 1:N
        push!(equality_constraints, c[i] * c[i])
        push!(equality_constraints, c_dag[i] * c_dag[i])
    end

    # Step 4: Create SimplifyAlgorithm (no special flags needed)
    simplify_algo = SimplifyAlgorithm(
        comm_gps=comm_gps,
        is_unipotent=false,
        is_projective=false
    )

    # Step 5: No inequality constraints
    inequality_constraints = empty(equality_constraints)

    # Step 6: Return algebra system
    return (
        variables=(c, c_dag),
        simplify_algo=simplify_algo,
        equality_constraints=equality_constraints,
        inequality_constraints=inequality_constraints,
        comm_gps=comm_gps
    )
end
```

---

### Phase 2: Export fermionic_algebra

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/NCTSSoS.jl`

**Action**: Add to module exports

**Location**: Around line 14-20

**Code**:
```julia
export pauli_algebra
export fermionic_algebra  # ADD THIS LINE
```

---

### Phase 3: Comprehensive Testing

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/test/algebra_constructors.jl`

**Action**: Add test suite for fermionic algebra

#### Test 3.1: Basic Structure Validation (N=1)

```julia
@testset "Fermionic Algebra Constructor" begin
    @testset "Basic Structure N=1" begin
        sys = fermionic_algebra(1)

        # Verify NamedTuple structure (matches Pauli pattern)
        @test hasfield(typeof(sys), :variables)
        @test hasfield(typeof(sys), :simplify_algo)
        @test hasfield(typeof(sys), :equality_constraints)
        @test hasfield(typeof(sys), :inequality_constraints)
        @test hasfield(typeof(sys), :comm_gps)

        # Verify variables
        c, c_dag = sys.variables
        @test length(c) == 1
        @test length(c_dag) == 1

        # Verify constraint count for N=1: 2(1)² + 3(1) = 5
        @test length(sys.equality_constraints) == 5
        @test isempty(sys.inequality_constraints)

        # Verify commutation groups
        @test length(sys.comm_gps) == 1
        @test length(sys.comm_gps[1]) == 2  # c[1] and c_dag[1]
    end
end
```

#### Test 3.2: Constraint Count Scaling

```julia
@testset "Constraint Scaling" begin
    # N=1: 2(1)² + 3(1) = 5
    sys1 = fermionic_algebra(1)
    @test length(sys1.equality_constraints) == 5

    # N=2: 2(4) + 3(2) = 14
    sys2 = fermionic_algebra(2)
    @test length(sys2.equality_constraints) == 14

    # N=3: 2(9) + 3(3) = 27
    sys3 = fermionic_algebra(3)
    @test length(sys3.equality_constraints) == 27

    # General formula verification
    for N in 1:5
        sys = fermionic_algebra(N)
        expected_count = 2*N^2 + 3*N
        @test length(sys.equality_constraints) == expected_count
    end
end
```

#### Test 3.3: Specific Constraints Verification (N=1)

```julia
@testset "Constraint Verification N=1" begin
    sys = fermionic_algebra(1)
    c, c_dag = sys.variables
    constraints = sys.equality_constraints

    # Should contain all 5 constraints:
    # 1. {c, c} = 2c² = 0
    @test c[1] * c[1] + c[1] * c[1] in constraints

    # 2. {c†, c†} = 2(c†)² = 0
    @test c_dag[1] * c_dag[1] + c_dag[1] * c_dag[1] in constraints

    # 3. {c, c†} = 1
    @test c[1] * c_dag[1] + c_dag[1] * c[1] - 1 in constraints

    # 4. c² = 0 (nilpotent constraint)
    @test c[1] * c[1] in constraints

    # 5. (c†)² = 0 (nilpotent constraint)
    @test c_dag[1] * c_dag[1] in constraints
end
```

#### Test 3.4: Specific Constraints Verification (N=2)

```julia
@testset "Constraint Verification N=2" begin
    sys = fermionic_algebra(2)
    c, c_dag = sys.variables
    constraints = sys.equality_constraints

    # Spot check key constraints

    # Anti-commutation {c₁, c₂} = 0
    @test c[1] * c[2] + c[2] * c[1] in constraints

    # Anti-commutation {c₁†, c₂†} = 0
    @test c_dag[1] * c_dag[2] + c_dag[2] * c_dag[1] in constraints

    # Canonical {c₁, c₁†} = 1
    @test c[1] * c_dag[1] + c_dag[1] * c[1] - 1 in constraints

    # Canonical {c₁, c₂†} = 0
    @test c[1] * c_dag[2] + c_dag[2] * c[1] in constraints

    # Nilpotent c₁² = 0
    @test c[1] * c[1] in constraints

    # Nilpotent (c₂†)² = 0
    @test c_dag[2] * c_dag[2] in constraints
end
```

#### Test 3.5: SimplifyAlgorithm Properties

```julia
@testset "SimplifyAlgorithm Properties" begin
    sys = fermionic_algebra(2)
    sa = sys.simplify_algo

    @test sa.is_unipotent == false
    @test sa.is_projective == false
    @test sa.n_gps == 1  # Single commutation group
end
```

#### Test 3.6: Integration with cpolyopt

```julia
@testset "Integration with cpolyopt" begin
    sys = fermionic_algebra(2)
    c, c_dag = sys.variables

    # Number operator: n₁ = c₁† c₁
    n1 = ComplexF64(1.0) * c_dag[1] * c[1]

    # Create optimization problem using algebra interface
    pop = cpolyopt(n1, sys)

    @test pop.objective == n1
    @test pop.is_unipotent == false
    @test pop.is_projective == false

    # Check that algebra constraints were included
    @test length(pop.eq_constraints) == length(sys.equality_constraints)

    # Verify fermionic constraints are present
    @test c[1] * c[1] in pop.eq_constraints  # c₁² = 0
end
```

#### Test 3.7: Additional User Constraints

```julia
@testset "Custom Constraints" begin
    sys = fermionic_algebra(2)
    c, c_dag = sys.variables

    # Number operators
    n1 = c_dag[1] * c[1]
    n2 = c_dag[2] * c[2]

    # Hamiltonian
    ham = ComplexF64(1.0) * n1 * n2

    # Add custom constraint: total particle number = 1
    N_total = n1 + n2 - 1

    pop = cpolyopt(ham, sys; eq_constraints=[N_total])

    # Should have algebra constraints + custom constraint
    @test length(pop.eq_constraints) == length(sys.equality_constraints) + 1
    @test N_total in pop.eq_constraints

    # Fermionic constraints still present
    @test c[1] * c[1] in pop.eq_constraints
end
```

#### Test 3.8: Error Handling

```julia
@testset "Input Validation" begin
    @test_throws AssertionError fermionic_algebra(0)
    @test_throws AssertionError fermionic_algebra(-1)
end
```

#### Test 3.9: Integration Tests with Solver (Optional)

**Run only with LOCAL_TESTING environment variable**

```julia
if haskey(ENV, "LOCAL_TESTING")
    @testset "Integration: Free Fermion Ground State" begin
        using MosekTools  # or Clarabel

        N = 2
        sys = fermionic_algebra(N)
        c, c_dag = sys.variables

        # Single-particle energies
        ε1, ε2 = -1.0, -0.5

        # Hamiltonian: H = ε₁n₁ + ε₂n₂
        n1 = c_dag[1] * c[1]
        n2 = c_dag[2] * c[2]
        ham = ComplexF64(ε1) * n1 + ComplexF64(ε2) * n2

        pop = cpolyopt(ham, sys)

        # Solve
        config = SolverConfig(order=1, optimizer=Mosek.Optimizer)
        result = cs_nctssos(pop, config)

        # Ground state should be vacuum (no particles): E = 0
        @test isapprox(result.objective, 0.0, atol=1e-6)
    end

    @testset "Integration: Single Particle Ground State" begin
        N = 2
        sys = fermionic_algebra(N)
        c, c_dag = sys.variables

        ε1, ε2 = -1.0, -0.5

        n1 = c_dag[1] * c[1]
        n2 = c_dag[2] * c[2]
        ham = ComplexF64(ε1) * n1 + ComplexF64(ε2) * n2

        # Constraint: exactly 1 particle
        N_total = n1 + n2 - 1
        pop = cpolyopt(ham, sys; eq_constraints=[N_total])

        config = SolverConfig(order=2, optimizer=Mosek.Optimizer)
        result = cs_nctssos(pop, config)

        # Ground state: particle in lowest energy mode
        # E = min(ε₁, ε₂) = -1.0
        @test isapprox(result.objective, ε1, atol=1e-4)
    end
end
```

---

### Phase 4: Documentation

**File**: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/fermionic_algebra_interface.jl`

**Action**: Create tutorial following pauli_algebra_interface.jl pattern

**Structure**:
1. Introduction to fermionic systems
2. Manual approach (show complexity)
3. Simplified approach with fermionic_algebra
4. Example: Free fermion system
5. Example: Fermi-Hubbard model
6. Advantages of fermionic_algebra
7. Common observables
8. References

**Key Points to Document**:
- All properties enforced via constraints (not simplification)
- Nilpotent constraints (c² = 0) are redundant but included for clarity
- Constraint count scales as O(N²) for N modes
- Compatible with cpolyopt(objective, algebra) interface
- Easy to add custom constraints

**Template** (abbreviated):
```julia
# # Fermionic Algebra Interface
#
# This tutorial demonstrates the `fermionic_algebra` constructor for
# working with fermionic quantum systems.
#
# ## Introduction
#
# Fermionic operators describe indistinguishable particles that obey
# the Pauli exclusion principle (e.g., electrons, quarks).
#
# ## The Simplified Approach

using NCTSSoS

# Create 2-mode fermionic system
sys = fermionic_algebra(2)
c, c_dag = sys.variables

# Number operators
n1 = c_dag[1] * c[1]
n2 = c_dag[2] * c[2]

# Free fermion Hamiltonian
ham = ComplexF64(-1.0) * n1 + ComplexF64(-0.5) * n2

# Create optimization problem
pop = cpolyopt(ham, sys)

# ## Example: Fermi-Hubbard Model
# [Include numerical example with solver]

# ## Advantages
# 1. Concise: 1 line vs 20+ lines of manual setup
# 2. Correct by construction: All constraints automatically included
# 3. Scales easily: Works for any N
# 4. Maintainable: Centralized algebra definition
```

---

## 3. Test-Driven Development Plan

### TDD Cycle

For each test section:
1. **Write failing test** for specific functionality
2. **Implement minimum code** to pass the test
3. **Explain what changed** and why
4. **Refactor** if needed
5. **Move to next test**

### Test Execution Order

1. **Test 3.1**: Basic structure N=1 (validates constructor works)
2. **Test 3.5**: SimplifyAlgorithm properties (validates setup)
3. **Test 3.3**: Constraint verification N=1 (validates constraint generation)
4. **Test 3.2**: Constraint scaling (validates generalization to N=2,3)
5. **Test 3.4**: Constraint verification N=2 (validates specific constraints)
6. **Test 3.6**: Integration with cpolyopt (validates interface)
7. **Test 3.7**: Custom constraints (validates composability)
8. **Test 3.8**: Error handling (validates input validation)
9. **Test 3.9**: Solver integration (validates end-to-end, if LOCAL_TESTING)

---

## 4. Implementation Checklist

### Phase 1: Core Implementation
- [ ] Implement `fermionic_algebra(N)` function in `src/algebra_constructors.jl`
- [ ] Add comprehensive docstring with examples
- [ ] Verify constraint generation for N=1,2,3

### Phase 2: Module Integration
- [ ] Export `fermionic_algebra` in `src/NCTSSoS.jl`
- [ ] Verify function is accessible after `using NCTSSoS`

### Phase 3: Testing
- [ ] Write and pass Test 3.1 (Basic structure N=1)
- [ ] Write and pass Test 3.5 (SimplifyAlgorithm properties)
- [ ] Write and pass Test 3.3 (Constraint verification N=1)
- [ ] Write and pass Test 3.2 (Constraint scaling)
- [ ] Write and pass Test 3.4 (Constraint verification N=2)
- [ ] Write and pass Test 3.6 (cpolyopt integration)
- [ ] Write and pass Test 3.7 (Custom constraints)
- [ ] Write and pass Test 3.8 (Error handling)
- [ ] (Optional) Write and pass Test 3.9 (Solver integration)
- [ ] Verify all tests pass
- [ ] Verify no regressions in existing tests

### Phase 4: Documentation
- [ ] Create `fermionic_algebra_interface.jl` tutorial
- [ ] Include physical motivation
- [ ] Show manual vs simplified approach
- [ ] Add free fermion example
- [ ] Add Fermi-Hubbard example (with numerical solution if possible)
- [ ] List advantages
- [ ] Add cross-references
- [ ] Update documentation index
- [ ] Build docs locally and verify

---

## 5. Expected Outcomes

After implementation:
- ✅ Users can create fermionic systems with `fermionic_algebra(N)`
- ✅ All fermionic properties enforced via equality constraints
- ✅ Nilpotent constraints (c² = 0) explicitly included
- ✅ Compatible with cpolyopt(objective, algebra) interface
- ✅ Easy to add custom constraints
- ✅ No modifications to FastPolynomials, PolyOpt, or ComplexPolyOpt
- ✅ Comprehensive test coverage
- ✅ Complete documentation with examples
- ✅ No regressions in existing functionality

---

## 6. Advantages of Constraint-Based Approach

### Simplicity
- **Implementation**: ~40 lines vs 500+ lines in previous plan
- **No core modifications**: Existing types unchanged
- **No workarounds**: Zero monomial issue irrelevant

### Correctness
- **Explicit constraints**: All algebra properties visible
- **SDP solver enforcement**: Constraints guaranteed during optimization
- **No simplification bugs**: No complex monomial manipulation

### Maintainability
- **Isolated code**: All fermionic logic in one function
- **Easy to extend**: Add new constraint types easily
- **Future-proof**: No dependencies on internal simplification logic

### Trade-offs
- **More constraints**: 2N additional constraints for nilpotency
  - But: Constraints are cheap for SDP solvers
  - Impact: Negligible for N < 10
- **No automatic simplification**: c² terms not simplified to 0
  - But: Constraint enforcement achieves same result
  - Impact: None on optimization results

---

## 7. Constraint Count Analysis

### Scaling

| N | Anti-comm | Canonical | Nilpotent | Total | O(N) Formula |
|---|-----------|-----------|-----------|-------|--------------|
| 1 | 2 | 1 | 2 | 5 | 2N² + 3N |
| 2 | 6 | 4 | 4 | 14 | 2N² + 3N |
| 3 | 12 | 9 | 6 | 27 | 2N² + 3N |
| 5 | 30 | 25 | 10 | 65 | 2N² + 3N |
| 10 | 110 | 100 | 20 | 230 | 2N² + 3N |

### Performance Considerations

**Recommended**: N ≤ 10 for practical computations
- Constraint count: ~230 for N=10
- SDP matrix size: Depends on relaxation order
- Typical relaxation order: 1-2 for fermionic systems

**Large N**: For N > 10, consider:
- Symmetry reduction (future work)
- Sparse constraint representation
- Problem-specific simplifications

---

## 8. Comparison with Pauli Algebra

| Aspect | Pauli | Fermionic |
|--------|-------|-----------|
| Variables | σx, σy, σz per site | c, c† per mode |
| Commutation | Anti-commute at same site | All anti-commute |
| Square property | σ² = 1 (unipotent) | c² = 0 (nilpotent) |
| Constraints/site | 6 | 2N + 3 (total formula) |
| Scaling | O(N) | O(N²) |
| SimplifyAlgorithm | is_unipotent=true | No special flags |
| Constraint encoding | Commutation relations | Anti-comm + nilpotent |

---

## 9. Implementation Notes

### Redundancy of Nilpotent Constraints

**Mathematical Note**: The nilpotent constraints c² = 0 are technically redundant:
- From {c,c} = 2c² = 0, we get c² = 0
- Similarly for c†

**Why Include Them?**:
1. **Explicitness**: Makes nilpotency clear in constraint set
2. **Numerical stability**: Some SDP solvers may benefit
3. **Documentation**: Self-documenting constraint structure
4. **Cost**: Only 2N additional constraints (negligible)

**Design Decision**: Include for clarity and robustness, despite redundancy.

### Variable Naming

**Choice**: `c_dag` instead of `c†`
- Julia identifiers don't support superscripts directly
- `c_dag` is clear and conventional in physics code
- Alternative: `cdag` (less clear), `c_dagger` (verbose)

**Recommendation**: Use `c_dag` throughout for consistency.

### Commutation Group Structure

**Single Group**: All operators in `[c[1], c_dag[1], c[2], c_dag[2], ...]`

**Implications**:
- No automatic reordering by SimplifyAlgorithm
- All algebraic relations via explicit constraints
- Matches the constraint-based philosophy

---

## 10. Future Extensions (Not in Current Scope)

### Symmetry Reduction
- Exploit particle number conservation
- Reduce constraint count for large N
- Requires analysis of problem structure

### Multiple Species
- Spin-up and spin-down fermions
- Requires additional variable sets
- More complex constraint structure

### Fermionic-Bosonic Mixed Systems
- Combine with bosonic operators
- Requires careful commutation structure
- Advanced use case

### Majorana Fermions
- Self-adjoint fermionic operators
- Different constraint structure
- Special case of fermionic algebra

---

## 11. Summary for Parent Agent

### Implementation Steps

1. **Add fermionic_algebra function** to `src/algebra_constructors.jl` (~40 lines)
2. **Export fermionic_algebra** in `src/NCTSSoS.jl` (1 line)
3. **Write comprehensive tests** in `test/algebra_constructors.jl` (~150 lines)
4. **Create tutorial documentation** (~100 lines)

**Total**: ~300 lines of new code, 0 lines of modifications to existing code

### Key Success Criteria

✅ **No core modifications**: FastPolynomials, PolyOpt, ComplexPolyOpt unchanged
✅ **Constraint-based**: All properties via equality constraints
✅ **Pattern matching**: Follows pauli_algebra interface exactly
✅ **TDD approach**: Write tests first, implement second
✅ **Documentation**: Clear examples and explanations

### Critical Differences from Previous Plan

| Aspect | Previous | Revised |
|--------|----------|---------|
| Core changes | Yes (SimplifyAlgorithm, PolyOpt) | No |
| Nilpotency | Simplification rule | Equality constraint |
| Zero monomials | Workaround needed | Not needed |
| Implementation size | ~500 lines changes | ~40 lines new code |
| Complexity | High | Low |
| Maintenance | Moderate | Easy |

---

**End of Revised Plan**

This plan is complete and ready for parent agent implementation via TDD.
The approach is dramatically simpler than the previous plan while achieving
the same user-facing functionality.
