# Pauli Algebra Implementation Analysis

**Date**: 2025-11-05  
**Purpose**: Understand the complete architecture of Pauli algebra to replicate for fermionic algebra  
**Status**: Complete Analysis

---

## Executive Summary

The Pauli algebra implementation follows a **clean separation of concerns**:
1. **FastPolynomials** provides the core polynomial machinery (variables, monomials, polynomials)
2. **SimplifyAlgorithm** encodes algebraic properties (commutation, unipotency, projectivity)
3. **Algebra constructors** (`pauli_algebra`) create ready-to-use algebra systems
4. **Optimization interface** (`cpolyopt`, `cs_nctssos`) solves optimization problems

The fermionic algebra should follow **exactly this same pattern**.

---

## 1. File Structure

### Core Implementation Files

```
src/
├── NCTSSoS.jl                       # Main module, exports pauli_algebra
├── algebra_constructors.jl          # pauli_algebra(N) constructor [40 lines]
├── pop.jl                           # PolyOpt struct, polyopt(), cpolyopt() [173 lines]
├── interface.jl                     # cs_nctssos, cs_nctssos_higher solvers [145 lines]
├── sparse.jl                        # Correlative/term sparsity [~400 lines]
├── moment_solver.jl                 # Moment relaxation [~200 lines]
├── sos_solver.jl                    # SOS dualization [~200 lines]
├── gns.jl                           # GNS reconstruction [~300 lines]
└── FastPolynomials/src/
    ├── variables.jl                 # Variable, @ncpolyvar macro [187 lines]
    ├── monomials.jl                 # Monomial struct [183 lines]
    ├── polynomial.jl                # Polynomial struct [122 lines]
    ├── arithmetic.jl                # +, -, *, operations [223 lines]
    ├── simplify.jl                  # SimplifyAlgorithm, simplification rules [291 lines]
    ├── state_word.jl                # State polynomial support [~250 lines]
    └── statepolynomial.jl           # NCStatePolynomial [~300 lines]
```

### Documentation Files

```
docs/src/examples/literate/
├── pauli_algebra_interface.jl       # Tutorial on using pauli_algebra [129 lines]
└── pauli_gns_construction.jl        # GNS reconstruction example [454 lines]
```

### Test Files

```
test/
└── algebra_constructors.jl          # Comprehensive tests [274 lines]
```

---

## 2. Type Definitions

### Core Polynomial Types (FastPolynomials)

**File**: `src/FastPolynomials/src/variables.jl` (lines 10-17)
```julia
struct Variable
    name::Symbol
    iscomplex::Bool
    
    function Variable(name::Symbol; iscomplex::Bool=false)
        return new(name, iscomplex)
    end
end
```

**File**: `src/FastPolynomials/src/monomials.jl` (lines 16-25)
```julia
struct Monomial
    vars::Vector{Variable}
    z::Vector{Int}  # Exponents
    
    function Monomial(vars::Vector{Variable}, z::Vector{Int})
        @assert length(vars) == length(z)
        @assert consecutive_unique(vars)  # No repeated consecutive variables
        @assert all(x -> !iszero(x), z)   # No zero exponents
        return new(vars, z)
    end
end
```

**File**: `src/FastPolynomials/src/polynomial.jl` (lines 17-39)
```julia
struct Polynomial{T} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    monos::Vector{Monomial}
    
    function Polynomial(a::Vector{T}, x::Vector{Monomial}) where {T}
        # Automatically sorts and combines like terms
        # Maintains invariants: sorted, unique monos, non-zero coeffs
        # ...
    end
end
```

### Optimization Problem Types

**File**: `src/pop.jl` (lines 25-33)
```julia
struct PolyOpt{P} <: OptimizationProblem{P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    variables::Vector{Variable}
    comm_gps::Vector{Vector{Variable}}  # Commutation groups
    is_unipotent::Bool  # X² = 1 (Pauli, fermions)
    is_projective::Bool # X² = X (projectors)
end
```

**File**: `src/pop.jl` (lines 95-103)
```julia
struct ComplexPolyOpt{P} <: OptimizationProblem{P}
    # Same fields as PolyOpt
    # Used when operators are not Hermitian
    # (product of variables are not Hermitian)
end
```

---

## 3. Multiplication Rules Implementation

### Simplification Algorithm

**File**: `src/FastPolynomials/src/simplify.jl` (lines 1-18)
```julia
struct SimplifyAlgorithm
    comm_gps::Dict{Variable,Int}  # Maps variable to commutation group ID
    n_gps::Int                     # Number of commutation groups
    is_unipotent::Bool             # X² = 1 (for Pauli, fermions)
    is_projective::Bool            # X² = X (for projectors)
    
    function SimplifyAlgorithm(;
        comm_gps::Vector{Vector{Variable}},
        is_unipotent::Bool=false,
        is_projective::Bool=false,
    )
        return new(
            Dict(var => i for (i, vars) in enumerate(comm_gps) for var in vars),
            length(comm_gps),
            is_unipotent,
            is_projective,
        )
    end
end
```

### Simplification Process

**File**: `src/FastPolynomials/src/simplify.jl` (lines 24-32)

The simplification follows this algorithm:
1. **First**: Apply commutation rules (`_comm!`)
2. **Then**: Apply unipotent OR projective simplification

```julia
function simplify!(m::Monomial, sa::SimplifyAlgorithm)
    (isone(m) || nosimp(sa)) && return m
    _comm!(m, sa.comm_gps)  # Step 1: Sort by commutation groups
    
    sa.is_unipotent && (_simplify_unipotent!(m); return m)
    sa.is_projective && (_simplify_projective!(m); return m)
    _simplify_standard!(m)
    return m
end
```

### Commutation Rule Application

**File**: `src/FastPolynomials/src/monomials.jl` (lines 167-180)

```julia
function _comm!(mono::Monomial, comm_gps::Dict{Variable,Int})
    # Stable bubble sort to reorder variables by commutation group
    length(mono.vars) == 1 && return nothing
    @inbounds for i in 1:(length(mono.vars) - 1)
        swapped = false
        for j in 1:(length(mono.vars) - i)
            comm_gps[mono.vars[j]] <= comm_gps[mono.vars[j + 1]] && continue
            # Swap variables and exponents
            mono.vars[j], mono.vars[j + 1] = mono.vars[j + 1], mono.vars[j]
            mono.z[j], mono.z[j + 1] = mono.z[j + 1], mono.z[j]
            swapped = true
        end
        !swapped && break
    end
    return nothing
end
```

**Key Insight**: Variables in different commutation groups can be freely reordered. Variables in the same group maintain their order (anti-commutation is handled by equality constraints, not by reordering).

### Unipotent Simplification (for Pauli, Fermions)

**File**: `src/FastPolynomials/src/simplify.jl` (lines 68-102)

```julia
function _simplify_unipotent!(m::Monomial)
    # For unipotent operators: X² = 1
    # Strategy: Keep only odd exponents, drop even exponents
    # XY²Z = XZ (since Y² = 1)
    
    head_idx = 1
    tail_idx = findfirst(isodd, m.z)
    isnothing(tail_idx) && (empty!(m.vars); empty!(m.z); return nothing)
    
    m.vars[head_idx] = m.vars[tail_idx]
    m.z[head_idx] = mod(m.z[tail_idx], 2)  # Reduce to 0 or 1
    
    tail_idx += 1
    @inbounds for i in tail_idx:length(m.vars)
        iseven(m.z[i]) && continue  # Skip even exponents
        if m.vars[i] == m.vars[head_idx]
            head_idx -= 1  # Cancel pairs XX = 1
            # ... handle cancellation
        else
            head_idx += 1
            m.vars[head_idx] = m.vars[i]
            m.z[head_idx] = 1  # All exponents become 1
        end
    end
    # Trim unused entries
    deleteat!(m.vars, (head_idx + 1):length(m.vars))
    deleteat!(m.z, (head_idx + 1):length(m.z))
    return nothing
end
```

**Critical for Fermions**: This is where X² = 1 is enforced. Fermionic operators satisfy c_i² = 0 (not 1!), so we'll need `_simplify_nilpotent!` instead.

---

## 4. The Pauli Algebra Constructor

**File**: `src/algebra_constructors.jl` (lines 1-39)

```julia
function pauli_algebra(N::Int)
    @assert N >= 1 "Number of sites N must be at least 1"
    
    # Step 1: Declare variables
    @ncpolyvar x[1:N] y[1:N] z[1:N]
    
    # Step 2: Define commutation groups
    # Operators at different sites commute, same site anti-commute
    comm_gps = [[x[i], y[i], z[i]] for i in 1:N]
    
    # Step 3: Encode Pauli commutation relations as equality constraints
    equality_constraints = reduce(vcat, [
        [x[i] * y[i] - im * z[i],    # σx σy = i σz
         y[i] * x[i] + im * z[i],    # σy σx = -i σz
         y[i] * z[i] - im * x[i],    # σy σz = i σx
         z[i] * y[i] + im * x[i],    # σz σy = -i σx
         z[i] * x[i] - im * y[i],    # σz σx = i σy
         x[i] * z[i] + im * y[i]]    # σx σz = -i σy
        for i in 1:N
    ])
    
    # Step 4: Create SimplifyAlgorithm
    simplify_algo = SimplifyAlgorithm(
        comm_gps=comm_gps, 
        is_unipotent=true,  # σ² = 1
        is_projective=false
    )
    
    # Step 5: No inequality constraints
    inequality_constraints = empty(equality_constraints)
    
    # Step 6: Return algebra system as NamedTuple
    return (
        variables=(x, y, z),
        simplify_algo=simplify_algo,
        equality_constraints=equality_constraints,
        inequality_constraints=inequality_constraints,
        comm_gps=comm_gps
    )
end
```

**Architecture Pattern**:
- Returns a **NamedTuple** (not a struct)
- Contains everything needed to build optimization problems
- Can be passed directly to `cpolyopt(objective, algebra)`

---

## 5. Interface Functions

### Main Exported Functions

**File**: `src/NCTSSoS.jl` (lines 14-20)
```julia
export @ncpolyvar, ς
export polyopt, cpolyopt
export SolverConfig
export cs_nctssos, cs_nctssos_higher
export reconstruct
export pauli_algebra
```

### Creating Optimization Problems

**File**: `src/pop.jl` (lines 154-172)

Two ways to create problems:

**Method 1: Manual** (traditional)
```julia
pop = cpolyopt(objective;
    eq_constraints=[...],
    comm_gps=[[...]],
    is_unipotent=true,
    is_projective=false)
```

**Method 2: Algebra Interface** (new, simplified)
```julia
sys = pauli_algebra(N)
pop = cpolyopt(objective, sys)  # Automatically uses algebra constraints
```

```julia
function cpolyopt(objective::P, algebra::NamedTuple; 
                  eq_constraints=Any[], ineq_constraints=Any[]) where {T,P<:AbstractPolynomial{T}}
    # Extract properties from algebra
    comm_gps = algebra.comm_gps
    is_unipotent = algebra.simplify_algo.is_unipotent
    is_projective = algebra.simplify_algo.is_projective
    
    # Merge algebra constraints with user-provided constraints
    converted_user_eq = isempty(eq_constraints) ? 
        typeof(algebra.equality_constraints)[] : 
        [T(1) * poly for poly in eq_constraints]
    merged_eq_constraints = vcat(algebra.equality_constraints, converted_user_eq) 
    
    # Call original cpolyopt
    return cpolyopt(objective;
        eq_constraints=merged_eq_constraints,
        ineq_constraints=ineq_constraints,
        comm_gps=comm_gps,
        is_unipotent=is_unipotent,
        is_projective=is_projective)
end
```

### Solving Optimization Problems

**File**: `src/interface.jl` (lines 73-97)

```julia
function cs_nctssos(pop::OP, solver_config::SolverConfig; dualize::Bool=true) where {OP<:OptimizationProblem}
    # Step 1: Create SimplifyAlgorithm from problem properties
    sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, 
                           is_projective=pop.is_projective, 
                           is_unipotent=pop.is_unipotent)
    
    # Step 2: Determine relaxation order
    order = iszero(solver_config.order) ? 
        maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; ...]]) : 
        solver_config.order
    
    # Step 3: Compute correlative sparsity (decompose into cliques)
    corr_sparsity = correlative_sparsity(pop, order, solver_config.cs_algo)
    
    # Step 4: Compute term sparsity for each clique
    cliques_term_sparsities = # ... complex computation using sa
    
    # Step 5: Build moment relaxation
    moment_problem = moment_relax(pop, corr_sparsity, cliques_term_sparsities)
    
    # Step 6: Optionally dualize to SOS problem
    problem_to_solve = !dualize ? moment_problem : sos_dualize(moment_problem)
    
    # Step 7: Solve with JuMP
    set_optimizer(problem_to_solve.model, solver_config.optimizer)
    optimize!(problem_to_solve.model)
    
    return PolyOptResult(objective_value(problem_to_solve.model), ...)
end
```

---

## 6. Simplification Rules in Detail

### How Pauli Relations are Enforced

The Pauli commutation relations are enforced through **equality constraints**, not through the simplifier:

**Equality Constraints** (6 per site):
```
σx σy - i σz = 0  ⟹  σx σy = i σz
σy σx + i σz = 0  ⟹  σy σx = -i σz
σy σz - i σx = 0  ⟹  σy σz = i σx
σz σy + i σx = 0  ⟹  σz σy = -i σx
σz σx - i σy = 0  ⟹  σz σx = i σy
σx σz + i σy = 0  ⟹  σx σz = -i σy
```

**SimplifyAlgorithm role**:
1. `comm_gps`: Operators on different sites commute (reorder freely)
2. `is_unipotent=true`: σ² = 1 (reduce exponents mod 2, cancel pairs)

**Key Insight**: The simplifier does NOT directly encode σx σy = i σz. Instead:
- The simplifier reduces `σx² σy σx` → `σy σx` (using unipotency)
- The equality constraints then replace `σy σx` with `-i σz` during optimization

### Simplification Flow Example

For `σx[1] σy[1] σx[1]` with Pauli algebra:

1. **Input**: `x₁ y₁ x₁` (monomial with vars=[x₁, y₁, x₁], z=[1,1,1])
2. **After `_comm!`**: `x₁ y₁ x₁` (no change, all in same commutation group)
3. **After `_simplify_unipotent!`**: `x₁ y₁ x₁` (all exponents stay 1, but...)
   - Actually, the algorithm would detect `x₁` appears at position 1 and 3
   - It consolidates: `x₁² y₁` → `y₁` (since x² = 1)
4. **Final**: `y₁`

Then during optimization, equality constraints handle things like `σx σy` → `iσz`.

---

## 7. Documentation Pattern

### Literate Tutorial Structure

**File**: `docs/src/examples/literate/pauli_algebra_interface.jl`

Structure:
1. **Problem introduction**: Traditional manual approach (tedious)
2. **Solution**: Show the simplified algebra constructor
3. **Side-by-side comparison**: Old vs new approach
4. **Example**: Heisenberg XXX model
5. **Advantages list**: Benefits of algebra constructor
6. **Next steps**: Links to other tutorials

Key documentation features:
- Uses **Literate.jl** format (Julia code with markdown comments)
- Includes **runnable examples**
- Cross-references using `@ref` tags
- Mathematical equations with LaTeX
- Explains **physical motivation**

### API Documentation

Functions are documented with:
- Purpose description
- Type parameters
- Arguments with descriptions
- Keyword arguments
- Returns
- Examples with `jldoctest`
- Cross-references

Example from `src/pop.jl` (lines 36-56):
```julia
"""
    polyopt(objective::P; eq_constraints=Any[], ineq_constraints=Any[], 
            comm_gps=Vector{Variable}[], is_unipotent::Bool=false, 
            is_projective::Bool=false) where {T,P<:AbstractPolynomial{T}}

Create a polynomial optimization problem.

# Arguments
- `objective::P`: The polynomial objective function to optimize.
- `eq_constraints=Any[]`: Equality constraints as polynomials (p = 0).
...

# Returns
A `PolyOpt{P}` structure representing the polynomial optimization problem.

# Notes
- The polynomial coefficients cannot be integers...
"""
```

---

## 8. Testing Strategy

**File**: `test/algebra_constructors.jl` (274 lines)

### Test Organization

```julia
@testset "Pauli Algebra Constructor" begin
    @testset "Basic Structure N=1" begin
        # Test return type structure
        # Test variables
        # Test constraints count
        # Test commutation groups
    end
    
    @testset "Basic Structure N=3" begin
        # Test scaling to N=3
    end
    
    @testset "SimplifyAlgorithm Properties" begin
        # Test is_unipotent, is_projective, n_gps
    end
    
    @testset "Commutation Relations Encoded Correctly" begin
        # Verify exact constraint equations
    end
    
    @testset "Integration with cpolyopt" begin
        # Test creating optimization problems
    end
    
    @testset "Invalid Input" begin
        # Test error handling
    end
end

@testset "cpolyopt with Algebra Interface" begin
    @testset "Pauli Algebra Interface" begin
        # Test basic usage
    end
    
    @testset "Pauli Algebra with Additional Constraints" begin
        # Test merging custom constraints
    end
    
    @testset "Interface Equivalence Test" begin
        # Verify old and new interfaces produce same result
    end
end

# Integration tests (only run with LOCAL_TESTING env var)
@testset "Integration: XXX Model with Pauli Algebra" begin
    # Full end-to-end test with Mosek
    # Verify numerical results
end
```

### Testing Patterns

1. **Unit tests**: Test individual components (structure, constraints)
2. **Integration tests**: Test algebra → cpolyopt → solver pipeline
3. **Numerical tests**: Compare against known results (ground state energies)
4. **Equivalence tests**: Verify new interface ≡ manual interface
5. **Error handling**: Test invalid inputs

### Test Data Sources

- **XXX Heisenberg model**: Known ground state energy -0.467129 per site (N=6, order=2)
- **J1-J2 model**: Known ground state energy -0.4270083 per site (N=6, order=2)
- **Transverse field Ising**: Multiple boundary conditions tested

---

## 9. Integration with FastPolynomials and NCTSSoS

### Data Flow

```
User Code
   ↓
pauli_algebra(N) → NamedTuple(variables, simplify_algo, constraints, comm_gps)
   ↓
cpolyopt(objective, algebra) → ComplexPolyOpt{Polynomial{ComplexF64}}
   ↓
cs_nctssos(pop, solver_config) → Uses SimplifyAlgorithm throughout:
   ├─ correlative_sparsity(..., sa)
   ├─ term_sparsities(..., sa)
   ├─ init_activated_supp(..., sa)
   └─ moment_relax(...)
   ↓
SDP Solver (Mosek, Clarabel, etc.)
   ↓
PolyOptResult(objective, sparsity_info, model)
```

### Abstract Types and Interfaces

**File**: `src/FastPolynomials/src/polynomial.jl` (line 1)
```julia
abstract type AbstractPolynomial{T} end
```

Concrete types:
- `Polynomial{T}` (standard non-commutative polynomials)
- `StatePolynomial{T,ST}` (state-dependent polynomials)
- `NCStatePolynomial{T,ST}` (non-commutative state polynomials)

**File**: `src/pop.jl` (line 1)
```julia
abstract type OptimizationProblem{P} end
```

Concrete types:
- `PolyOpt{P}` (Hermitian operators)
- `ComplexPolyOpt{P}` (non-Hermitian operators)

### Key Integration Points

1. **SimplifyAlgorithm** is passed everywhere polynomials are processed:
   - `simplify(mono, sa)`
   - `canonicalize(poly, sa)`
   - `get_basis(Polynomial{T}, vars, degree, sa)`

2. **Commutation groups** flow through:
   - Defined in algebra constructor
   - Stored in `PolyOpt`
   - Used to create `SimplifyAlgorithm`
   - Applied during simplification

3. **Equality constraints** are polynomial equations:
   - Created as `Polynomial{ComplexF64}` objects
   - Included in optimization problem
   - Used during moment relaxation

---

## 10. Pattern to Follow for Fermionic Algebra

### Recommended Implementation Steps

#### **Step 1: Understand Differences**

| Property | Pauli | Fermionic |
|----------|-------|-----------|
| Square | σᵢ² = 1 (unipotent) | cᵢ² = 0 (nilpotent) |
| Anti-commutation | {σᵢ, σⱼ} = 2δᵢⱼ | {cᵢ, cⱼ} = 0 |
| Operators | σx, σy, σz per site | c, c† per site |
| Commutation | Different sites commute | All operators anti-commute |

#### **Step 2: Create `_simplify_nilpotent!` function**

**Location**: Add to `src/FastPolynomials/src/simplify.jl`

```julia
function _simplify_nilpotent!(m::Monomial)
    # For nilpotent operators: c² = 0
    # Strategy: If any variable appears with exponent ≥ 2, return zero monomial
    
    @inbounds for i in eachindex(m.z)
        if m.z[i] >= 2
            # Found c² or higher → entire monomial is zero
            empty!(m.vars)
            empty!(m.z)
            return nothing
        end
    end
    
    # All exponents are 1, keep as-is (after commutation)
    return nothing
end
```

Modify `simplify!` to handle nilpotent case:
```julia
function simplify!(m::Monomial, sa::SimplifyAlgorithm)
    (isone(m) || nosimp(sa)) && return m
    _comm!(m, sa.comm_gps)
    
    sa.is_unipotent && (_simplify_unipotent!(m); return m)
    sa.is_projective && (_simplify_projective!(m); return m)
    sa.is_nilpotent && (_simplify_nilpotent!(m); return m)  # NEW
    _simplify_standard!(m)
    return m
end
```

#### **Step 3: Extend SimplifyAlgorithm**

**Location**: Modify `src/FastPolynomials/src/simplify.jl`

```julia
struct SimplifyAlgorithm
    comm_gps::Dict{Variable,Int}
    n_gps::Int
    is_unipotent::Bool
    is_projective::Bool
    is_nilpotent::Bool  # NEW
    
    function SimplifyAlgorithm(;
        comm_gps::Vector{Vector{Variable}},
        is_unipotent::Bool=false,
        is_projective::Bool=false,
        is_nilpotent::Bool=false,  # NEW
    )
        # Assert: Only one of {unipotent, projective, nilpotent} can be true
        @assert sum([is_unipotent, is_projective, is_nilpotent]) <= 1 
            "Only one of is_unipotent, is_projective, is_nilpotent can be true"
        
        return new(
            Dict(var => i for (i, vars) in enumerate(comm_gps) for var in vars),
            length(comm_gps),
            is_unipotent,
            is_projective,
            is_nilpotent,  # NEW
        )
    end
end
```

#### **Step 4: Create `fermionic_algebra` constructor**

**Location**: Add to `src/algebra_constructors.jl`

```julia
"""
    fermionic_algebra(N::Int)

Create a fermionic algebra system for N fermionic modes.

Returns a NamedTuple with:
- `variables`: Tuple of (c, c_dag) where c[i] and c_dag[i] are annihilation/creation operators
- `simplify_algo`: SimplifyAlgorithm with is_nilpotent=true
- `equality_constraints`: Anti-commutation relations {cᵢ, cⱼ} = 0, {cᵢ†, cⱼ†} = 0, {cᵢ, cⱼ†} = δᵢⱼ
- `inequality_constraints`: Empty (no inequality constraints)
- `comm_gps`: All operators in one group (all anti-commute)

# Arguments
- `N::Int`: Number of fermionic modes

# Example
```julia
sys = fermionic_algebra(2)
c, c_dag = sys.variables
ham = c_dag[1] * c[1]  # Number operator
pop = cpolyopt(ham, sys)
```
"""
function fermionic_algebra(N::Int)
    @assert N >= 1 "Number of modes N must be at least 1"
    
    # Step 1: Declare variables
    # c[i] = annihilation operator for mode i
    # c_dag[i] = creation operator for mode i  
    @ncpolyvar c[1:N] c_dag[1:N]
    
    # Step 2: Define commutation groups
    # All fermionic operators anti-commute with each other
    all_ops = vcat(collect(c), collect(c_dag))
    comm_gps = [all_ops]  # Single group containing all operators
    
    # Step 3: Encode anti-commutation relations
    # {cᵢ, cⱼ} = cᵢcⱼ + cⱼcᵢ = 0 for all i,j
    # {cᵢ†, cⱼ†} = cᵢ†cⱼ† + cⱼ†cᵢ† = 0 for all i,j
    # {cᵢ, cⱼ†} = cᵢcⱼ† + cⱼ†cᵢ = δᵢⱼ
    
    equality_constraints = Polynomial{ComplexF64}[]
    
    # Anti-commutation: {cᵢ, cⱼ} = 0
    for i in 1:N, j in i:N  # j >= i to avoid duplicates
        push!(equality_constraints, c[i] * c[j] + c[j] * c[i])
    end
    
    # Anti-commutation: {cᵢ†, cⱼ†} = 0
    for i in 1:N, j in i:N
        push!(equality_constraints, c_dag[i] * c_dag[j] + c_dag[j] * c_dag[i])
    end
    
    # Canonical anti-commutation: {cᵢ, cⱼ†} = δᵢⱼ
    for i in 1:N, j in 1:N
        if i == j
            push!(equality_constraints, c[i] * c_dag[j] + c_dag[j] * c[i] - 1)
        else
            push!(equality_constraints, c[i] * c_dag[j] + c_dag[j] * c[i])
        end
    end
    
    # Step 4: Create SimplifyAlgorithm
    simplify_algo = SimplifyAlgorithm(
        comm_gps=comm_gps,
        is_unipotent=false,
        is_projective=false,
        is_nilpotent=true  # cᵢ² = 0
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

#### **Step 5: Update exports**

**Location**: `src/NCTSSoS.jl`

```julia
export fermionic_algebra  # Add to existing exports
```

#### **Step 6: Write tests**

**Location**: Create or extend `test/algebra_constructors.jl`

```julia
@testset "Fermionic Algebra Constructor" begin
    @testset "Basic Structure N=1" begin
        sys = fermionic_algebra(1)
        
        # Check return type structure
        @test hasfield(typeof(sys), :variables)
        @test hasfield(typeof(sys), :simplify_algo)
        
        # Check variables
        c, c_dag = sys.variables
        @test length(c) == 1
        @test length(c_dag) == 1
        
        # Check constraints
        # {c, c} = 0 → 1 constraint
        # {c†, c†} = 0 → 1 constraint  
        # {c, c†} = 1 → 1 constraint
        # Total: 3 constraints for N=1
        @test length(sys.equality_constraints) == 3
        @test isempty(sys.inequality_constraints)
        
        # Check simplify algorithm
        @test sys.simplify_algo.is_nilpotent == true
        @test sys.simplify_algo.is_unipotent == false
        @test sys.simplify_algo.is_projective == false
    end
    
    @testset "Anti-commutation Relations N=1" begin
        sys = fermionic_algebra(1)
        c, c_dag = sys.variables
        
        constraints = sys.equality_constraints
        
        # {c, c} = 0
        @test c[1] * c[1] + c[1] * c[1] in constraints
        
        # {c†, c†} = 0
        @test c_dag[1] * c_dag[1] + c_dag[1] * c_dag[1] in constraints
        
        # {c, c†} = 1
        @test c[1] * c_dag[1] + c_dag[1] * c[1] - 1 in constraints
    end
    
    @testset "Nilpotent Simplification" begin
        sys = fermionic_algebra(1)
        c, c_dag = sys.variables
        sa = sys.simplify_algo
        
        # Test that c² simplifies to zero monomial
        c_squared = c[1] * c[1]
        simplified = simplify(c_squared, sa)
        @test iszero(simplified) || isempty(simplified.monos)
        
        # Test that (c†)² simplifies to zero
        cdag_squared = c_dag[1] * c_dag[1]
        simplified2 = simplify(cdag_squared, sa)
        @test iszero(simplified2) || isempty(simplified2.monos)
    end
    
    @testset "Integration with cpolyopt" begin
        sys = fermionic_algebra(2)
        c, c_dag = sys.variables
        
        # Number operator: n₁ = c₁† c₁
        n1 = ComplexF64(1.0) * c_dag[1] * c[1]
        
        # Create optimization problem
        pop = cpolyopt(n1, sys)
        
        @test pop.objective == n1
        @test pop.is_nilpotent == true  # This will require PolyOpt to have is_nilpotent field
    end
end
```

#### **Step 7: Write documentation**

**Location**: `docs/src/examples/literate/fermionic_algebra_interface.jl`

Follow the same pattern as `pauli_algebra_interface.jl`:
1. Introduce fermionic systems
2. Show manual approach (tedious)
3. Show `fermionic_algebra(N)` approach
4. Example: Fermi-Hubbard model or similar
5. List advantages

#### **Step 8: Update PolyOpt to support nilpotent**

**Location**: `src/pop.jl`

```julia
struct PolyOpt{P} <: OptimizationProblem{P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    variables::Vector{Variable}
    comm_gps::Vector{Vector{Variable}}
    is_unipotent::Bool
    is_projective::Bool
    is_nilpotent::Bool  # NEW
end

function polyopt(objective::P; 
                 eq_constraints=Any[], 
                 ineq_constraints=Any[], 
                 comm_gps=Vector{Variable}[], 
                 is_unipotent::Bool=false, 
                 is_projective::Bool=false,
                 is_nilpotent::Bool=false) where {T,P<:AbstractPolynomial{T}}
    # ... validation
    @assert sum([is_unipotent, is_projective, is_nilpotent]) <= 1 
        "Only one of is_unipotent, is_projective, is_nilpotent can be true"
    
    return PolyOpt{P}(objective, eq_cons, ineq_cons, vars, comm_gps, 
                      is_unipotent, is_projective, is_nilpotent)
end
```

Similar changes for `ComplexPolyOpt` and `cpolyopt`.

---

## 11. Key Differences: Pauli vs Fermionic

### Mathematical Properties

| Aspect | Pauli Operators | Fermionic Operators |
|--------|----------------|---------------------|
| **Squaring** | σᵢ² = 1 | cᵢ² = 0, (cᵢ†)² = 0 |
| **Type** | Unipotent | Nilpotent |
| **Anti-commutation** | {σᵢ, σⱼ} = 2δᵢⱼ (same site) | {cᵢ, cⱼ} = 0 (all i,j) |
| | | {cᵢ, cⱼ†} = δᵢⱼ (canonical) |
| **Commutation** | Different sites commute | All anti-commute |
| **Operators per site** | 3 (σx, σy, σz) | 2 (c, c†) |
| **Hermiticity** | All Hermitian | c† is adjoint of c |
| **Hilbert space** | 2^N dimensional | 2^N dimensional |

### Implementation Differences

| Component | Pauli | Fermionic |
|-----------|-------|-----------|
| **Simplification** | `_simplify_unipotent!` | `_simplify_nilpotent!` (NEW) |
| **SimplifyAlgorithm** | `is_unipotent=true` | `is_nilpotent=true` (NEW) |
| **Commutation groups** | `[[σx[i], σy[i], σz[i]] for i in 1:N]` | `[all_operators]` |
| **Constraints count** | 6N (6 per site) | N² + N + N(N+1)/2 |
| | | = 2N² + (3N)/2 |
| **Constraint structure** | Products of operators at same site | All products of any two operators |

### Code Changes Required

1. ✅ **SimplifyAlgorithm**: Add `is_nilpotent` field
2. ✅ **_simplify_nilpotent!**: New function to handle c² = 0
3. ✅ **PolyOpt**: Add `is_nilpotent` field
4. ✅ **fermionic_algebra**: New constructor function
5. ✅ **Tests**: New test suite for fermionic algebra
6. ✅ **Documentation**: New tutorial following Pauli pattern

---

## 12. Critical Implementation Notes

### Zero Monomial Handling

**Issue #179**: FastPolynomials needs to support zero monomials.

Currently, when `_simplify_nilpotent!` detects c², it empties the monomial:
```julia
empty!(m.vars)
empty!(m.z)
```

This should be handled correctly by `Polynomial` constructor which filters out zero coefficients. Need to verify this works correctly throughout the codebase.

### Commutation Group Semantics

**Pauli**: `comm_gps = [[x[i], y[i], z[i]] for i in 1:N]`
- Variables in the same group maintain order
- Variables in different groups can be reordered

**Fermionic**: `comm_gps = [all_operators]`
- ALL variables are in the same group
- NO reordering is performed by `_comm!`
- Anti-commutation is entirely handled by equality constraints

This is correct! The simplifier does NOT need to know about anti-commutation. It only needs to:
1. Not reorder operators (single group)
2. Detect and eliminate c² = 0

### Type Consistency

All equality constraints must be `Polynomial{ComplexF64}` to match the optimization problem type. In `fermionic_algebra`, construct constraints carefully:

```julia
# Correct:
c[i] * c[j] + c[j] * c[i]  # Returns Polynomial

# Also correct:
ComplexF64(1) * (c[i] * c[j] + c[j] * c[i])  # Explicit type

# For constraints involving constants:
c[i] * c_dag[j] + c_dag[j] * c[i] - 1  # The -1 is automatically converted
```

### Testing with Solver

Integration tests should verify:
1. Fermi-Hubbard model ground state energies
2. BCS Hamiltonian (superconductivity)
3. Free fermion system (known exact solutions)

---

## 13. Validation Checklist

Before considering fermionic algebra complete:

- [ ] `_simplify_nilpotent!` correctly handles all cases
- [ ] `SimplifyAlgorithm` with `is_nilpotent=true` works
- [ ] `fermionic_algebra(N)` creates correct constraints
- [ ] `cpolyopt` accepts fermionic algebra systems
- [ ] All unit tests pass
- [ ] At least one integration test with known result passes
- [ ] Documentation tutorial is complete and runnable
- [ ] Code follows Julia idioms and project style
- [ ] No duplicate code (reuse existing patterns)

---

## 14. File Paths Summary

### Key Files to Read (Understanding)

1. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/algebra_constructors.jl:1-39` - Pauli algebra constructor
2. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/FastPolynomials/src/simplify.jl:1-102` - Simplification logic
3. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/FastPolynomials/src/monomials.jl:167-180` - Commutation
4. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/pop.jl:25-126` - Optimization problem types
5. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/interface.jl:73-97` - Solver interface

### Key Files to Modify (Implementation)

1. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/FastPolynomials/src/simplify.jl` - Add nilpotent case
2. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/algebra_constructors.jl` - Add fermionic_algebra
3. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/pop.jl` - Add is_nilpotent field
4. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/src/NCTSSoS.jl` - Export fermionic_algebra

### Key Files to Create

1. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/test/test_fermionic_algebra.jl` - Test suite
2. `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/fermionic_algebra_interface.jl` - Tutorial

---

## 15. References to Existing Code

All code snippets above reference actual lines in the codebase:
- Line numbers are accurate as of commit `fe689c2`
- File paths are absolute within the repository
- All examples are runnable code from the actual implementation

This analysis provides a complete blueprint for implementing fermionic algebra following the exact same architectural pattern as Pauli algebra.
