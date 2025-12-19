# Research Report: State Polynomial Implementation in NCTSSoS.jl

## Executive Summary

NCTSSoS.jl has a **substantial and largely complete** implementation for state polynomial optimization. The core type system, arithmetic operations, and solver pipeline integration are all functional. However, there are some important **known limitations** and **missing features** that affect certain use cases.

### Key Findings

- **Working**: StatePolynomial, NCStatePolynomial types, solver pipeline (cs_nctssos), basic tests passing
- **Known Limitations**: Compound StateWords (squared expectations like `<A>²`), term sparsity (MMD) for state polynomials
- **Missing**: Direct moment solving for StatePolyOpt (dualize=false not implemented)

---

## 1. Current State of Implementation

### 1.1 FastPolynomials Layer

#### State Types (`state_types.jl`)
**Status: COMPLETE**

Two state types are implemented:
```julia
abstract type StateType end
struct Arbitrary <: StateType end     # General quantum state: <M>
struct MaxEntangled <: StateType end  # Trace state: tr(M)
```

#### StateWord (`state_word.jl`)
**Status: COMPLETE** (850 lines)

`StateWord{ST,A,T}` represents products of expectations `<M1><M2>...<Mk>`:
- Stores sorted vector of canonicalized monomials
- Involution invariant: each monomial is canonicalized to `min(m, adjoint(m))`
- Complete arithmetic: multiplication, comparison, hashing
- Precomputed hash for fast equality checks

`NCStateWord{ST,A,T}` represents `<M1><M2>...<Mk> * Onc`:
- Combines commutative StateWord with non-commutative operator monomial
- Full arithmetic support
- `expval(ncsw)` converts to StateWord by adding nc_word as expectation

Convenience constructors:
- `ς(m::Monomial)` → `StateWord{Arbitrary}`
- `tr(m::Monomial)` → `StateWord{MaxEntangled}`

Basis generation:
- `get_state_basis(registry, d; state_type)` → `Vector{NCStateWord}`

#### StatePolynomial (`state_polynomial.jl`)
**Status: COMPLETE** (903 lines)

`StatePolynomial{C,ST,A,T}` - polynomial in StateWords:
- Automatic sorting, deduplication, zero-coefficient removal
- Full arithmetic: addition, subtraction, scalar multiplication, polynomial multiplication
- Accessors: `coefficients()`, `monomials()`, `terms()`, `degree()`, `variables()`

`NCStatePolynomial{C,ST,A,T}` - polynomial in NCStateWords:
- Same invariants and operations as StatePolynomial
- Additional: `variable_indices()`, `maxdegree()` for solver integration
- `expval(ncsp)` converts to StatePolynomial

#### Utility Functions (`utils.jl`)
**Status: COMPLETE**

Key functions for optimization:
- `_neat_dot3(a::NCStateWord, m, b::NCStateWord)` → `NCStatePolynomial`
  - Handles middle argument as Monomial, NCStateWord, or StateWord
  - Applies algebra-specific simplification rules
- `neat_dot(a::Monomial, b::Monomial)` → `Monomial` (adjoint(a) * b)
- `ς(p::Polynomial)` → `StatePolynomial` (convert polynomial to state expectations)

### 1.2 NCTSSoS Optimization Layer

#### Problem Definition (`pop.jl`)
**Status: COMPLETE**

```julia
struct StatePolyOpt{A<:AlgebraType, ST<:StateType, P<:NCStatePolynomial} <: OptimizationProblem{A, P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    registry::VariableRegistry{A}
end
```

Factory function `polyopt(::NCStatePolynomial, registry; eq_constraints, ineq_constraints)` creates `StatePolyOpt`.

#### Correlative Sparsity (`sparse.jl`)
**Status: COMPLETE**

`StateCorrelativeSparsity{A,ST,T,P,M}` mirrors `CorrelativeSparsity` for state polynomials:
- `correlative_sparsity(::StatePolyOpt, order, elim_algo)` builds sparsity structure
- `get_state_correlative_graph()` extracts variable correlations
- `assign_state_constraint()` distributes constraints to cliques
- Uses `get_state_basis()` for NCStateWord bases

Term sparsity functions have NCStateWord overloads:
- `init_activated_supp(::NCStatePolynomial, ...)`
- `term_sparsities(::Vector{NCStateWord}, ...)`
- `get_term_sparsity_graph(::Vector{NCStateWord}, ...)`
- `iterate_term_sparse_supp(::NCStatePolynomial, ...)`

#### Moment Solver (`moment_solver.jl`)
**Status: MOSTLY COMPLETE**

`StateMomentProblem{A,ST,T,M,P}` - symbolic moment problem for state polynomials:
- `moment_relax(::StatePolyOpt, corr_sparsity, term_sparsities)` builds the problem
- `_build_state_constraint_matrix()` constructs polynomial-valued matrices

**Gap**: `solve_moment_problem()` for `StateMomentProblem` not implemented (throws error for `dualize=false`)

#### SOS Solver (`sos_solver.jl`)
**Status: COMPLETE**

`sos_dualize(::StateMomentProblem)` converts to SOS dual:
- `_sos_dualize_state()` handles state polynomial optimization
- Creates matrix variables for each constraint
- Matches coefficients using `expval()` to convert NCStateWord to StateWord

Helper: `_get_state_Cαj()` extracts coefficient matrices for state polynomial constraints.

#### Interface (`interface.jl`)
**Status: COMPLETE**

`StatePolyOptResult{T,A,ST,TI,P,M}` stores optimization results.

`cs_nctssos(::StatePolyOpt, solver_config; dualize=true)`:
- Computes correlative sparsity
- Builds partial objectives per clique
- Generates term sparsities
- Calls `moment_relax()` then `sos_dualize()`
- Returns `StatePolyOptResult`

---

## 2. Test Coverage and Known Issues

### 2.1 Passing Tests

| Test | Description | Expected Value |
|------|-------------|----------------|
| 7.2.0 | CHSH-like Bell inequality | -2.8284... |
| 7.2.2 | Covariance optimization (3x3) | -5.0 |
| Example 6.2.0 | Trace version of 7.2.0 | -2.8284... |
| Example 6.2.2 | Trace version of 7.2.2 | -5.0 |

### 2.2 Known Failing Tests

| Test | Issue | Root Cause |
|------|-------|------------|
| 7.2.1 | Returns ~0 instead of -4.0 | Squared expectations limitation |
| 7.2.0 Sparse | Returns 0 with MMD term sparsity | Term sparsity graph construction |
| Example 6.1 | Returns incorrect value | Compound StateWords from trace products |
| Example 6.2.1 | Returns incorrect value | Same as 7.2.1 |

### 2.3 Root Cause Analysis: Squared Expectations

**Problem**: Objectives with terms like `<xy>²` (squared expectations) are not handled correctly.

**Technical Details**:

1. Current basis (`get_state_basis`) generates NCStateWords of form `<I>*M` (identity StateWord, operator Monomial)

2. When computing moment matrix entries via `_neat_dot3`:
   - `_neat_dot3(<I>*a, I, <I>*b)` → `<I>*(a†b)` (single expectation)
   - Cannot generate compound StateWords like `<A><B>` (products of expectations)

3. The objective has terms `<xy><xy>` which require compound StateWords with 2+ state_monos in StateWord. These are NOT in the basis.

4. Result: Terms not in basis are ignored, leading to an incorrectly relaxed problem.

**Attempted Fix (Reverted)**:
Adding `<M>*I` form alongside `<I>*M` form caused double-counting because both have same `expval()` but produce different results in `_neat_dot3`.

### 2.4 Term Sparsity Issue

MMD term sparsity algorithm creates incorrect block decomposition for state polynomial optimization. The term sparsity graph construction (`get_term_sparsity_graph`) doesn't properly handle NCStateWord types, leading to blocks that don't capture the full problem structure.

**Workaround**: Use `ts_algo=NoElimination()` for state polynomial optimization.

---

## 3. List of Missing Features and Gaps

### 3.1 Critical Gaps (Breaking Tests)

1. **Compound StateWord Support** (High Priority)
   - Objectives with squared/multiplied expectations (`<A>²`, `<A><B>`) give wrong results
   - Requires fundamental redesign of basis generation or moment matrix structure
   - Affects: Tests 7.2.1, Example 6.1, Example 6.2.1

2. **Term Sparsity for State Polynomials** (Medium Priority)
   - MMD algorithm doesn't work correctly
   - `get_term_sparsity_graph` for NCStateWord needs fixing
   - Affects: Test 7.2.0 Sparse

### 3.2 Missing Features (Non-Breaking)

3. **Direct Moment Solving** (`dualize=false`)
   - `solve_moment_problem(::StateMomentProblem, optimizer)` not implemented
   - Currently throws: "Direct moment solving (dualize=false) is not yet implemented for StatePolyOpt"
   - Only SOS dualization path works

4. **Higher-Order Iteration**
   - `cs_nctssos_higher(::StatePolyOpt, prev_res, solver_config)` not implemented
   - Can't do iterative refinement for state polynomial optimization

5. **GNS Reconstruction for State Polynomials**
   - `reconstruct()` in `gns.jl` only works with regular Polynomials
   - No equivalent for extracting state representations from state polynomial solutions

### 3.3 Minor Gaps

6. **StatePolynomial ↔ NCStatePolynomial Conversion**
   - `StatePolynomial * Monomial` → `NCStatePolynomial` works
   - No direct `NCStatePolynomial → StatePolynomial` conversion (only via `expval` which loses nc_word info)

7. **Documentation**
   - No user-facing documentation for state polynomial optimization
   - Missing examples in `docs/src/examples/`

---

## 4. Prioritized Implementation Plan

### Phase 1: Fix Critical Test Failures (HIGH PRIORITY)

#### Task 1.1: Compound StateWord Support
**Complexity**: HIGH (architectural change)
**Impact**: Fixes tests 7.2.1, Example 6.1, Example 6.2.1

**Options**:
1. **StateWord-indexed moment matrix**: Change basis from NCStateWord to StateWord
   - Variables for `<xy>` separate from `<x>`, `<y>`
   - Different SDP formulation needed
   
2. **Objective-aware basis generation**: Analyze objective to determine required NCStateWord forms
   - Generate `<M>*I` forms only when compound StateWords are needed
   - More complex but preserves current architecture

3. **Dual basis approach**: Use both `<I>*M` and `<M>*I` with proper scaling
   - Avoids double-counting by normalizing coefficients
   
**Recommendation**: Start with Option 3 as it's least invasive, fall back to Option 1 if needed.

#### Task 1.2: Fix Term Sparsity for State Polynomials
**Complexity**: MEDIUM
**Impact**: Fixes Test 7.2.0 Sparse, enables sparsity exploitation

**Implementation**:
- Review `get_term_sparsity_graph` for NCStateWord
- Ensure edges are added correctly based on NCStateWord structure
- Test with simple cases before complex problems

### Phase 2: Complete Solver Pipeline (MEDIUM PRIORITY)

#### Task 2.1: Implement Direct Moment Solving
**Complexity**: MEDIUM
**Impact**: Enables alternative solving path, useful for debugging

**Implementation**:
- Extend `_solve_real_moment_problem` for `StateMomentProblem`
- Create `_substitute_state_poly` to handle NCStatePolynomial → JuMP expression conversion
- Uses existing `expval()` to map NCStateWords to moment variables

#### Task 2.2: Higher-Order Iteration
**Complexity**: LOW
**Impact**: Enables iterative refinement

**Implementation**:
- Add `cs_nctssos_higher(::StatePolyOpt, ...)` following existing pattern
- Reuse term sparsity from previous result

### Phase 3: GNS and Advanced Features (LOW PRIORITY)

#### Task 3.1: GNS Reconstruction for State Polynomials
**Complexity**: HIGH
**Impact**: Enables state extraction from solutions

**Requires**: Understanding of state polynomial moment matrices

#### Task 3.2: Documentation and Examples
**Complexity**: LOW
**Impact**: User accessibility

**Tasks**:
- Add `state_poly_opt.jl` to `docs/src/examples/literate/`
- Document `StatePolyOpt`, `ς()`, `tr()` usage
- Add API documentation for state polynomial types

---

## 5. Architectural Considerations

### 5.1 Type Hierarchy

```
AbstractMonomial
├── Monomial{A,T}
├── StateWord{ST,A,T}
└── NCStateWord{ST,A,T}

AbstractPolynomial{T} (internal alias)
├── Polynomial{A,T,C}
├── StatePolynomial{C,ST,A,T}
└── NCStatePolynomial{C,ST,A,T}

OptimizationProblem{A,P}
├── PolyOpt{A,P}
└── StatePolyOpt{A,ST,P}
```

### 5.2 Key Design Decisions

1. **NCStateWord structure**: Separates commutative (StateWord) and non-commutative (Monomial) parts. This is mathematically correct for state polynomial optimization but creates the compound StateWord limitation.

2. **Involution invariant**: StateWord canonicalizes each monomial to `min(m, adjoint(m))`. This ensures `<M>` and `<M†>` are treated as equivalent.

3. **Basis generation**: `get_state_basis` generates `<I>*M` form NCStateWords. This is natural for operator optimization but insufficient for expectation-product problems.

4. **SOS dualization**: Uses StateWord-indexed variables after `expval()` conversion. This is correct but loses information about the nc_word structure.

### 5.3 Challenges for Compound StateWord Support

The fundamental issue is mathematical, not just implementation:

1. **Current approach**: NCStateWord-indexed moment matrix
   - Entries: `M[i,j] = <basis_i† * op * basis_j>` = single expectation
   - Cannot represent products of expectations

2. **The `<A>²` problem**: 
   - `<xy>²` is the SQUARE of an expectation value
   - NOT the same as `<(xy)²>` = `<xyxy>`
   - Requires different variable indexing scheme

3. **Potential solution**: State polynomial moment matrix indexed by StateWords
   - Variable for each unique StateWord (including compounds)
   - Different constraint structure
   - Significant architectural change

---

## 6. Sources

### Implementation Files
- `src/FastPolynomials/src/state_types.jl` (66 lines)
- `src/FastPolynomials/src/state_word.jl` (850 lines)
- `src/FastPolynomials/src/state_polynomial.jl` (903 lines)
- `src/FastPolynomials/src/utils.jl` (266 lines)
- `src/pop.jl` (309 lines)
- `src/sparse.jl` (868 lines)
- `src/moment_solver.jl` (686 lines)
- `src/sos_solver.jl` (493 lines)
- `src/interface.jl` (289 lines)

### Test Files
- `test/state_poly_opt.jl` (79 lines)
- `test/state_moment_solver.jl` (224 lines)
- `test/trace_poly_opt.jl` (80 lines)
- `test/fastpoly_test/statepolynomial.jl` (329 lines)
- `test/fastpoly_test/state_word.jl` (248 lines)

### Documentation
- `.claude/tasks/statepolyopt-solver/context.md` - Known issues and attempted fixes
- `.claude/tasks/statepolyopt-solver/research.md` - Original implementation plan
