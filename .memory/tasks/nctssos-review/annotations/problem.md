# Review: src/optimization/problem.jl

**Lines**: 310
**Purpose**: Define optimization problem types and constructors

## Structure

```
OptimizationProblem{A<:AlgebraType, P}  (abstract)
├── PolyOpt{A,T,P}         # Standard polynomial optimization
└── StatePolyOpt{A,T,ST,P}  # State polynomial (quantum expectations)
```

## Key Components

### 1. `PolyOpt` (lines 54-59)
- Fields: `objective`, `eq_constraints`, `ineq_constraints`, `registry`
- Type params carry algebra type `A` → enables dispatch on algebraic structure
- Clean design: no manual flags (`comm_gps`, `is_unipotent`), algebra type encodes behavior

### 2. `polyopt()` constructor (lines 102-118)
- Integer coefficient check (JuMP requirement)
- Constraint deduplication via `unique!(copy(...))`
- Type inference from polynomial and registry

### 3. `_is_complex_problem()` trait (lines 168-191)
- Maps algebra type → real vs complex moment relaxation
- Complex: Pauli, Fermionic, Bosonic (phases from products)
- Real: NonCommutative, Projector, Unipotent (no phases)

### 4. `StatePolyOpt` (lines 237-242)
- Extends to NCStatePolynomial for quantum state expectations
- Used for Bell inequality optimization (CHSH example in docs)

## Quality Assessment

| Aspect | Rating | Notes |
|--------|--------|-------|
| Documentation | ✓✓✓ | Excellent docstrings with examples |
| Type safety | ✓✓✓ | Parametric types enforce consistency |
| API design | ✓✓ | Clean, but two `polyopt` methods might confuse |
| Error handling | ✓✓ | Integer coefficient check, but no other validation |

## Observations

**Strengths**:
- Algebra type as first-class parameter eliminates flag sprawl
- Registry coupling ensures symbol resolution
- Trait-based dispatch (`_is_complex_problem`) is idiomatic Julia

**Minor concerns**:
1. **No constraint validation** - constraints aren't checked for matching algebra/registry
2. **No degree checks** - could silently accept objectives incompatible with solver
3. **Duplicate `polyopt` methods** - one for Polynomial, one for NCStatePolynomial.
   Single function with Union type might be cleaner.

## Code Quality: Good

No correctness issues found. Well-structured with clear separation of concerns.
