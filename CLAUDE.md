# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NCTSSoS.jl is a Julia package for solving sparse noncommutative polynomial optimization problems using the structured moment-SOHS (Sum of Hermitian Squares) hierarchy. Successor to Python NCTSSOS.

## Build & Test Commands

```bash
make test-minimal      # Fast smoke test (~30s) - run first
make test-polynomials  # Algebra tests only, no solver
make test              # Full suite with Mosek
make test-ci           # CI suite (COSMO solver, no physics)
make test-file FILE=test/relaxations/sparsity.jl  # Single file
make servedocs         # Live documentation server
```

Test flags via Pkg.test:
- `--minimal` / `--polynomials` / `--quality` / `--solvers` / `--physics`
- `--local` enables Mosek solver (required for physics tests)

## Architecture

### Type Hierarchy
```
AlgebraType (singleton dispatch)
├── NonCommutativeAlgebra   # Generic, no simplification rules
├── PauliAlgebra            # σᵢ² = I, anticommutation
├── FermionicAlgebra        # {aᵢ, aⱼ†} = δᵢⱼ
├── BosonicAlgebra          # [cᵢ, cⱼ†] = δᵢⱼ
├── ProjectorAlgebra        # Pᵢ² = Pᵢ
└── UnipotentAlgebra        # U² = I

Monomial{A<:AlgebraType, T<:Integer} → Term → Polynomial
```

### Solver Pipeline
1. `polyopt(obj; eq_constraints, ineq_constraints)` → PolyOpt
2. Sparsity detection → clique decomposition (CliqueTrees)
3. Moment/SOS matrix construction
4. SDP solve via JuMP (Clarabel, COSMO, Mosek)
5. `cs_nctssos(prob; config)` → PolyOptResult

### Key Modules
- `src/types/` - Core algebraic types (algebra.jl → monomial.jl → polynomial.jl)
- `src/simplification/` - Algebra-specific reduction rules (dispatch on AlgebraType)
- `src/algorithms/` - Canonicalization (symmetric_canon, cyclic_canon) and basis generation
- `src/optimization/` - SDP relaxation: sparsity.jl, moment.jl, sos.jl, interface.jl
- `src/states/` - StatePolynomial for quantum expectations with trace operations

## Code Conventions

### Fermionic Operators
Signed integers encode creation vs annihilation: positive = creation (a†), negative = annihilation (a). Example: `[1, -1, 2]` represents a₁†a₁a₂†.

### Pauli Variables
Organized by site with special commutation handling. Use `create_pauli_variables(n)` to get site-indexed operators.

### Simplification
All simplification goes through `simplify()` which dispatches on AlgebraType. Never manually manipulate monomial words - use the simplification pipeline.

### Canonicalization
- `symmetric_canon(m)` - word + reverse for Hermitian equivalence
- `cyclic_canon(m)` - trace equivalence (cyclic permutations)

## Public API Entry Points

```julia
# Variable creation
create_pauli_variables(n), create_fermionic_variables(n), create_bosonic_variables(n)
create_projector_variables(n), create_unipotent_variables(n), create_noncommutative_variables(n)

# Problem setup
prob = polyopt(objective; eq_constraints=[...], ineq_constraints=[...])
result = cs_nctssos(prob; config=SolverConfig(...))

# Polynomial operations
degree(p), monomials(p), coefficients(p), terms(p), variables(p)
simplify(p), simplify!(p), canonicalize(p)
```

## Testing

Test structure:
- `test/polynomials/` - Algebra tests (no solver needed)
- `test/relaxations/` - SOS/moment component tests
- `test/problems/` - Integration tests with solvers
- `test/quality/` - Aqua + ExplicitImports checks
- `test/oracles/` - Reference values from NCTSSOS Python
