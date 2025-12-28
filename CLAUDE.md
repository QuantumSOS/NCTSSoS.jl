# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NCTSSoS.jl is a Julia package for solving sparse noncommutative polynomial optimization problems using structured moment-SOHS (Sum of Hermitian Squares) hierarchy. It's the successor to NCTSSOS and is designed for quantum operator algebra optimization.

## Build & Test Commands

```bash
# Initialize project (precompile dependencies)
make init

# Run all tests (includes extended tests)
LOCAL_TESTING=true make test

# Run tests without extended tests (faster, CI mode)
make test

# Run only polynomial algebra tests
julia --project -e 'using Pkg; Pkg.test()' -- test/polynomials/runtests.jl

# Run a single test file
julia --project -e 'using Pkg; Pkg.test()' -- test/pop.jl

# Generate documentation examples
make examples

# Serve docs locally with live reload
make servedocs

# Run benchmarks comparing to main branch
make bench TARGET=main
```

## Architecture

### Source Structure

```
src/
├── NCTSSoS.jl              # Module entry point
├── types/                  # Core data structures
│   ├── algebra.jl          # AlgebraType hierarchy (6 concrete types)
│   ├── registry.jl         # VariableRegistry: symbols ↔ indices
│   ├── monomial.jl         # Monomial{A,T}
│   ├── term.jl             # Term{M,C}
│   ├── polynomial.jl       # Polynomial{A,T,C}
│   └── composed.jl         # ComposedMonomial for tensor products
├── simplification/         # Algebra-specific simplification rules
│   ├── pauli.jl            # σᵢ² = I, anticommutation
│   ├── fermionic.jl        # {aᵢ, aⱼ†} = δᵢⱼ, normal ordering
│   ├── bosonic.jl          # [cᵢ, cⱼ†] = δᵢⱼ, normal ordering
│   ├── projector.jl        # Pᵢ² = Pᵢ idempotency
│   ├── unipotent.jl        # U² = I involutory
│   └── noncommutative.jl   # No simplification (generic)
├── algorithms/             # Pure algorithms
│   ├── canonicalization.jl # symmetric_canon, cyclic_canon
│   └── basis.jl            # get_ncbasis, get_state_basis
├── states/                 # State polynomial layer
│   ├── types.jl            # StateWord, StatePolynomial types
│   ├── word.jl             # StateWord algebra
│   └── polynomial.jl       # StatePolynomial operations
├── optimization/           # SDP framework
│   ├── problem.jl          # PolyOpt definition
│   ├── sparsity.jl         # CorrelativeSparsity, TermSparsity
│   ├── moment.jl           # Moment relaxation
│   ├── sos.jl              # SOS relaxation
│   ├── gns.jl              # GNS construction
│   ├── elimination.jl      # Elimination strategies
│   └── interface.jl        # Public solver API (cs_nctssos)
└── util/
    └── helpers.jl          # Utility functions
```

### Type Hierarchy

- **AlgebraType**: `NonCommutativeAlgebra`, `PauliAlgebra`, `FermionicAlgebra`, `BosonicAlgebra`, `ProjectorAlgebra`, `UnipotentAlgebra`
- **Core Types**: `Monomial{A,T}` → `Term{M,C}` → `Polynomial{A,T,C}`
- **Optimization**: `PolyOpt` → `CorrelativeSparsity` → `TermSparsity` → SDP

### Key Data Flow

```
Variables → Monomials → simplify() → Terms → Polynomial
                                              ↓
PolyOpt → correlative_sparsity() → CorrelativeSparsity
                                              ↓
         term_sparsities() → TermSparsity blocks
                                              ↓
         moment_relax() → SDP → JuMP model → solve
```

## Code Conventions

- Monomials use integer indices internally; `VariableRegistry` handles symbol mapping
- Polynomial invariants: sorted terms, no duplicates, no zero coefficients
- Pauli algebra uses `ComplexF64` coefficients; others use `Float64`
- `@ncpolyvar` macro creates legacy-compatible variables
- Use `create_*_variables()` functions for type-safe variable creation

## Dependencies

- **JuMP**: Optimization modeling
- **CliqueTrees/ChordalGraph**: Sparsity graph algorithms
- **Clarabel/COSMO/MosekTools**: SDP solvers (test dependencies)

## Test Structure

- `test/polynomials/`: Polynomial algebra unit tests (types, simplification, arithmetic)
- `test/*.jl`: Integration tests (moment solver, SOS solver, interface)
- Physics model tests: `heisenberg.jl`, `xy_model.jl`, `bose_hubbard.jl` (run with `LOCAL_TESTING=true`)
- Quality checks: `Aqua.jl`, `ExplicitImports.jl`
