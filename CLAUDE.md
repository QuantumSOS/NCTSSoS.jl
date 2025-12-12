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

# Run only FastPolynomials tests
make test-FastPoly

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

### Two-Layer Design

1. **FastPolynomials** (`src/FastPolynomials/`): High-performance polynomial library
   - Type-parameterized algebra system with compile-time dispatch
   - Six algebra types: `NonCommutativeAlgebra`, `PauliAlgebra`, `FermionicAlgebra`, `BosonicAlgebra`, `ProjectorAlgebra`, `UnipotentAlgebra`
   - Core types: `Monomial{A,T}` → `Term{M,C}` → `Polynomial{A,T,C}`
   - `VariableRegistry` maps symbols ↔ integer indices with algebra-specific encoding

2. **NCTSSoS** (root `src/`): Optimization framework
   - `PolyOpt`/`ComplexPolyOpt`: Problem definitions with constraints
   - `CorrelativeSparsity`: Clique decomposition for block-structured SDP
   - `TermSparsity`: Further sparsity exploitation within cliques
   - `cs_nctssos()`: Main solver entry point

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

### Simplification System

Each algebra type has its own simplification rules in `src/FastPolynomials/src/simplification/`:
- Pauli: `σᵢ² = I`, anticommutation, cyclic products → complex phases
- Fermionic: `{aᵢ, aⱼ†} = δᵢⱼ` anticommutation → normal ordering
- Bosonic: `[cᵢ, cⱼ†] = δᵢⱼ` commutation → normal ordering with corrections
- Projector: `Pᵢ² = Pᵢ` idempotency
- Unipotent: `U² = I` involutory

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

- `test/fastpoly_test/`: FastPolynomials unit tests
- `test/*.jl`: Integration tests (moment solver, SOS solver, interface)
- Tests under `LOCAL_TESTING` env var take longer (moment problems)
