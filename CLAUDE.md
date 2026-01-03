# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NCTSSoS.jl solves sparse noncommutative polynomial optimization via structured moment-SOHS hierarchy. Successor to NCTSSOS.

## Build & Test Commands

```bash
make init                # Install deps, precompile
make test                # Full suite with Mosek (--local)
make test-ci             # CI suite with COSMO
make test-polynomials    # Core algebra only (no solver)
make test-quality        # Aqua, ExplicitImports, Doctest

# Direct Pkg.test
julia --project -e 'using Pkg; Pkg.test(test_args=["--polynomials"])'
julia --project -e 'using Pkg; Pkg.test(test_args=["--local"])'
julia --project -e 'using Pkg; Pkg.test(test_args=["--relaxations", "--problems"])'

make servedocs           # Serve docs locally
```

## Generating Oracle Test Results

Oracles = reference results from original NCTSSOS package for validation.

**Run in separate NCTSSOS repo** (requires MosekTools):
```bash
cd ~/NCTSSOS
julia --project /path/to/NCTSSoS.jl/test/oracles/scripts/nctssos_chsh.jl
```

**Structure:**
- `test/oracles/scripts/nctssos_*.jl` - Scripts running NCTSSOS
- `test/oracles/results/*_oracles.jl` - Generated oracle dicts
- `test/oracles/scripts/oracle_utils.jl` - Shared utilities

**Oracle format:** `(opt, sides, nuniq)` where:
- `opt` = optimal objective value
- `sides` = moment matrix block sizes
- `nuniq` = unique moment indices (affine constraints)

## Architecture

### Algebra Types (`src/types/algebra.jl`)

Singleton types for dispatch:
- `NonCommutativeAlgebra` - generic NC (no simplification)
- `PauliAlgebra` - σ² = I, cyclic products
- `FermionicAlgebra` - anticommutation
- `BosonicAlgebra` - commutation
- `ProjectorAlgebra` - P² = P
- `UnipotentAlgebra` - U² = I

### Polynomial Stack

1. `Monomial{A,T}` - word of variable indices
2. `Term{M,C}` - coefficient × monomial
3. `Polynomial{A,T,C}` - collection of terms

State polynomials in `src/states/`: `NCStateWord`, `NCStatePolynomial`

### Variable Creation

```julia
reg, (σx, σy, σz) = create_pauli_variables(1:n)
reg, (a,) = create_fermionic_variables(1:n)
reg, (x,) = create_noncommutative_variables(:x, 1:n)
```

### Optimization Pipeline

1. `polyopt(objective, registry; eq_constraints=[], ineq_constraints=[])` → `PolyOpt`
2. `SolverConfig(optimizer=Mosek.Optimizer; order=3, cs_algo=MF(), ts_algo=MMD())`
3. `cs_nctssos(pop, solver_config)` - main solver
4. `cs_nctssos_higher(pop, prev_result, solver_config)` - higher-order refinement

### Sparsity

- **Correlative**: decomposes into cliques (`cs_algo`)
- **Term**: reduces block sizes (`ts_algo`)
- Options: `MF()`, `MMD()`, `NoElimination()`

## Test Structure

```
test/
├── polynomials/     # Core algebra (no solver)
├── relaxations/     # SOS, sparsity, GNS
├── problems/        # Physics: bell_inequalities, condensed_matter, etc.
├── oracles/         # NCTSSOS reference values
└── setup.jl         # COSMO default, Mosek with --local
```

## Solver Notes

- **CI default**: COSMO (open-source)
- **Full suite**: Mosek with `--local` (required for physics tests)
