# Test Coverage Matrix

Maps source files to their corresponding test files.

## Core Types (`src/types/`)

| Source File | Test File(s) | Notes |
|-------------|--------------|-------|
| `algebra.jl` | `polynomials/algebra_types.jl` | 6 algebra singletons |
| `monomial.jl` | `polynomials/monomials.jl` | Construction, comparison, iteration |
| `term.jl` | `polynomials/term.jl` | Coefficient handling |
| `polynomial.jl` | `polynomials/polynomial.jl`, `polynomials/arithmetic.jl` | Full arithmetic suite |
| `registry.jl` | `polynomials/variables.jl` | Variable creation |
| `composed.jl` | `polynomials/composed_monomial.jl` | State/trace monomials |

## Simplification (`src/simplification/`)

| Source File | Test File(s) | Notes |
|-------------|--------------|-------|
| `noncommutative.jl` | `polynomials/simplify.jl` | No simplification rules |
| `pauli.jl` | `polynomials/simplify.jl` | σ² = I, cyclic products |
| `unipotent.jl` | `polynomials/simplify.jl` | U² = I |
| `projector.jl` | `polynomials/simplify.jl` | P² = P |
| `fermionic.jl` | `polynomials/simplify.jl`, `problems/fermionic/` | Anticommutation |
| `bosonic.jl` | `polynomials/simplify.jl`, `problems/condensed_matter/bose_hubbard.jl` | Commutation |

## Algorithms (`src/algorithms/`)

| Source File | Test File(s) | Notes |
|-------------|--------------|-------|
| `basis.jl` | `polynomials/basis.jl`, `polynomials/state_basis.jl` | Monomial basis generation |
| `canonicalization.jl` | `polynomials/canonicalization.jl` | Symmetric ordering |

## States (`src/states/`)

| Source File | Test File(s) | Notes |
|-------------|--------------|-------|
| `types.jl` | `polynomials/statepolynomial.jl` | NCStateWord, NCStatePolynomial types |
| `word.jl` | `polynomials/state_word.jl` | State word operations |
| `polynomial.jl` | `polynomials/statepolynomial.jl` | State polynomial arithmetic |

## Optimization (`src/optimization/`)

| Source File | Test File(s) | Notes |
|-------------|--------------|-------|
| `problem.jl` | `relaxations/interface.jl` | PolyOpt, StatePolyOpt construction |
| `interface.jl` | `relaxations/interface.jl` | cs_nctssos, SolverConfig |
| `sparsity.jl` | `relaxations/sparsity.jl`, `problems/trace_polynomial/trace_polynomial.jl` | CS/TS algorithms |
| `elimination.jl` | `relaxations/sparsity.jl` | MF, MMD, NoElimination |
| `sos.jl` | `relaxations/sos.jl` | SOS dualization, Cαj |
| `moment.jl` | `relaxations/interface.jl` | Moment formulation |
| `gns.jl` | `relaxations/gns.jl` | GNS reconstruction |

## Utility (`src/util/`)

| Source File | Test File(s) | Notes |
|-------------|--------------|-------|
| `helpers.jl` | `polynomials/utils.jl` | Helper functions |

---

## Problem Coverage

| Problem Domain | Test File(s) | Oracle File | Solver |
|----------------|--------------|-------------|--------|
| CHSH Bell | `problems/bell_inequalities/chsh.jl` | `chsh_oracles.jl` | Both |
| I3322 Bell | `problems/bell_inequalities/i3322.jl` | `i3322_oracles.jl` | Both |
| General Bell | `problems/bell_inequalities/bell_inequalities.jl` | - | Both |
| NC Examples | `problems/nc_polynomial/nc_examples.jl` | `example1_oracles.jl`, `corr_sparsity_oracles.jl` | Both |
| State Poly | `problems/state_polynomial/state_polynomial.jl` | `state_poly_oracles.jl` | Both |
| Trace Poly | `problems/trace_polynomial/trace_polynomial.jl` | `trace_poly_oracles.jl` | Both |
| Heisenberg Star | `problems/condensed_matter/heisenberg_star.jl` | `heisenberg_star_oracles.jl` | Mosek |
| Heisenberg Chain | `problems/condensed_matter/heisenberg_chain.jl` | - | Mosek |
| Heisenberg (XXX) | `problems/condensed_matter/heisenberg.jl` | - | Mosek |
| XY Model | `problems/condensed_matter/xy_model.jl` | - | Mosek |
| Ising | `problems/condensed_matter/ising.jl` | - | Mosek |
| PXP | `problems/condensed_matter/pxp.jl` | - | Mosek |
| Bose-Hubbard | `problems/condensed_matter/bose_hubbard.jl` | - | Mosek |
| Fermionic | `problems/fermionic/fermionic.jl` | - | Mosek |
| Fermionic Chain | `problems/fermionic/fermionic_chain.jl` | - | Mosek |
| Bilocal Networks | `problems/quantum_networks/bilocal_networks.jl` | - | Mosek |
| Benchmarks | `problems/benchmarks/ncpop_benchmarks.jl` | - | Mosek |

---

## Test Categories

| Category | Command | Solver | Description |
|----------|---------|--------|-------------|
| Polynomials | `make test-polynomials` | None | Core algebra tests |
| Quality | `make test-quality` | None | Aqua, ExplicitImports, Doctest |
| Relaxations | `make test-solvers` | COSMO/Mosek | SOS, sparsity, GNS |
| Problems | `make test-physics` | Mosek | Physics models |
| Full Suite | `make test` | Mosek | All tests |
| CI Suite | `make test-ci` | COSMO | No physics |
| Single File | `make test-file FILE=...` | Mosek | Any test file |

---

## Running Tests

```bash
# Full suite (requires Mosek)
make test

# CI suite (COSMO, no physics)
make test-ci

# Individual categories
make test-polynomials
make test-quality
make test-solvers
make test-physics

# Single test file
make test-file FILE=test/relaxations/sparsity.jl
make test-file FILE=test/problems/bell_inequalities/chsh.jl
```
