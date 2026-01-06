# NCTSSoS.jl Test Suite Review Plan

## Objective

Systematic walkthrough of all test cases to understand coverage, patterns, and gaps.

---

## Phase 1: Polynomials (Core Algebra)

**Goal**: Understand type system and simplification rules
**Files**: 16 | **Lines**: ~7,200 | **Solver**: None

### 1.1 Type System Foundation
- [ ] `algebra_types.jl` (528 lines) - 6 algebra singletons, dispatch behavior
- [ ] `variables.jl` (347 lines) - Registry creation, `create_*_variables` API

### 1.2 Monomial & Term Layer
- [ ] `monomials.jl` (564 lines) - Word representation, indexing
- [ ] `term.jl` (346 lines) - Coefficient × monomial operations
- [ ] `composed_monomial.jl` (500 lines) - Multi-variable composition

### 1.3 Polynomial Operations
- [ ] `polynomial.jl` (862 lines) - Core polynomial structure
- [ ] `arithmetic.jl` (137 lines) - Add/multiply/scalar ops
- [ ] `compare.jl` (115 lines) - Equality, ordering

### 1.4 Simplification (Critical)
- [ ] `canonicalization.jl` (524 lines) - Normal ordering
- [ ] `simplify.jl` (1,776 lines) - **Largest file**: algebra-specific rules
  - NonCommutative (no simplification)
  - Pauli (σ² = I, cyclic)
  - Fermionic (anticommutation)
  - Bosonic (commutation)
  - Projector (P² = P)
  - Unipotent (U² = I)

### 1.5 Basis Construction
- [ ] `basis.jl` (237 lines) - `get_ncbasis` generation

### 1.6 State Polynomials
- [ ] `state_word.jl` (358 lines) - NCStateWord representation
- [ ] `statepolynomial.jl` (336 lines) - NCStatePolynomial type
- [ ] `state_basis.jl` (171 lines) - State basis construction

### 1.7 Utilities & Quality
- [ ] `utils.jl` (106 lines) - Degree, extraction helpers
- [ ] `matrix_oracles.jl` (151 lines) - Oracle matrices
- [ ] `allocations.jl` (97 lines) - Memory efficiency

---

## Phase 2: Relaxations (Algorithm Components)

**Goal**: Understand SOS dualization and sparsity decomposition
**Files**: 4 | **Lines**: ~936 | **Solver**: Required

### 2.1 Problem Interface
- [ ] `interface.jl` (283 lines) - PolyOpt/StatePolyOpt constructors, constraints

### 2.2 Sparsity Decomposition
- [ ] `sparsity.jl` (438 lines) - **Key file**
  - Correlative sparsity graph construction
  - Clique decomposition (MF, MMD)
  - Term sparsity reduction
  - Constraint assignment to cliques

### 2.3 SOS Components
- [ ] `sos.jl` (54 lines) - Cαj coefficient extraction

### 2.4 GNS Reconstruction
- [ ] `gns.jl` (140 lines) - Moment matrix reconstruction (mostly disabled)

---

## Phase 3: Problems (Integration Tests)

**Goal**: Verify optimization correctness against NCTSSOS oracles
**Files**: 17 | **Lines**: ~3,000 | **Solver**: Required (some Mosek-only)

### 3.1 Bell Inequalities (Quantum Nonlocality)
- [ ] `chsh.jl` (222 lines) - CHSH bound (-2√2), 3 formulations
- [ ] `i3322.jl` (113 lines) - I_3322 inequality
- [ ] `bell_inequalities.jl` (174 lines) - General framework

### 3.2 NC Polynomial Examples
- [ ] `nc_examples.jl` (346 lines) - NCTSSOS paper examples
  - Example 1: unconstrained 3-var
  - Example 2: constrained 2-var
  - Large-scale CS+TS (n=10)

### 3.3 State Polynomial (7.2.x)
- [ ] `state_polynomial.jl` (155 lines) - ς() operator tests

### 3.4 Trace Polynomial (6.x)
- [ ] `trace_polynomial.jl` (161 lines) - tr() operator tests
- [ ] `sparsity_variants.jl` (44 lines) - Algorithm comparison

### 3.5 Benchmarks
- [ ] `ncpop_benchmarks.jl` (31 lines) - Classical optimization (Rosenbrock)

### 3.6 Quantum Networks (--local only)
- [ ] `bilocal_networks.jl` (161 lines) - Alice-Bob-Charlie scenarios

### 3.7 Condensed Matter (DISABLED)
- [ ] `heisenberg_star.jl` (192 lines)
- [ ] `heisenberg_chain.jl` (28 lines)
- [ ] `pxp.jl` (146 lines)
- [ ] `xy_model.jl` (37 lines)
- [ ] `ising.jl` (31 lines)

### 3.8 Fermionic (DISABLED)
- [ ] `fermionic.jl` (210 lines)
- [ ] `fermionic_chain.jl` (60 lines)

---

## Phase 4: Quality Checks

**Goal**: Verify code hygiene
**Files**: 3 | **Lines**: ~53 | **Solver**: None

- [ ] `Aqua.jl` - Ambiguities, type stability, piracy
- [ ] `ExplicitImports.jl` - Import hygiene
- [ ] `Doctest.jl` - Docstring examples

---

## Phase 5: Oracles (Reference Data)

**Goal**: Understand validation strategy
**Files**: 3 data + 14 scripts

### Oracle Data Files
- [ ] `example2_oracles.jl` (470 lines)
- [ ] `state_poly_oracles.jl` (1,168 lines)
- [ ] `heisenberg_star_oracles.jl` (929 lines)

### Oracle Generation (run in NCTSSOS repo)
- [ ] Review `oracle_utils.jl` - Shared utilities
- [ ] Review generation scripts pattern

---

## Review Checklist Per File

For each test file, evaluate:

1. **Coverage**: What functionality is tested?
2. **Edge cases**: Are boundary conditions covered?
3. **Algebra types**: Which algebras are tested?
4. **Sparsity variants**: Dense/CS/TS/CS+TS?
5. **Oracle validation**: Reference values present?
6. **Gaps**: What's missing or disabled?

---

## Key Questions to Answer

### Architecture
- [ ] How do the 6 algebra types differ in simplification?
- [ ] What's the relationship between NC, State, and Trace polynomials?
- [ ] How does sparsity decomposition flow through the solver?

### Coverage Gaps
- [ ] Why is GNS mostly disabled?
- [ ] Why are condensed matter tests commented out?
- [ ] What's blocking fermionic tests?

### Validation
- [ ] Are oracle values sufficient for regression detection?
- [ ] What tolerance is used for floating-point comparison?
- [ ] How are moment matrix sizes validated?

---

## Execution Commands

```bash
# Run specific phases
make test-polynomials    # Phase 1 only
make test-relaxations    # Phase 2 only
make test-problems       # Phase 3 only
make test-quality        # Phase 4 only

# Full suite
make test                # With Mosek (--local)
make test-ci             # CI mode (COSMO)

# Single file exploration
julia --project -e 'using Pkg; Pkg.test(test_args=["--polynomials"])'
```

---

## Notes

- Start with Phase 1 (no solver dependency)
- Phase 3 requires Mosek for full coverage
- Disabled tests (condensed matter, fermionic) need investigation
- Oracle generation requires separate NCTSSOS installation
