# Progress Log: fastpoly-integration

## Session: 2025-12-13 - Phase 2 Complete (Registry-Based Algebra Constructors)

**Agent:** polyglot-implementation-engineer
**Feature:** Phase 2.1 (algebra_constructors.jl refactoring)

### Actions
- Refactored `algebra_constructors.jl` to use FastPolynomials `create_*_variables` functions
- Replaced `@ncpolyvar` macro with registry-based variable creation
- Removed manual equality constraints (algebra simplification handles Pauli rules automatically)
- Removed `is_unipotent`, `is_projective`, `comm_gps` (inferred from algebra type)
- Added six algebra constructors:
  - `pauli_algebra(N)` - returns `(registry, sx, sy, sz)`
  - `fermionic_algebra(N)` - returns `(registry, a, a_dag)`
  - `bosonic_algebra(N)` - returns `(registry, c, c_dag)`
  - `projector_algebra(prefix, N)` and `projector_algebra(prefixes, N)` - returns `(registry, projectors)`
  - `unipotent_algebra(prefix, N)` and `unipotent_algebra(prefixes, N)` - returns `(registry, variables)`
  - `noncommutative_algebra(prefix, N)` and `noncommutative_algebra(prefixes, N)` - returns `(registry, variables)`
- Updated `NCTSSoS.jl`:
  - Added imports for `create_*_variables` functions from FastPolynomials
  - Exported new algebra constructors

### Files Modified
- `src/algebra_constructors.jl` (complete rewrite)
- `src/NCTSSoS.jl` (imports and exports)

### Verification
- Module compiles without errors
- `make test-FastPoly` passes (1141 tests)
- All algebra constructors tested and working:
  - Correct registry types (PauliAlgebra, FermionicAlgebra, etc.)
  - Correct monomial types with proper index types
  - Multi-prefix variants work correctly

### Outcome
Phase 2 complete. All algebra constructors now use FastPolynomials `create_*_variables` API:
- No `@ncpolyvar` usage in algebra_constructors.jl
- No manual equality constraints for algebra rules (handled by simplification)
- Registry-based API: `(registry = ..., vars = ...)`

### Next Steps
- Phase 3: Interface refactoring (interface.jl, sos_solver.jl, gns.jl)
- Phase 4: Legacy code removal
- Phase 5: Test migration

**Commit:** `cfcf4ba` - refactor(algebra): registry-based algebra constructors

---

## Session: 2025-12-12 - Phase 1.3 Complete (Unified Symbolic MomentProblem)

**Agent:** polyglot-implementation-engineer
**Feature:** Phase 1.3 (moment_solver.jl + complex_moment_solver.jl merge)

### Actions
- Merged `moment_solver.jl` and `complex_moment_solver.jl` into unified symbolic design
- Created `MomentProblem{A<:AlgebraType, T<:Integer, M<:Monomial{A,T}, P<:Polynomial{A,T}}`:
  - `objective::P` - polynomial objective function
  - `constraints::Vector{Tuple{Symbol, Matrix{P}}}` - symbolic constraint matrices
  - `total_basis::Vector{M}` - union of basis monomials
- Implemented unified `moment_relax()` that dispatches on `_is_complex_problem(A)`:
  - Real algebras (NonCommutative, Projector, Unipotent): `:PSD` cone
  - Complex algebras (Pauli, Fermionic, Bosonic): `:HPSD` cone
- Added `_build_constraint_matrix()` for symbolic constraint construction
- Added `solve_moment_problem()` interface for direct solving (with `_solve_real_moment_problem` and `_solve_complex_moment_problem` helpers)
- Updated `sos_solver.jl` with unified `sos_dualize()`:
  - `_sos_dualize_real()` for real algebras
  - `_sos_dualize_hermitian()` for complex algebras with 2n x 2n embedding
- Added backward compatibility:
  - `complex_moment_relax()` as alias for `moment_relax()`
  - `ComplexMomentProblem` type alias
- Updated `interface.jl` to use unified `moment_relax()`
- Deleted `src/complex_moment_solver.jl` (merged into moment_solver.jl)

### Files Modified
- `src/moment_solver.jl` (major refactor - unified symbolic MomentProblem)
- `src/sos_solver.jl` (refactored for unified type)
- `src/interface.jl` (simplified dispatch logic)
- `src/NCTSSoS.jl` (removed complex_moment_solver.jl include)

### Files Deleted
- `src/complex_moment_solver.jl` (merged into moment_solver.jl)

### Verification
- Module compiles without errors
- `make test-FastPoly` passes (1141 tests)

### Outcome
Phase 1.3 complete. The new `MomentProblem{A,T,M,P}` is a purely symbolic representation:
- No JuMP model in the struct (moved to `solve_moment_problem()`)
- Constraints as `Vector{Tuple{Symbol, Matrix{P}}}` with cone types (:Zero, :PSD, :HPSD)
- Unified `moment_relax()` handles both real and complex via `_is_complex_problem(A)` trait
- `sos_dualize()` handles both real and Hermitian dualization

### Next Steps
- Phase 2: Refactor algebra_constructors.jl
- Re-enable tests for completed phases

**Commit:** `1a3c257` - refactor(moment): unified symbolic MomentProblem{A,T,M,P}

---

## Session: 2025-12-12 - Checkpoint (Phases 1.1-1.2 Complete)

**Agent:** orchestrator
**Feature:** N/A (checkpoint)

### Actions
- Completed Phase 1.1: pop.jl refactored with algebra type parameters
- Completed Phase 1.2: sparse.jl refactored with registry-based cliques
- Updated plan.md with progress markers
- Discussed test re-enablement strategy with user

### Outcome
Phase 1 is 2/3 complete. Core types refactored with algebra parameters. Module compiles and FastPolynomials tests pass (1141 tests).

### Current State
- Phase 1.1 (pop.jl): ✓ Complete
- Phase 1.2 (sparse.jl): ✓ Complete
- Phase 1.3 (moment solvers): Pending
- Tests: Commented out, will re-enable after Phase 2

### Next Steps
- Phase 1.3: Merge and refactor moment_solver.jl + complex_moment_solver.jl
- Then Phase 2: Algebra constructors
- Then re-enable tests (pop.jl, algebra_constructors.jl)

**Commit:** (checkpoint push)

---

## Session: 2025-12-12 - Phase 1.2 Complete (sparse.jl refactored)

**Agent:** polyglot-implementation-engineer
**Feature:** Phase 1.2 (sparse.jl registry-based cliques and index mapping)

### Actions
- Added `subregistry()` function to FastPolynomials (`src/FastPolynomials/src/variable_registry.jl`)
- Added `variable_indices()` for Monomial type (`src/FastPolynomials/src/polynomial.jl`)
- Refactored `CorrelativeSparsity` struct with new type signature:
  - `CorrelativeSparsity{A<:AlgebraType, T<:Integer, P<:Polynomial{A,T}, M<:Monomial{A,T}}`
  - `cliques::Vector{Vector{T}}` (stores indices instead of Variables)
  - `registry::VariableRegistry{A,T}` (for symbol lookups)
- Updated `get_correlative_graph` to use `variable_indices()` and position-based mapping
- Updated `assign_constraint` to work with index-based cliques
- Replaced `get_basis(clique, order)` with `subregistry()` + `get_ncbasis()` pattern
- Added `extract_monomials_from_basis()` helper for basis polynomial -> monomial extraction
- Added `clique_variables()` function for backward compatibility
- Updated `PolyOptResult` type signature in interface.jl
- Updated `cs_nctssos` to use `variable_indices()` for clique objective computation
- Updated `moment_relax` and `complex_moment_relax` signatures to match new CorrelativeSparsity type

### Files Modified
- `src/FastPolynomials/src/variable_registry.jl` (+subregistry function)
- `src/FastPolynomials/src/polynomial.jl` (+variable_indices for Monomial)
- `src/FastPolynomials/src/FastPolynomials.jl` (exports)
- `src/NCTSSoS.jl` (imports)
- `src/sparse.jl` (major refactor)
- `src/interface.jl` (PolyOptResult type, cs_nctssos logic)
- `src/moment_solver.jl` (type signature)
- `src/complex_moment_solver.jl` (type signature)

### Verification
- `make test-FastPoly` passes (1141 tests)
- Module compiles without errors

### Outcome
Phase 1.2 complete. The new `CorrelativeSparsity{A,T,P,M}` type with:
- Index-based cliques (Vector{Vector{T}})
- Registry reference for symbol lookups
- Position-based graph construction using `variable_indices()`
- Basis generation via `subregistry()` + `get_ncbasis()`

### Next Steps
- Phase 1.3: Merge and refactor moment_solver.jl + complex_moment_solver.jl

**Commit:** `a347588` - refactor(sparse): registry-based cliques and index mapping

---

## Session: 2025-12-12 - Phase 1.1 Complete (pop.jl refactored)

**Agent:** polyglot-implementation-engineer
**Feature:** F001 (Phase 1.1 pop.jl refactoring)

### Actions
- Refactored `src/pop.jl` with new type signatures:
  - `OptimizationProblem{A<:AlgebraType, P}` (was `OptimizationProblem{P}`)
  - `PolyOpt{A<:AlgebraType, P<:Polynomial{A}}` (was `PolyOpt{P}`)
  - Added `registry::VariableRegistry{A}` field (replaces `variables::Vector{Variable}`)
  - Removed `comm_gps`, `is_unipotent`, `is_projective` fields
- Added algebra trait function `_is_complex_problem(::Type{A})` for dispatch
- Added backward compatibility:
  - `ComplexPolyOpt` as alias for `PolyOpt`
  - `cpolyopt` as alias for `polyopt`
  - `getproperty` accessor for `variables`, `is_unipotent`, `is_projective`, `comm_gps`
- Updated type signatures in dependent files (minimal changes to compile):
  - `src/sparse.jl`: `correlative_sparsity` function signature
  - `src/interface.jl`: `cs_nctssos`, `cs_nctssos_higher` function signatures
  - `src/moment_solver.jl`: `moment_relax` function signature
  - `src/complex_moment_solver.jl`: renamed to `complex_moment_relax`
- Updated `src/NCTSSoS.jl` with new imports (VariableRegistry, symbols, indices, index_type, algebra_type)

### Files Modified
- `src/pop.jl` (major refactor)
- `src/NCTSSoS.jl` (imports)
- `src/sparse.jl` (type signature only)
- `src/interface.jl` (type signature + trait dispatch)
- `src/moment_solver.jl` (type signature only)
- `src/complex_moment_solver.jl` (type signature + rename function)

### Verification
- `make test-FastPoly` passes (1141 tests)
- Module compiles without errors

### Outcome
Phase 1.1 complete. The new `PolyOpt{A, P}` type with algebra parameter and registry-based API is in place. Backward compatibility accessors allow existing code paths to work during transition.

### Next Steps
- Phase 1.2: Full refactor of sparse.jl
- Phase 1.3: Full refactor and merge of moment_solver.jl + complex_moment_solver.jl

**Commit:** `3b38219` - refactor(pop): add algebra type parameters and registry-based API

---

## Session: 2025-12-12 - Phases 4-5 Planning Approved

**Agent:** orchestrator
**Feature:** N/A (planning)

### Actions
- Reviewed Phase 4 (Legacy Code Removal) plan
- Reviewed Phase 5 (Test Migration) plan
- Clarified 3 design decisions via user Q&A

### Decisions Made
1. **Legacy API**: Remove completely - no deprecation warnings, clean break
2. **Test Migration**: Comment out ALL tests in `runtests.jl`, re-enable incrementally as phases complete
3. **AbstractPolynomial**: Keep union type for flexibility/future extensions

### Outcome
All 5 phases planned and approved. Ready for implementation.

### Next Steps
- Begin Phase 1 Step 1.1: Refactor pop.jl
- Comment out tests in runtests.jl (Phase 5.0)

**Commit:** `059a7b5` - docs(fastpoly-integration): complete planning (all 5 phases approved)

---

## Session: 2025-12-12 - Planning Complete (Phases 1-3 Approved)

**Agent:** orchestrator
**Feature:** N/A (planning)

### Actions
- Created task folder structure with 20 features
- Conducted detailed step-by-step plan review with user
- Approved design decisions through Q&A

### Key Design Decisions Approved

**Phase 1: Core Type Refactoring**
- Step 1.1 (pop.jl): Add `A<:AlgebraType` to types, merge `PolyOpt`+`ComplexPolyOpt`, use `VariableRegistry`, remove `comm_gps`/`is_unipotent`/`is_projective`
- Step 1.2 (sparse.jl): Add `A` to `CorrelativeSparsity`, use `subregistry()`, `variable_indices()`, position-based graph mapping, signed indices are distinct
- Step 1.3 (moment_solver.jl): Single symbolic `MomentProblem{A,T,M,P,K}`, basis is `Vector{Polynomial}`, merge with complex_moment_solver.jl

**Phase 2: Algebra Constructors**
- Remove `@ncpolyvar`, use `create_*_variables`
- Remove manual equality constraints (algebra handles simplification)
- Return simplified NamedTuple: `(registry=..., vars=...)`

**Phase 3: Interface Refactoring**
- Add `ProblemKind` hierarchy (`RealProblem`, `HermitianProblem`)
- Dispatch `sos_dualize` on `K` type parameter
- `gns.jl`: `Vector{Variable}` → `VariableRegistry`, return `Vector{Matrix}`

### Outcome
Phases 1-3 planning complete and approved

### Next Steps
- Review Phase 4: Legacy code removal
- Review Phase 5: Test migration
- Then proceed to implementation

**Commit:** (pending - plan files)

---

## Session: 2025-12-12 - Initialization

**Agent:** orchestrator
**Feature:** N/A (task setup)

### Actions
- Created task folder structure
- Defined 20 features in features.json
- Clarified requirements with user

### Outcome
Task initialized

---
