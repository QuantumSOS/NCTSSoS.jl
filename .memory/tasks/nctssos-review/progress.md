# Review Progress Log

## Session 1: 2026-01-07

**Focus**: Setup and planning

**Actions**:
- Explored codebase structure (3 agents in parallel)
- Created review plan with 6 phases
- Established session tracking structure

**Files reviewed**: None (planning only)

**Outcome**: Plan approved, ready to start Phase 1

**Next step**: Review `src/optimization/problem.jl` (Phase 1.1)

---

## Session 2: 2026-01-08

**Focus**: Phase 1 - Optimization Pipeline

**Files reviewed**:
| File | Lines | Key Findings |
|------|-------|--------------|
| `problem.jl` | 310 | Clean types, good docs, algebra type as first-class param |
| `interface.jl` | 467 | Clear pipeline, some code duplication, unused var |
| `sparsity.jl` | 984 | Correlative + term sparsity, matches NCTSSOS algorithm |
| `moment.jl` | 841 | Symbolic representation, Hermitian embedding correct |
| `sos.jl` | 460 | SOS dualization, sparse coefficient extraction |

**Key Architecture Insights**:
```
PolyOpt → correlative_sparsity() → cliques
       → init_activated_supp() → term_sparsities() → blocks
       → moment_relax() → MomentProblem (symbolic)
       → sos_dualize() → SOSProblem (JuMP)
       → optimize!() → PolyOptResult
```

**Quality Summary**:
- Documentation: ✓✓✓ Excellent throughout
- Type safety: ✓✓✓ Parametric types enforce consistency
- Correctness: ✓✓✓ Algorithms match NCTSSOS paper
- Code structure: ✓✓ Duplication between Polynomial/StatePolynomial

**Issues Found**:
1. Code duplication (~200-300 lines per file) between PolyOpt and StatePolyOpt paths
2. `Matrix{Any}` type instability in complex moment solving (noted in code)
3. No solver status check after `optimize!()` in interface.jl
4. Unused variable `M = NCStateWord{...}` in interface.jl:453

**Annotations created**:
- `.memory/tasks/nctssos-review/annotations/problem.md`
- `.memory/tasks/nctssos-review/annotations/interface.md`
- `.memory/tasks/nctssos-review/annotations/sparsity.md`
- `.memory/tasks/nctssos-review/annotations/moment.md`
- `.memory/tasks/nctssos-review/annotations/sos.md`

**Outcome**: Phase 1 complete

**Next step**: Phase 2 - Foundation types (algebra.jl, monomial.jl, etc.)

---

## Session 3: 2026-01-08 (Interactive Redo)

**Focus**: Phase 1 redo - Interactive review with user

**Actions**:
- Re-reviewed `problem.jl` interactively
  - Discussed StatePolyOpt necessity → noted as MAINT-001 (duplication issue)
- Started `interface.jl` review
  - Identified redundant `_count_unique_moment_elements` computation
  - User requested immediate fix → created branch `ys/dry-moment-stats`
  - Refactored: n_unique now computed during sos_dualize/solve_moment_problem
  - Removed ~70 lines of redundant nested loops
  - All tests passing (minimal: 24, polynomials: 2005, relaxations: 94)
  - Merged back to `ys/type-system-revamp`

**Commits**:
- `f461966` refactor: DRY n_unique_moment_elements computation

**Issues Found/Fixed**:
- MAINT-001: StatePolyOpt/PolyOpt duplication (noted, not fixed)
- FIXED: Redundant _count_unique_moment_elements recomputation
- FIXED: Unused variable `M = NCStateWord{...}` removed

**Files modified**:
- `src/optimization/interface.jl` (-67 lines)
- `src/optimization/moment.jl` (+9 lines)
- `src/optimization/sos.jl` (+3 lines)

**Outcome**: DRY refactor complete, review paused

**Next step**: Continue `interface.jl` review, then `sparsity.jl`

---

## Session 4: 2026-01-08

**Focus**: Extract testable functions from `cs_nctssos` + add unit tests

**Extractions completed**:
1. `compute_relaxation_order(pop, user_order)` - unified order computation (DRY: 2→1)
2. `solve_sdp(moment_problem, optimizer; dualize)` + `_check_solver_status` - **fixes missing error handling bug** (DRY: 3→1)
3. `project_to_clique(poly, clique_indices)` - 2 dispatched methods for Polynomial/NCStatePolynomial (DRY: 2→1 per type)

**Bugs fixed**:
- Solver status was never checked after `optimize!()` - now throws on INFEASIBLE/UNBOUNDED/NUMERICAL_ERROR
- `compute_relaxation_order` now returns at least 1 for trivial polynomials (degree 0) via `max(1, ...)`

**Unit tests added** (`test/relaxations/interface.jl`):
- `compute_relaxation_order`: auto-compute, user override, constraints, trivial polynomial
- `project_to_clique`: partial/empty/full projection, NCStatePolynomial
- `_check_solver_status`: acceptable status constants, integration test

**Tests**: All passing
- Minimal: 24/24
- Polynomials: 2005/2005
- Relaxations: 114/114 (+20 new unit tests)

**Commits**:
- `1212753` refactor(interface): extract testable functions
- `f9caecb` test: add unit tests for extracted interface functions
- `5d1d8b1` docs: update review progress

**Next step**: Phase 1.3 - Correlative Sparsity (`sparsity.jl`)

**Status**: COMPLETE

---

## Session 5: 2026-01-08

**Focus**: Extract `compute_sparsity` for debugging sparsity without solving

**Problem**: User wanted to inspect `initial_activated_supps` and `cliques_term_sparsities` for debugging, but these were only accessible after solving the SDP.

**Solution**: Introduced `SparsityResult` and `StateSparsityResult` structs with a `compute_sparsity()` function that returns sparsity info without running the solver.

**Changes**:
1. New types: `SparsityResult{A,TI,P,M}`, `StateSparsityResult{A,ST,TI,P,M}`
2. New function: `compute_sparsity(pop, solver_config)` - returns sparsity before solving
3. Refactored `PolyOptResult` and `StatePolyOptResult` to store nested `sparsity` field
4. Simplified `cs_nctssos`: now calls `compute_sparsity` → `moment_relax` → `solve_sdp`
5. Updated `cs_nctssos_higher` to use new field paths
6. Exported: `SparsityResult`, `StateSparsityResult`, `compute_sparsity`

**API Change** (breaking):
```julia
# Before
result.corr_sparsity
result.cliques_term_sparsities

# After
result.sparsity.corr_sparsity
result.sparsity.initial_activated_supps  # NEW - now accessible
result.sparsity.cliques_term_sparsities
```

**New Debug Usage**:
```julia
# Inspect sparsity without solving
sparsity = compute_sparsity(pop, config)
@show sparsity.initial_activated_supps
@show sparsity.cliques_term_sparsities

# Full solve (unchanged API)
result = cs_nctssos(pop, config)
```

**Tests**: All passing (minimal: 24, relaxations: 114)

**Next step**: Phase 1.3 - Correlative Sparsity (`sparsity.jl`)

**Status**: COMPLETE

---

## Session 6: 2026-01-08

**Focus**: Phase 1.3/1.4 - Correlative & Term Sparsity deep dive

**Files reviewed**: `src/optimization/sparsity.jl` (985 lines)

**Architecture traced**:
```
correlative_sparsity(pop, order, elim_algo)
   → get_correlative_graph() → SimpleGraph (variables connected if in same monomial)
   → clique_decomp(G, elim_algo) → cliques (via CliqueTrees.jl)
   → assign_constraint() → clq_cons, global_cons
   → get_ncbasis() per clique → clq_mom_mtx_bases
   → CorrelativeSparsity

term_sparsities(activated_supp, cons, bases, ...)
   → get_term_sparsity_graph() → edges where bᵢ†·supp·bⱼ ∈ activated_supp
   → clique_decomp → blocks
   → TermSparsity{block_bases}
```

**Key Insights**:
1. Correlative sparsity: Variables that never appear together → don't interact → can be in separate SDP blocks
2. Term sparsity: Within a clique, basis pairs whose products are outside support → zero matrix entries → block structure
3. Objective uses per-monomial edges (finer sparsity), constraints use whole-polynomial edges (ensures assignability)

**Bug Found & Fixed** (BUG-001):
- Location: `get_state_correlative_graph()` lines 661-666
- Issue: Used per-monomial constraint handling (should be whole constraint)
- Impact: State polynomial constraints could become global unnecessarily
- Fix: Changed to `variable_indices(con)` for whole constraint
- Commit: `5dcbae1` fix(sparsity): use whole constraint for state correlative graph

**Review Notes** (to investigate):
- REVIEW-001: `abs(idx)` in get_positions (line 115) - treats a₁† and a₁ as same variable
- REVIEW-002: `init_activated_supp` includes diagonals for Polynomial but not State (monosquare=false)

**Refactoring Opportunities**:
- MAINT-002: `simplified_monomials()` helper to avoid `Polynomial(simplify(...))` overhead
- MAINT-003: Unify Polynomial/State sparsity functions (~400 lines duplication)

**Tests**: All passing after bug fix (minimal: 24, relaxations: 94)

**Next step**: Monomial hierarchy revamp (prerequisite for MAINT-002/MAINT-003)

**Status**: REVIEW COMPLETE, REFACTORING BLOCKED (see `plan-sparsity-refactor.md`)

---

## Checklist Progress

### Phase 1: Optimization Pipeline (Interactive Redo)
- [x] 1.1 Problem Definition (`problem.jl`) - reviewed StatePolyOpt duplication issue
- [x] 1.2 Main Solver Entry (`interface.jl`) - extracted testable functions, fixed solver status bug
- [x] 1.3 Correlative Sparsity (`sparsity.jl`) - fixed BUG-001, identified unification opportunity
- [x] 1.4 Term Sparsity (`sparsity.jl:460-600`) - traced algorithm, noted monosquare asymmetry
- [ ] 1.5 Moment Relaxation (`moment.jl`)
- [ ] 1.6 SOS Dualization (`sos.jl`)

### Phase 2: Foundation
- [ ] 2.1 AlgebraType Hierarchy (`algebra.jl`)
- [ ] 2.2 Monomial + Term + Polynomial
- [ ] 2.3 Variable Registry (`registry.jl`)

### Phase 3: Simplification Rules
- [ ] 3.1 Dispatch Architecture
- [ ] 3.2 Algebra-Specific Rules (6 algebras)

### Phase 4: Canonicalization + Algorithms
- [ ] Canonicalization (`canonicalization.jl`)
- [ ] Basis generation (`basis.jl`)

### Phase 5: State Polynomials
- [ ] StateSymbol, StateWord, StatePolynomial

### Phase 6: Test Coverage Audit
- [ ] Run test suites
- [ ] Coverage analysis
- [ ] Test quality review
