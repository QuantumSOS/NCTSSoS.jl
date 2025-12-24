# NCTSSoS.jl Core Review Plan

## Objective
Systematic review of the NCTSSoS.jl optimization framework (excluding the `FastPolynomials` submodule). The goal is to ensure code quality, performance optimization, architectural soundness, and adherence to Julia best practices.

## Scope
**Target Directory**: `src/`
**Excluded**: `src/FastPolynomials/`

## Review Phases

### Phase 1: Architecture & Data Structures
**Focus**: Module structure, type definitions, and problem representation.
*   **`src/NCTSSoS.jl`**
    *   [x] Analyze module exports and dependency management.
    *   [x] Check for global state or constants.
*   **`src/pop.jl`** (Polynomial Optimization Problem)
    *   [x] Review `POP` struct definition.
    *   [x] Examine how objective functions and constraints are stored.
    *   [x] Verify type stability of core structures.

### Phase 2: User Interface & API
**Focus**: Usability, multiple dispatch, and input validation.
*   **`src/interface.jl`** *(Reviewed - Complete)*
    *   [x] Evaluate public API surface.
    *   [x] Check keyword argument handling and defaults.
    *   [x] Specific review of error messages and validation logic.
    *   **Fixes Applied**:
        *   [x] **Critical**: Fixed runtime crash on empty objective - added `init=zero(pop.objective)` to `reduce` (lines 99-106)
        *   [x] **Warning**: Replaced splat with `reduce(sorted_union, ...)` to avoid stack issues (line 158)

### Phase 3: Core Solver Logic
**Focus**: SDP construction, solver interaction, and numerical correctness.
*   **`src/sos_solver.jl`** *(Reviewed - Complete)*
    *   [x] Audit SDP relaxation construction (Hierarchy levels).
    *   [x] Review interaction with JuMP or low-level solver interfaces.
    *   **Fixes Applied**:
        *   [x] **Critical**: Fixed repeated sorting in `_sos_dualize_real` - sort basis once before loop (lines 200-211)
*   **`src/moment_solver.jl`** *(Reviewed - Complete)*
    *   [x] Analyze moment matrix construction.
    *   [x] Check dual variable handling and extraction.
    *   **Fixes Applied**:
        *   [x] **Warning**: Fixed typo `rol_idx` -> `row_idx` (lines 145, 645)

### Phase 4: Optimization & Sparsity
**Focus**: Algorithmic efficiency, graph theory implementations, and linear algebra.
*   **`src/sparse.jl`**
    *   [ ] Review Correlative Sparsity and Term Sparsity algorithms.
    *   [ ] Check clique decomposition efficiency.
*   **`src/elimination.jl`**
    *   [ ] Verify variable elimination/reduction logic.
    *   [ ] Check for potential bottlenecks in symbolic manipulation.

### Phase 5: Specialized Features
**Focus**: Quantum state features and GNS construction.
*   **`src/gns.jl`**
    *   [ ] Review GNS construction logic for state extraction.
    *   [ ] Ensure compatibility with `FastPolynomials` algebra types.

## General Review Criteria
1.  **Performance**: Identification of type instabilities, unnecessary allocations, or suboptimal algorithms.
2.  **Style**: Adherence to Blue Style and standard Julia idioms.
3.  **Integration**: Correct usage of `FastPolynomials` types and methods.
4.  **Testing**: Cross-reference with `test/` to ensure critical paths are covered.

## Progress Tracking
- [x] Phase 1 Complete
- [x] Phase 2 Complete
- [x] Phase 3 Complete
- [ ] Phase 4 Complete
- [ ] Phase 5 Complete
