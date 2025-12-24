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
*   **`src/interface.jl`** *(Reviewed - Changes Requested)*
    *   [x] Evaluate public API surface.
    *   [x] Check keyword argument handling and defaults.
    *   [x] Specific review of error messages and validation logic.
    *   **Pending Fixes**:
        *   [ ] **High**: Refactor lines 99-105 - inefficient partial objective calculation (allocates zero terms)
        *   [ ] **Medium**: Optimize order calculation (lines 93, 246) - use `Iterators.flatten` instead of array concatenation
        *   [ ] **Low**: Add solver status checking before accessing `objective_value` (line 123)
        *   [ ] **Low**: Consider typing `optimizer` field in `SolverConfig` (line 61)

### Phase 3: Core Solver Logic
**Focus**: SDP construction, solver interaction, and numerical correctness.
*   **`src/sos_solver.jl`**
    *   [ ] Audit SDP relaxation construction (Hierarchy levels).
    *   [ ] Review interaction with JuMP or low-level solver interfaces.
*   **`src/moment_solver.jl`**
    *   [ ] Analyze moment matrix construction.
    *   [ ] Check dual variable handling and extraction.

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
- [ ] Phase 2 Complete
- [ ] Phase 3 Complete
- [ ] Phase 4 Complete
- [ ] Phase 5 Complete
