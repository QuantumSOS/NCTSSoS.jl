# NCTSSoS.jl Core Source Code Review Report

## Executive Summary
The core optimization framework (`src/`) is well-structured, modular, and effectively leverages the `FastPolynomials` backend. The architecture separates problem definition (`pop.jl`), solver interface (`interface.jl`), and algorithmic logic (`sparse.jl`, `sos_solver.jl`), promoting maintainability.

## Detailed component Analysis

### 1. Architecture & Data Structures (`NCTSSoS.jl`, `pop.jl`)
*   **Strengths**:
    *   Clean type hierarchy with `OptimizationProblem`, `PolyOpt`, and `StatePolyOpt`.
    *   Effective use of parametric types (`A<:AlgebraType`, `P<:Polynomial`) for type stability.
    *   `_is_complex_problem` trait offers a robust mechanism for dispatching between Real and Hermitian SDP relaxations.
*   **Observations**:
    *   The `polyopt` constructor performs `unique!(copy(eq_constraints))`. For problems with thousands of constraints, this could be a minor bottleneck if polynomial hashing is expensive.
    *   Explicitly non-exported types (`Monomial`, `Polynomial`) prevent namespace pollution but require users to be aware of `NCTSSoS.FastPolynomials`.

### 2. Interface & Solver API (`interface.jl`)
*   **Strengths**:
    *   `SolverConfig` provides a flexible configuration object for passing optimizer settings and algorithm choices.
    *   The separation of `cs_nctssos` (first pass) and `cs_nctssos_higher` (refinement) is a good design for iterative optimization.
    *   `reduce(+, ...)` usage for partial objectives is memory-efficient.

### 3. Core Solvers (`sos_solver.jl`, `moment_solver.jl`)
*   **Strengths**:
    *   **Hermitian Embedding**: The `_sos_dualize_hermitian` function correctly implements the standard `[Re -Im; Im Re]` embedding for complex SDPs.
    *   **Fermionic Parity**: `_add_parity_constraints!` in `moment_solver.jl` is a critical correctness feature, automatically enforcing superselection rules.
    *   **Dualization**: `get_Cαj` efficiently extracts coefficients for the SOS dual problem using dictionary mapping.
*   **Risks**:
    *   `solve_moment_problem` instantiates a symbolic moment matrix. For large relaxation orders ($d \ge 3$), this matrix can grow very quickly. The SOS dual approach (`sos_dualize`) is preferred and correctly prioritized in docs.

### 4. Sparsity & Optimization (`sparse.jl`, `elimination.jl`)
*   **Strengths**:
    *   **Graph Theory**: `get_correlative_graph` and `clique_decomp` robustly handle problem decomposition.
    *   **Term Sparsity**: The iterative TS logic (`term_sparsities`) is implemented correctly following the CS-TSSOS methodology.
*   **Optimization Opportunities**:
    *   **Bottleneck**: `get_term_sparsity_graph` contains a nested loop over basis elements (O(N²)) that performs symbolic multiplication and simplification (`_neat_dot3`) inside the loop.
        *   *Recommendation*: This is a prime candidate for parallelization (`Threads.@threads`) or memoization if bases are reused.

### 5. GNS Reconstruction (`gns.jl`)
*   **Strengths**:
    *   Implements the "Flat Extension" check (`rank(H) == rank(H_block)`), providing user feedback if reconstruction is mathematically invalid.
    *   Clear separation of "localizing matrix" construction.
*   **Observations**:
    *   SVD is performed on dense matrices. For extremely large moment matrices, an iterative eigensolver might be needed, though GNS is typically a post-processing step on smaller blocks.

## Recommendations
1.  **Performance**: Investigate parallelizing the edge-checking loop in `src/sparse.jl:get_term_sparsity_graph`.
2.  **Safety**: Add a check in `polyopt` to ensure `objective` and `constraints` share the same `VariableRegistry`.
3.  **Documentation**: Add a "Performance Tips" section to docs explaining the memory implications of `solve_moment_problem` vs `dualize=true`.
