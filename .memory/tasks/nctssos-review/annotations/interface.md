# Review: src/optimization/interface.jl

**Lines**: 467
**Purpose**: Main solver entry point and result types

## Structure

```
PolyOptResult{T,A,TI,P,M}       # Standard result
StatePolyOptResult{T,A,ST,TI,P,M}  # State polynomial result
SolverConfig                    # optimizer + order + cs_algo + ts_algo

cs_nctssos(pop, config)         # Main solver
cs_nctssos_higher(pop, prev, config)  # Higher-order refinement
```

## Key Pipeline (cs_nctssos, lines 177-221)

```
1. Auto-compute order from max degree: ceil(max_deg / 2)
2. correlative_sparsity(pop, order, cs_algo) → clique decomposition
3. Compute partial objectives per clique (filter terms by variable indices)
4. init_activated_supp() → initial activated supports
5. term_sparsities() → compute term sparsity per clique
6. moment_relax() → create symbolic moment problem
7. Either:
   - dualize=true: sos_dualize() → solve SOS dual
   - dualize=false: solve_moment_problem() directly
8. Return PolyOptResult with metrics
```

## Result Struct Fields

- `objective::T` - optimal value
- `corr_sparsity` - clique structure
- `cliques_term_sparsities` - term sparsity blocks per clique
- `model` - JuMP model (for debugging/extraction)
- `moment_matrix_sizes` - block sizes
- `n_unique_moment_matrix_elements` - SDP variable count

## Quality Assessment

| Aspect | Rating | Notes |
|--------|--------|-------|
| Documentation | ✓✓✓ | Clear docstrings with descriptions |
| Type safety | ✓✓✓ | Parametric types throughout |
| Code structure | ✓✓ | Some duplication between PolyOpt/StatePolyOpt versions |
| Error handling | ✓ | Limited - relies on downstream errors |

## Observations

**Strengths**:
1. Clear pipeline separation: sparsity → relaxation → solve
2. Higher-order refinement (`cs_nctssos_higher`) builds on previous results
3. Metrics (`moment_matrix_sizes`, `n_unique_moment_matrix_elements`) useful for debugging

**Concerns**:
1. **Code duplication** - `cs_nctssos` for PolyOpt and StatePolyOpt nearly identical (~80% overlap). Could factor shared logic.

2. **Order=0 handling** (line 179-184):
   ```julia
   max_deg = maximum(degree(poly) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints])
   ```
   If all polys are empty/trivial, this might fail or produce unexpected results.

3. **Partial objective computation** (lines 190-197):
   ```julia
   cliques_objective = map(corr_sparsity.cliques) do clique_indices
       clique_set = Set(clique_indices)
       reduce(+, ...)
   ```
   Creates Set per clique. Could be more efficient for large problems, but likely fine.

4. **No solver status check** - After `optimize!(model)`, doesn't check termination status before returning `objective_value()`. Silent failures if solver fails.

5. **Unused variable** in state version (line 453):
   ```julia
   M = NCStateWord{ST,A,T}  # Never used after assignment
   ```

## Data Flow Insight

```
PolyOpt → correlative_sparsity() → cliques
       ↘
         init_activated_supp() → term_sparsities() → block bases
                                                  ↓
                           moment_relax() → MomentProblem
                                                  ↓
                           sos_dualize() → SOSProblem
                                                  ↓
                           optimize!() → PolyOptResult
```

## Code Quality: Good

Clean pipeline, but some minor issues (solver status, duplication, unused var).
