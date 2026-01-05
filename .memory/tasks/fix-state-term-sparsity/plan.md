# Plan: Fix State Polynomial Term Sparsity SOS Dualization

## Problem Summary

State polynomial optimization with term sparsity produces incorrect results (e.g., bilocal network test: expected -4.0, actual -1.76M). The root cause is that the SOS dualization loses StateWord contributions when term sparsity creates multiple smaller blocks.

## Root Cause Analysis

The issue is in `_sos_dualize_state` (src/optimization/sos.jl:252-330):

1. **Constraint matrices are built per term sparsity block** (moment.jl:678-682):
   - Each block has dimensions `length(block_basis) × length(block_basis)`
   - Matrix (row, col) indices are local to the block: `1:length(block_basis)`

2. **Coefficient extraction searches the full basis** (sos.jl:310-323):
   - `_get_state_Cαj(sorted_unsym_basis, mat)` searches in `mp.total_basis`
   - Returns `(basis_idx, row, col)` where `row, col` are valid for the matrix
   - But `basis_idx` points to position in full basis, not block basis

3. **The mismatch**:
   - `sorted_unsym_basis[basis_idx]` gives the NCStateWord
   - But multiple blocks may contain the same StateWord at different (row, col) positions
   - Only one contribution gets recorded, losing others

## Solution

Refactor `StateMomentProblem` to store block basis information with each constraint matrix, then modify `_sos_dualize_state` to:

1. Loop over blocks and compute StateWord contributions directly from (row, col) pairs
2. Accumulate contributions across all blocks like NCTSSOS does

### Step 1: Modify `StateMomentProblem` struct

**File:** `src/optimization/moment.jl:590-594`

Change constraint storage from:
```julia
struct StateMomentProblem{A<:AlgebraType, ST<:StateType, T<:Integer, M<:NCStateWord{ST,A,T}, P<:NCStatePolynomial}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}
    total_basis::Vector{M}
end
```

To:
```julia
struct StateMomentProblem{A<:AlgebraType, ST<:StateType, T<:Integer, M<:NCStateWord{ST,A,T}, P<:NCStatePolynomial}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}, Vector{M}}}  # Add block_basis
    total_basis::Vector{M}
end
```

### Step 2: Update `moment_relax` for state polynomials

**File:** `src/optimization/moment.jl:674-694`

When building constraint matrices, store the block basis:
```julia
for ts_sub_basis in term_sparsity.block_bases
    cone = poly in pop.eq_constraints ? :Zero : psd_cone
    mat = _build_state_constraint_matrix(poly, ts_sub_basis, cone)
    # Store block_basis with constraint
    push!(constraints, (mat[1], mat[2], ts_sub_basis))
end
```

### Step 3: Rewrite `_sos_dualize_state`

**File:** `src/optimization/sos.jl:252-330`

Replace the coefficient extraction logic (lines 308-324) with direct block iteration:

```julia
# Process each constraint matrix with its block basis
for (i, (cone, mat, block_basis)) in enumerate(mp.constraints)
    dim = size(mat, 1)
    # Iterate over upper triangle of this block's matrix
    for row in 1:dim
        for col in row:dim
            poly = mat[row, col]
            for (coeff, ncsw) in zip(coefficients(poly), monomials(poly))
                # Convert NCStateWord to StateWord
                sw = symmetric_canon(expval(ncsw))
                sw_idx = searchsortedfirst(state_basis, sw)
                if sw_idx <= n_basis && state_basis[sw_idx] == sw
                    # Multiply off-diagonal by 2 for symmetric contribution
                    effective_coeff = (row == col) ? coeff : 2 * coeff
                    add_to_expression!(fα_constraints[sw_idx], -effective_coeff, dual_variables[i][row, col])
                end
            end
        end
    end
end
```

### Step 4: Remove `_get_state_Cαj` function

**File:** `src/optimization/sos.jl:332-366`

This function is no longer needed after the refactor.

## Files to Modify

1. `src/optimization/moment.jl:590-594` - Update `StateMomentProblem` struct to include block basis
2. `src/optimization/moment.jl:674-694` - Update `moment_relax` to store block bases with constraints
3. `src/optimization/sos.jl:252-330` - Rewrite `_sos_dualize_state` to iterate over blocks directly
4. `src/optimization/sos.jl:332-366` - Remove `_get_state_Cαj` (no longer needed)

## Testing Strategy

1. Run `make test` (full suite with Mosek)
2. Key tests to verify:
   - `test/problems/state_polynomial/` - state polynomial examples
   - `test/problems/quantum_networks/bilocal_networks.jl` - term sparsity test (expected -4.0)
3. Dense tests should continue passing (same logic, just single block)

## Constraints

- DO NOT modify expected test results
- Must maintain backward compatibility with dense (no term sparsity) case
- Must match NCTSSOS oracle values
