# Term Sparsity Debug Progress

## Problem
The term sparsity implementation for NCStateWord (state polynomial optimization) in NCTSSoS is not working correctly. The bilocal network test (`test/problems/quantum_networks/bilocal_networks.jl`) fails with:
- Expected objective: -4.0
- Actual objective: -1.76M (completely wrong)
- Expected blocks: 610 with sizes [16×2, 5×12, 3×24, 1×572]
- Actual blocks: 552 with sizes [9, 7×13, 2×many, 1×many]

## Root Cause Analysis (Iteration 1)

### Issue 1: FIXED - Over-inclusive activated support
The original `init_activated_supp` for NCStateWord included ALL pairwise products `_neat_dot3(bi, I, bj)`, which made the term sparsity graph fully connected and destroyed all sparsity.

**Fix applied** (src/optimization/sparsity.jl lines 815-835):
- Changed to only include diagonal entries `bi†*I*bi` (monosquare terms)
- Plus objective monomials and constraint monomials

### Issue 2: FIXED - Missing canonicalization in graph edge check
The graph edge check was comparing raw NCStateWords instead of canonicalized forms.

**Fix applied** (src/optimization/sparsity.jl lines 873-910):
- Added `symmetric_canon(ncsw)` before checking membership in activated support

### Issue 3: STILL PRESENT - Graph structure mismatch
After fixes, the graph has:
- 649 vertices, 321 edges
- 551 connected components (max size 15)
- MMD produces 552 cliques (max size 9)

But NCTSSOS produces:
- 610 blocks (max size 16)
- Different clique structure

### Key Observation
The objective is still completely wrong (-1.76M), even though:
1. The graph has correct edges (verified: edge between identity and ⟨x₁y₁z₁⟩ exists)
2. Linear objective terms ARE in the activated support
3. The term sparsity blocks are being generated

This suggests the moment relaxation built from these blocks is still incorrect.

## Hypotheses for Further Investigation

### Hypothesis A: Basis size mismatch
NCTSSoS has 649 basis elements, NCTSSOS has 637. This could affect the block structure.

### Hypothesis B: Block overlap handling
In NCTSSOS, blocks can overlap (sum of block sizes > basis size: 736 > 637).
In NCTSSoS, blocks appear non-overlapping (sum = 649).

### Hypothesis C: Different block decomposition algorithm
NCTSSOS uses `chordal_cliques!` which may produce different results than NCTSSoS's `clique_decomp`.

### Hypothesis D: Constraint matrix construction issue
The `_build_state_constraint_matrix` function might be constructing matrices incorrectly for the given blocks.

## Files Modified
- `src/optimization/sparsity.jl`:
  - `init_activated_supp` for NCStateWord (lines 815-835)
  - `get_term_sparsity_graph` for NCStateWord (lines 873-910)

## Iteration 1 - Additional Findings

### Verified Working
1. **Term sparsity graph construction**: Graph has correct edges (e.g., edge between ⟨xyz1⟩ and ⟨xyz2⟩ exists)
2. **Block decomposition**: MMD correctly puts connected basis elements in same block
3. **Constraint matrix construction**: `_build_state_constraint_matrix` produces correct entries (M[1,2] = ⟨xyz1⟩⟨xyz2⟩)

### Test Case Analysis
For simple cross-term objective `0.25 * ⟨xyz1⟩⟨xyz2⟩`:
- Dense objective: -0.25 (correct)
- TS objective: -8.84 (WRONG - too negative)
- Block 10 correctly contains [⟨xyz2⟩, ⟨xyz1⟩]
- Constraint matrix for block 10 is a 2×2 matrix with correct entries

### Remaining Issue
The SOS dualization or JuMP model construction must be incorrect. The constraint matrices are correct, but the optimization result is wrong.

Possible causes:
1. Incorrect coefficient extraction in `_get_state_Cαj`
2. Incorrect basis indexing when mapping NCStateWords to StateWords
3. Double-counting or missing constraints in the assembly

## Next Steps
1. Add debug output to `_sos_dualize_state` to trace coefficient extraction
2. Compare the final JuMP model between Dense and TS cases
3. Check if all StateWords in the objective are correctly matched to constraint coefficients
