# Plan: Mermin Square GNS Reconstruction Example

## Objective
Create a literate example in `docs/src/examples/literate/` that demonstrates GNS reconstruction for the Mermin Square game, extracting finite-dimensional matrix representations of the quantum operators A[i,j] and B[i,j] from the moment matrix computed by the SDP relaxation.

## Background: Mermin Square Game
The Mermin Square is a 3×3 grid of observables:
- Row operators: A[i,1], A[i,2], A[i,3] for i=1,2,3
- Column operators: B[1,j], B[2,j], B[3,j] for j=1,2,3

**Game constraints:**
- Each row product equals +1: A[i,1] * A[i,2] * A[i,3] = +1
- Each column product equals -1: B[1,j] * B[2,j] * B[3,j] = -1
- Operators anti-commute across different positions
- Compatibility constraint: A[i,j] = B[i,j] (operators at the same grid position are identical)

This game demonstrates quantum contextuality - it's classically impossible but quantum mechanically achievable.

## Proposed Approach

### Key Data Structures
1. **Polynomial variables**: Complex non-commuting polynomial variables A[1:3,1:3], B[1:3,1:3]
2. **Constraints**:
   - `game_cons`: Product constraints for rows (A) and columns (B)
   - `comm_cons`: Commutativity constraints between operators
   - `entry_cons`: Compatibility constraints A[i,j] = B[i,j]
3. **Optimization problem**: Complex polynomial optimization with unipotency
4. **Moment matrix**: Extracted from the dual SDP solution

### Core Logic
1. Set up the Mermin Square polynomial optimization problem with all constraints
2. Use custom `cs_nctssos_with_entry` function to inject entry constraints into the moment relaxation
3. Solve the SDP relaxation to obtain a moment matrix (Hankel matrix)
4. Extract the moment matrix from the JuMP model constraints
5. Use `reconstruct()` to obtain matrix representations of all 18 operators
6. Verify the reconstructed matrices satisfy:
   - Game constraints (row/column products)
   - Anti-commutation relations
   - Compatibility (A[i,j] ≈ B[i,j])
   - Operator properties (squares, Hermiticity if applicable)

## Implementation Plan

### Step 1: File Structure and Introduction
Create `docs/src/examples/literate/mermin_square_gns.jl` with:
- Title and introduction to the Mermin Square game
- Mathematical description of the game constraints
- Explanation of quantum vs classical strategies

### Step 2: Package Imports (Clean Version)
```julia
using NCTSSoS
using NCTSSoS.FastPolynomials
using JuMP
using SCS
```
Remove all `Pkg.add` commands from the provided code.

### Step 3: Helper Function
Include the `cs_nctssos_with_entry` function definition (provided by user) that:
- Takes an optimization problem and entry constraints
- Builds the moment relaxation with correlative sparsity
- Injects entry constraints as HPSD constraints
- Returns a PolyOptResult with the solved model

### Step 4: Problem Setup
- Define dimensions (n=3, T1=ComplexF64)
- Declare non-commuting polynomial variables A[1:3,1:3], B[1:3,1:3]
- Set objective (trivial: one(T1) * one(A[3,3]))
- Define game_cons, comm_cons, entry_cons as provided
- Create the complex polynomial optimization problem with `is_unipotent=true`

### Step 5: Solve SDP Relaxation
- Configure solver (SCS with appropriate settings)
- Create SolverConfig with order=2
- Call `cs_nctssos_with_entry(pop, solver_config, entry_cons; dualize=true)`
- Extract the model from the result

### Step 6: Extract Hankel Matrix
From the user's code snippet:
```julia
cons = all_constraints(model, include_variable_in_set_constraints=true)
H = value.(cons[1])  # First constraint should be the moment matrix
```
Need to:
- Identify which constraint contains the Hankel matrix
- Extract its value
- Document the structure/indexing

### Step 7: Build Basis and Call reconstruct()
- Determine the basis: use `get_basis(variables, degree)` where degree=2 (from solver_config)
- Identify all unique variables for reconstruction (flatten A and B arrays)
- Call `reconstruct(H, vars, H_deg; atol=1e-6)` with appropriate degree

### Step 8: Reshape and Verify Results
- Reshape the vector of reconstructed matrices back into A_recon[1:3,1:3] and B_recon[1:3,1:3]
- Verify algebraic properties:
  - Row products: `A_recon[i,1] * A_recon[i,2] * A_recon[i,3] ≈ I`
  - Column products: `B_recon[1,j] * B_recon[2,j] * B_recon[3,j] ≈ -I`
  - Compatibility: `A_recon[i,j] ≈ B_recon[i,j]`
  - Squares: `A_recon[i,j]^2 ≈ I` (if they're involutions)
- Print verification results with tolerance checks

### Step 9: Documentation and Analysis
Add markdown comments explaining:
- The dimension of the reconstructed Hilbert space (size of matrices)
- Physical interpretation of the results
- Connection to quantum contextuality
- Any interesting observations about the matrix structure

## Open Questions / Decisions Needed

1. **Hankel matrix extraction**: The user's code shows `value.(cons[1])` - need to verify this is indeed the moment matrix. May need to filter constraints by type or inspect the structure.

2. **Basis construction**: The degree is 2 (from solver_config), but need to confirm:
   - Should we use all variables [A[:]..., B[:]...] or exploit A[i,j] = B[i,j]?
   - What's the total degree for the basis?

3. **Variable ordering**: When calling `reconstruct()`, need consistent ordering between:
   - The basis used to construct H in the solver
   - The variables passed to reconstruct()

4. **Tolerance values**:
   - `atol` for reconstruct (singular value cutoff)
   - Verification tolerances for algebraic relations

   May need to experiment to find appropriate values.

5. **SCS solver settings**: Should we customize SCS parameters for better accuracy?

## Expected Outcome

A complete literate example file that:
- Explains the Mermin Square game clearly
- Demonstrates the full workflow: problem setup → SDP solving → GNS reconstruction → verification
- Produces finite-dimensional matrix representations (likely 4×4 or 8×8 based on quantum realization)
- Verifies all game constraints are satisfied by the reconstructed operators
- Serves as a template for other quantum game/contextuality examples
