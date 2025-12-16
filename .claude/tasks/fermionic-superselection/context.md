# Task: Fermionic Superselection Rules Implementation

## Request
Research and implement parity superselection constraints for fermionic moment relaxations in polynomial optimization. The user believes that moment matrix entries corresponding to operators with odd parity (odd number of creation/annihilation operators) must be set to zero.

## Project State
- NCTSSoS.jl has a working fermionic algebra implementation with Wick's theorem simplification
- `FermionicAlgebra` uses complex coefficients (Float64 coefficients, not ComplexF64)
- Moment relaxation is handled by `moment_relax()` in `src/moment_solver.jl`
- Fermionic problems use Hermitian PSD (HPSD) cones
- Test file `test/fermionic_chain.jl` was already failing before this implementation (pre-existing issue)

## Decisions
- [x] ~~Implement parity filtering in moment basis~~ **INCORRECT APPROACH - ABANDONED**
- [x] **CORRECT APPROACH**: Keep ALL monomials in basis, add zero constraints for odd-parity moment matrix entries
- [x] Add objective validation to reject odd-parity objectives (KEEP THIS)
- [x] Modify moment_relax to add parity constraints instead of filtering basis
- [ ] Spinful fermions - no special handling needed (count all operators)

## Progress
- [x] Initial research on parity superselection theory
- [x] Current implementation analysis
- [x] **USER CORRECTION**: Understand why filtering is wrong
- [x] **REVISED RESEARCH**: Complete analysis of correct constraint-based approach
- [ ] Remove incorrect basis filtering
- [ ] Implement parity constraint injection
- [ ] Update tests to reflect new behavior
- [ ] Verify fermionic_chain.jl passes

## Why Previous Implementation Was Wrong

The previous implementation **filtered odd-parity monomials from the basis**, which is INCORRECT because:

1. **Moment matrix entry** `M[i,j] = ⟨basis[i]† · operator · basis[j]⟩`
2. **Even if** `basis[i]` and `basis[j]` have odd parity, the **total operator** can have even parity!
3. **Example**: `basis[i] = a[1]` (odd), `operator = I` (even), `basis[j] = a[1]` (odd)
   - Total: `a_dag[1] · I · a[1] = a_dag[1] * a[1]` (even parity!)
   - This represents the number operator and should be allowed
4. **Filtering removes needed basis elements** → primal infeasible

### Test Evidence
- `fermionic_chain.jl` fails with "primal infeasible" due to missing basis elements
- Ground state energy: -2.48e10 (nonsense) instead of -4.0 (correct)

## Files That Need Modification

### TO BE REMOVED (Incorrect Filtering)
1. `src/moment_solver.jl:152-156` - Remove basis filtering

### TO BE KEPT (Still Correct)
1. `has_even_parity()` helper - KEEP (still needed for checking)
2. `src/pop.jl` objective validation - KEEP (objective must have even parity)

### TO BE ADDED (Correct Implementation)
1. `src/moment_solver.jl` - Add `_has_odd_parity_only(poly)` helper
2. `src/moment_solver.jl` - Add `_add_parity_constraints!(mp)` function
3. `src/moment_solver.jl` - Call constraint injection at end of `moment_relax`

### TO BE UPDATED (Tests)
1. `test/fermionic_parity_test.jl` - Update to verify constraints instead of filtering
2. Add test verifying odd-parity monomials ARE in basis
3. Add test verifying zero constraints ARE added

## Handoff Summary
- **Research Completed**: Identified why current filtering approach is WRONG
- **Key Finding**: Basis filtering removes needed monomials → primal infeasible
- **Correct Approach**: Keep ALL basis monomials, add zero constraints for odd-parity moment entries
- **Test Failure Root Cause**: fermionic_chain.jl fails BECAUSE of incorrect basis filtering (not pre-existing)
- **Implementation Plan**: Documented in research.md with step-by-step code changes
- **Ready for Engineer**: Complete specification for fixing the implementation
- **Expected Outcome**: After fix, fermionic_chain.jl should pass with E₀ ≈ -4.0
