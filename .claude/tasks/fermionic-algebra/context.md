# Fermionic Algebra Implementation Task

## Task Overview
Implement fermionic algebra for the NCTSSoS.jl package, following the pattern established by the Pauli algebra implementation.

## Current Status
- **Phase**: COMPLETE - All phases finished, tested, and committed
- **Date Created**: 2025-11-05
- **Last Updated**: 2025-11-13 (Hubbard model tests added)
- **Ready for**: Pull Request to main branch

## Latest Update (2025-11-13 Evening - Second Session)
- ✅ Computed exact baseline ground state energies via direct exact diagonalization
- ✅ Created `scripts/compute_exact_hubbard.jl` for baseline computation
- ✅ Updated testing plan with actual computed values
- **Exact Values**:
  - U=0: E₀ = -1.4142 (= -√2, free fermions)
  - U=1: E₀ = -1.0000
  - U=4: E₀ = -1.0000
  - U=8: E₀ = -1.0000
- **Status**: Baseline values computed, ready for commit

## Completed Actions
1. ✅ Created GitHub issue #179 for FastPolynomials zero monomial support
2. ✅ Comprehensive analysis of Pauli algebra implementation (see `pauli_analysis.md`)
3. ✅ Initial implementation plan created (see `plan.md` - archived in git history)
4. ✅ Plan revision completed (2025-11-12) - Constraint-based approach
5. ✅ **Phase 0 Complete (2025-11-12)**: AbstractAlgebra type hierarchy implemented
   - Defined `AbstractAlgebra` abstract type
   - Created `PauliAlgebra` and `FermionicAlgebra` structs
   - Updated `pauli_algebra()` to return `PauliAlgebra` instance
   - Added comprehensive tests (18 new tests, all 477 tests passing)
   - Fixed module include order for proper type availability
6. ✅ **Phase 1 Complete (2025-11-12)**: fermionic_algebra constructor implemented
   - Implemented `fermionic_algebra(N)` function in `src/algebra_constructors.jl`
   - All fermionic algebra properties encoded as polynomial equality constraints
   - Exported `fermionic_algebra` from NCTSSoS.jl
   - Added comprehensive test suite (47 fermionic algebra tests)
   - All tests passing (basic structure, constraint verification, scaling, cpolyopt integration, custom constraints, error handling)
7. ✅ **Phase 4 Complete (2025-11-12)**: Documentation tutorial created
   - Created `docs/src/examples/literate/fermionic_algebra_interface.jl` tutorial
   - Follows same pattern as `pauli_algebra_interface.jl`
   - Covers: manual vs simplified approach, free fermion example, adding constraints, implementation notes
   - Added to documentation index in `docs/make.jl`
   - Generated markdown file at `docs/src/examples/generated/fermionic_algebra_interface.md`
8. ✅ **Final Commit (2025-11-13)**: All work committed and tested
   - Commit: `614fef2` - "feat: add fermionic_algebra constructor with constraint-based implementation"
   - Fixed test failure: Removed Literate from [deps] in Project.toml (Aqua.jl stale dependency issue)
   - All 524 tests passing (477 existing + 47 new fermionic algebra tests)
   - Session documented in `.claude/history/2025-11-13_fermionic-algebra-complete.md`

## User-Requested Simplifications (2025-11-12)

**Key Changes**:
1. **Remove zero monomial/nilpotency support**: Do NOT modify FastPolynomials, PolyOpt, or ComplexPolyOpt
2. **Direct polynomial constraint approach**: Implement c² = 0 and c†² = 0 as polynomial constraints (not via simplification)
3. **Create algebra struct**: Define a struct for fermionic algebra (similar to Pauli algebra if it has one)
4. **Follow pauli_algebra pattern**: Directly implement fermionic_algebra interface like pauli_algebra

## Key Requirements

### Fermionic Algebra Implementation
**Pattern to Follow**: Imitate the Pauli algebra implementation structure exactly

**Fermionic Algebra Properties**:
- Creation operators: c†
- Annihilation operators: c
- Anti-commutation relations: {c_i, c_j} = 0, {c†_i, c†_j} = 0, {c_i, c†_j} = δ_ij
- Nilpotent property (as constraints): c² = 0, (c†)² = 0

**Implementation Approach**:
- Create a struct for fermionic algebra
- Generate polynomial constraints for all anti-commutation relations
- Add polynomial constraints for c² = 0 and c†² = 0
- No modifications to FastPolynomials or optimization types needed

## Dependencies
1. Understanding of Pauli algebra implementation structure
2. Understanding of fermionic operator algebra mathematics

## Files to Explore
- Pauli algebra implementation files (for struct and interface pattern)
- Existing algebra constructors
- Test structure for algebra implementations

## Revision Instructions for Sub-Agent

**Task**: Revise plan.md with simplified approach

**Starting Point**: Create a struct for fermionic algebra first

**Requirements**:
1. Analyze Pauli algebra struct (if exists) and interface
2. Design fermionic algebra struct
3. Implement fermionic_algebra(N) constructor function
4. Generate constraints for:
   - Anti-commutation relations: {c_i, c_j} = 0, {c†_i, c†_j} = 0, {c_i, c†_j} = δ_ij
   - Nilpotent relations: c_i² = 0, (c†_i)² = 0
5. Return appropriate type compatible with cpolyopt
6. Write tests following existing patterns
7. Write documentation following pauli_algebra style

**What NOT to do**:
- Do NOT modify FastPolynomials
- Do NOT modify PolyOpt or ComplexPolyOpt
- Do NOT add nilpotent/zero monomial support to existing types

**Next Steps**:
1. Sub-agent revises plan.md
2. User reviews and approves revised plan
3. Parent agent implements via TDD

## Plan Revision Summary (2025-11-12)

### Architectural Decision: Constraint-Based Approach

The revised plan adopts a **pure constraint-based approach** that eliminates all core type modifications:

**Key Changes from Original Plan**:
1. **No SimplifyAlgorithm modifications**: No `is_nilpotent` flag, no `_simplify_nilpotent!` function
2. **No PolyOpt modifications**: No new fields in PolyOpt or ComplexPolyOpt
3. **Nilpotency as constraints**: c_i² = 0 and (c†_i)² = 0 added as polynomial equality constraints
4. **Follows pauli_algebra pattern exactly**: Returns NamedTuple with standard fields
5. **Dramatically simpler**: ~40 lines of new code vs ~500 lines of modifications

**How It Works**:
- All fermionic algebra properties encoded as polynomial equality constraints
- Anti-commutation relations: {c_i, c_j} = 0, {c†_i, c†_j} = 0, {c_i, c†_j} = δ_ij
- Nilpotent constraints: c_i² = 0, (c†_i)² = 0 (redundant but explicit)
- SDP solver enforces all constraints during optimization
- No need for zero monomial support or workarounds

**Constraint Count**: 2N² + 3N for N modes
- N=1: 5 constraints
- N=2: 14 constraints
- N=3: 27 constraints
- Recommended N ≤ 10 for practical performance

**Benefits**:
- ✅ Much simpler implementation
- ✅ No core type modifications
- ✅ No workarounds needed
- ✅ Easier to maintain and extend
- ✅ Self-documenting (all algebra properties visible as constraints)

## Notes
- This implementation maintains complete consistency with pauli_algebra
- Documentation follows the same pattern as Pauli algebra
- Tests cover all fermionic algebra properties
- Constraint-based approach: No core type modifications required
- **Hubbard model tests**: Demonstrate real physics applications
- **Limitation identified**: Hopping terms require Jordan-Wigner transformation
- **Future work**: Exact verification using ExactDiagonalization.jl or XDiag.jl
