# Fermionic Algebra Implementation Task

## Task Overview
Implement fermionic algebra for the NCTSSoS.jl package, following the pattern established by the Pauli algebra implementation.

## Current Status
- **Phase**: Planning complete - Awaiting user approval to proceed with implementation
- **Date Created**: 2025-11-05
- **Last Updated**: 2025-11-05

## Completed Actions
1. ✅ Created GitHub issue #179 for FastPolynomials zero monomial support
2. ✅ Comprehensive analysis of Pauli algebra implementation (see `pauli_analysis.md`)
3. ✅ Complete implementation plan created (see `plan.md`)

## Key Requirements

### 1. FastPolynomials Zero Monomial Support Issue
**Priority**: CRITICAL - Must be resolved before fermionic algebra implementation

**Problem**:
Current implementation of `FastPolynomials` doesn't support zero monomials. This is needed because fermionic algebra has the property that `c^2 = 0` (creation/annihilation operators square to zero).

**Actions Needed**:
- Create GitHub issue documenting this limitation
- Research how to extend FastPolynomials to support zero monomials
- Implement solution (or work around if necessary)

### 2. Fermionic Algebra Implementation
**Pattern to Follow**: Imitate the Pauli algebra implementation structure

**Fermionic Algebra Properties**:
- Creation operators: c†
- Annihilation operators: c
- Anti-commutation relations: {c_i, c_j} = 0, {c†_i, c†_j} = 0, {c_i, c†_j} = δ_ij
- Key property: c^2 = 0, (c†)^2 = 0

## Dependencies
1. Understanding of Pauli algebra implementation structure
2. Resolution of FastPolynomials zero monomial issue
3. Understanding of fermionic operator algebra mathematics

## Files to Explore
- Pauli algebra implementation files
- FastPolynomials implementation
- Test structure for algebra implementations

## Implementation Plan Summary

**Detailed plan location**: `.claude/tasks/fermionic-algebra/plan.md`

### Architecture Decisions

1. **Nilpotent Extension**: Add `is_nilpotent::Bool` to SimplifyAlgorithm, PolyOpt, ComplexPolyOpt
2. **Zero Monomial Workaround**: Use empty vars/z vectors to represent zero (temporary solution for issue #179)
3. **Variable Naming**: `c[i]` (annihilation), `c_dag[i]` (creation)
4. **Commutation Structure**: Single group containing all operators (all anti-commute)
5. **Constraint Count**: 2N² + N equality constraints for N modes

### Implementation Phases

**Phase 1**: Core Simplification Logic (FastPolynomials)
- Add `is_nilpotent` field to SimplifyAlgorithm
- Implement `_simplify_nilpotent!` function (detects X² and returns zero)
- Update simplify! dispatch to handle nilpotent case

**Phase 2**: Optimization Problem Types
- Extend PolyOpt and ComplexPolyOpt with `is_nilpotent` field
- Update constructors and algebra interface
- Enforce mutual exclusivity (only one of {unipotent, projective, nilpotent})

**Phase 3**: Fermionic Algebra Constructor
- Implement `fermionic_algebra(N)` following Pauli pattern
- Generate all anti-commutation relation constraints
- Return NamedTuple compatible with cpolyopt(obj, algebra)

**Phase 4**: Comprehensive Testing
- Unit tests: structure, constraints, nilpotent simplification
- Integration tests: with cpolyopt and cs_nctssos
- Numerical tests: free fermions, Fermi-Hubbard model (LOCAL_TESTING only)

**Phase 5**: Documentation
- Tutorial following Pauli algebra interface pattern
- API documentation with examples
- Cross-references and physical motivation

### Files to Modify/Create

**Modify**:
- `src/FastPolynomials/src/simplify.jl` (~50 lines changed + 20 new)
- `src/pop.jl` (~30 lines changed)
- `src/NCTSSoS.jl` (1 line: export)
- `src/interface.jl` (verify SimplifyAlgorithm construction)

**Create**:
- `src/algebra_constructors.jl` (append ~80 lines for fermionic_algebra)
- `test/algebra_constructors.jl` (append ~300 lines)
- `docs/src/examples/literate/fermionic_algebra_interface.jl` (new file ~200 lines)

**Total Impact**: ~500 lines new code, ~80 lines modifications

### Critical Implementation Notes

1. **TDD Required**: Write failing test first, implement minimum code, explain, refactor
2. **Zero Monomial Handling**: Empty monomial after nilpotent simplification represents zero (workaround)
3. **Constraint Scaling**: O(N²) constraints - document performance limitations for large N
4. **Mutual Exclusivity**: Enforce that only one of {is_unipotent, is_projective, is_nilpotent} is true
5. **Single Commutation Group**: All fermionic operators in one group (no reordering by simplifier)

### Next Steps
1. User reviews and approves plan.md
2. Parent agent implements via TDD following plan.md phases
3. Each step: write test → implement → explain → verify → next step

## Notes
- This implementation should maintain consistency with existing algebra implementations
- Documentation should follow the same pattern as Pauli algebra
- Tests should cover all fermionic algebra properties
