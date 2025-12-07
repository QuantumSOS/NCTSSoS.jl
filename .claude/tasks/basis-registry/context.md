# Task: basis-registry

## Request
Rewrite `basis.jl` to use `VariableRegistry` for type consistency and correct variable indices.

## Current Branch
`fastpoly-basis-registry` (branched from `fastpolynomial-redesign`)

## Project State
- FastPolynomials review in progress (8/18 files done)
- 25 testsets skipped due to Int/UInt64 type mismatch between basis.jl and VariableRegistry
- Plan created at `.claude/tasks/basis-registry/plan.md`

## Key Decisions
1. **New API:** `get_ncbasis(registry::VariableRegistry{A,T}, d)` - registry-aware
2. **VariableRegistry{A,T}:** Add algebra type A as type parameter (Option A chosen) ✓ VALIDATED
3. **Core algorithm:** ~~Generate ALL words → simplify each → collect unique canonical monomials~~ **REJECTED** ❌
   - **NEW**: Generate ONLY canonical words directly (canonical-direct generation)
4. **Fermionic/Bosonic:** ~~Don't pre-filter to normal order~~ **REJECTED** ❌
   - **NEW**: Generate words directly in normal order; avoid non-canonical words entirely

## Validation Results (2025-12-07)
**Status**: Plan requires major revision before implementation

**Findings**:
1. ❌ **Core Algorithm INCORRECT**: Generate-all-then-simplify violates NCSOS canonical-only requirement
   - Creates exponential redundancy
   - Mixes degree levels (breaks moment matrix structure)
   - Violates linear independence
2. ❌ **Degree Mixing**: Simplification of a₁a₁† → 1 - a₁†a₁ produces degree-0 + degree-2 terms
   - Both should NOT be included in degree-2 basis
   - Only canonical degree-2 term (a₁†a₁) belongs in degree-2 basis
3. ✓ **VariableRegistry{A,T}**: Validated as correct approach, no issues found

**Key Sources**:
- Wang & Magron (2021): NCTSSOS.jl theory - basis must be canonical-only
- Wittek (2015): Ncpol2sdpa - direct canonical generation, not simplify-based
- Local tests: Confirmed a₁a₁† → degree-0 + degree-2 split

## NEXT ACTION
**CRITICAL**: Revise plan.md to use canonical-direct generation algorithm before implementation.
See `/Users/yushengzhao/projects/NCTSSoS-main/.claude/tasks/basis-registry/validation.md` for detailed corrections.

## Handoff Summary
- **Completed**: Plan validation with lead-researcher
- **Key Finding**: Original algorithm fundamentally flawed - must generate canonical words directly, not via simplification
- **Decision Made**: VariableRegistry{A,T} approved; core algorithm needs complete redesign
- **Blocker**: Cannot proceed with implementation until plan.md is corrected
- **Next Step**: User must review validation.md and approve revised algorithm approach
