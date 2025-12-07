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
2. **VariableRegistry{A,T}:** Add algebra type A as type parameter (Option A chosen)
3. **Core algorithm:** Generate ALL words → simplify each → collect unique canonical monomials
4. **Fermionic/Bosonic:** Don't pre-filter to normal order; simplify captures all canonical forms

## NEXT ACTION (after context clear)
**Use `lead-researcher` agent to validate this plan before implementation.**

Questions to validate:
1. Is the core algorithm (generate all → simplify → dedupe) correct for NCTSSOS-style optimization?
2. Are there edge cases in fermionic/bosonic simplification we're missing?
3. Does adding algebra type to VariableRegistry have any downstream implications?

## Handoff Summary
- **Completed**: Plan created and reviewed with user
- **Key Finding**: Fermionic basis needs all words (not just normal-ordered) because non-normal words simplify to SUMS, not zero
- **Decision Made**: VariableRegistry gets algebra type parameter {A,T}
- **Next Step**: Validate plan with lead-researcher, then implement Step 1
