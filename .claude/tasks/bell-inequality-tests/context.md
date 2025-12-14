# Task: bell-inequality-tests

## Request
Bring back `test/bell_ineq.jl` using the new FastPolynomials API. The original test file was removed during the API redesign and needs to be updated to work with the new type-parameterized algebra system.

## Project State
- Branch: `feat/bell-inequality-tests` (based on `fastpolynomial-redesign`)
- Original file preserved in: `.claude/tasks/bell-inequality-tests/original_bell_ineq.jl`
- The project uses a two-layer architecture:
  1. FastPolynomials: Type-parameterized algebra system
  2. NCTSSoS: Optimization framework

## Key Changes from Old to New API
- Old: `@ncpolyvar X[1:5] Y[1:5]` macro for variable creation
- New: Type-parameterized `create_*_variables()` functions with `VariableRegistry`
- Old: `polyopt(p, comm_gps=[...], is_projective=true)`
- New: `polyopt(objective, registry)` with new constraint/configuration system

## Decisions
- [ ] Determine what algebra type Bell inequality tests should use
- [ ] Determine how to handle `comm_gps` (commuting groups) in new API
- [ ] Determine how to handle `is_projective` flag in new API

## Progress
- [x] Original file retrieved from git history
- [ ] Research complete
- [ ] Implementation plan created
- [ ] Implementation
- [ ] QA Review

## Handoff Summary
Original file saved. Needs research into:
1. What the Bell inequality tests are mathematically testing
2. How `comm_gps` and `is_projective` map to the new API
3. Current API patterns for similar tests (heisenberg.jl, etc.)
