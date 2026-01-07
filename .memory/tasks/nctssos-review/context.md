# NCTSSoS.jl Review Context

## Goal
Comprehensive codebase familiarization + quality audit of NCTSSoS.jl, a Julia package for sparse noncommutative polynomial optimization using moment-SOHS hierarchy.

## Approach
- **Use-case driven**: Trace CHSH inequality through full solver pipeline
- **Audit focus**: Correctness, Maintainability, Test Coverage
- **Phase order**: Optimization pipeline first, then types/simplification as needed

## Key Decisions

| Decision | Rationale | Date |
|----------|-----------|------|
| Start with optimization | User preference: top-down understanding | 2026-01-07 |
| Skip Python comparison | Focus on Julia code quality only | 2026-01-07 |

## Architecture Summary
*(To be filled as review progresses)*

## Key Findings
*(To be filled during review)*

## Handoff Summary

**Last session**: Initial setup
**Completed**: Plan creation, session tracking structure
**Next step**: Begin Phase 1.1 - Review `src/optimization/problem.jl`
