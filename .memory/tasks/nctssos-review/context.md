# NCTSSoS.jl Review Context

## Goal
Comprehensive codebase familiarization + quality audit of NCTSSoS.jl, a Julia package for sparse noncommutative polynomial optimization using moment-SOHS hierarchy.

## Approach
- **Interactive review**: I present files, you ask questions/direct focus
- **Collaborative**: I explain context + bigger picture before each file
- **Use-case driven**: Trace CHSH inequality through full solver pipeline
- **Audit focus**: Correctness, Maintainability, Test Coverage
- **Phase order**: Optimization pipeline first, then types/simplification as needed

## Review Protocol
1. **I present** the next file: what it does, why, bigger picture
2. **You decide**: dive in, skip, or redirect focus
3. **We examine together**: you guide what you want to understand
4. **I summarize findings** for your approval before moving on

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
