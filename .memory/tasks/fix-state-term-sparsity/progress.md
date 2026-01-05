# Progress Log: Fix State Polynomial Term Sparsity

## 2026-01-04 - Ralph Loop Iteration 1

**Agent:** Orchestrator (Opus 4.5)

**Actions:**
1. Explored test structure and understood ralph-loop context
2. Read term_sparsity_debug_progress.md - identified prior work on this issue
3. Explored SOS dualization code for state polynomials
4. Identified root cause: `_get_state_Cαj` searches full basis but uses block-local (row,col) indices
5. Designed fix: store block basis with constraints, iterate directly over blocks
6. Created plan.md with detailed implementation steps

**Outcome:** Plan completed and documented

**Next step:** Implement the 4 code changes:
1. Update `StateMomentProblem` struct
2. Update `moment_relax` to store block bases
3. Rewrite `_sos_dualize_state` to iterate over blocks directly
4. Remove `_get_state_Cαj`
