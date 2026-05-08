# H4/Nk sector-pruning diagnostic

- generated: 2026-05-07T16:57:05.068
- input: `output/h4_pruning`
- red `n_unexpected` means target-sector orphan keys survived the symmetry classifier; that is the excess-kernel proxy.

| Nk | blocking | target sector | n_moments | n_orphans | n_expected_zero | n_unexpected | status_unpruned | bpsdp_iters_unpruned | walltime_unpruned | status_pruned | bpsdp_iters_pruned | walltime_pruned | obj_diff | artifact |
|----|----------|---------------|-----------|-----------|-----------------|--------------|-----------------|----------------------|-------------------|---------------|--------------------|-----------------|----------|----------|
| 2 | momentum | `N=0,K=0,2Sz=0,P=even` | 170273 | 46624 | 46624 | 0 | — | — | — | ITERATION_LIMIT | 200 | 248.332 | — | `h4_nk2` |

## Notes

- Nk=2: objective equivalence not checked (need both unpruned and pruned runs).
