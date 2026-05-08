# Pivot-coverage audit

Generated: 2026-05-07 22:01:42
Git commit: `unknown`

Purpose: measure whether cached `MomentLinearData.pivots` covers every canonical moment key, or whether `MomentLinearData.free_keys` must carry orphan moments.

| case | build s | moments | pivots | orphans | orphan % | category | PSD/HPSD blocks | max block | zero constraints |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| `h2_chain_nk2` | 11.032 | 2913 | 2913 | 0 | 0.0000 | No orphans | 18 | 65 | 6 |

## Details

### `h2_chain_nk2`

primary H2/Nk=2 PQG BPSDP target; no explicit N̂-N constraint; spin blocks + paper spin constraints

- constraints: 24
- orphan category: No orphans
- orphan reason counts:
  - none

## Gate decision
All audited cases have no orphans. Production lowering can treat `free_keys` as empty for these workloads.
