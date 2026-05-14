# H4/Nk=2 solver benchmark

Lowering: `formulation=:psd_blocks`, `representation=:complex`, `orphan_policy=:aux_psd_free`.

## Problem

| quantity | value |
|---|---:|
| `spin_orbital_modes` | 32 |
| `total_electrons` | 8 |
| `hamiltonian_monomials` | 23752 |
| `moment_eq_polynomials` | 3 |
| `total_canonical_moments` | 123649 |
| `direct_real_moment_variables` | 247298 |
| `real_lift_psd_scalar_rows` | 1545187 |
| HPSD block sizes | `[240, 256, 240, 256, 513, 512]` |

## Results

| solver | status | raw status | lower s | solve s | iter | objective | primal res | dual res | gap | error |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `scs_indirect` | `ITERATION_LIMIT` | `solved (inaccurate - reached max_iters)` | 7.583091005 | 324.918413851 | 500 | -3.7447951801265003 |  |  | 0.003292024433422025 |  |
