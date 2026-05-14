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
| `bpsdp` | `ITERATION_LIMIT` | `max_iter` | 8.047456652 | 528.243073212 | 500 | unavailable: Result index of attribute MathOptInterface.ObjectiveValue(1) out of bounds. There are currently 0 solution(s) in the model. | 1.733148405699826 | 2.6474347789010716e-5 | 0.0015530719813927618 |  |
