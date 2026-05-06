# Pivot-coverage audit

Generated: 2026-05-06 12:40:45
Git commit: `a7ad246`

Purpose: measure whether existing runtime pivot discovery covers every canonical moment key, or whether planned `MomentLinearData.free_keys` must carry orphan moments.

| case | build s | moments | pivots | orphans | orphan % | category | PSD/HPSD blocks | max block | zero constraints |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| `h2_chain_nk2` | 8.564 | 9393 | 2913 | 6480 | 68.9875 | Many orphans | 18 | 65 | 983 |
| `h4_chain_nk2_proxy_small` | 0.603 | 9393 | 2913 | 6480 | 68.9875 | Many orphans | 18 | 65 | 983 |
| `pauli_3qubit_ground_state` | 6.195 | 64 | 64 | 0 | 0.0000 | No orphans | 1 | 37 | 0 |
| `fermionic_4site_hubbard` | 5.206 | 2517 | 2517 | 0 | 0.0000 | No orphans | 1 | 137 | 7346 |
| `bosonic_2mode_truncated` | 6.852 | 70 | 70 | 0 | 0.0000 | No orphans | 1 | 15 | 27 |
| `cs_failure_E10` | 6.937 | 129 | 129 | 0 | 0.0000 | No orphans | 406 | 3 | 0 |

## Details

### `h2_chain_nk2`

primary H2/Nk=2 PQG BPSDP target; spin blocks + paper spin constraints

- constraints: 1001
- orphan category: Many orphans
- orphan reason counts:
  - zero constraints only: 6480
- sample orphan keys:
  - `Int8[-1, 1, 10, 11]`
  - `Int8[-1, 1, 10, 12]`
  - `Int8[-1, 1, 10, 13]`
  - `Int8[-1, 1, 10, 14]`
  - `Int8[-1, 1, 10, 15]`
  - `Int8[-1, 1, 10, 16]`
  - `Int8[-1, 1, 11, 12]`
  - `Int8[-1, 1, 11, 13]`
  - `Int8[-1, 1, 11, 14]`
  - `Int8[-1, 1, 11, 15]`
  - `Int8[-1, 1, 11, 16]`
  - `Int8[-1, 1, 12, 13]`

### `h4_chain_nk2_proxy_small`

2-k, 4-orbital half-filled PQG proxy for H4; structural audit only

- constraints: 1001
- orphan category: Many orphans
- orphan reason counts:
  - zero constraints only: 6480
- sample orphan keys:
  - `Int8[-1, 1, 10, 11]`
  - `Int8[-1, 1, 10, 12]`
  - `Int8[-1, 1, 10, 13]`
  - `Int8[-1, 1, 10, 14]`
  - `Int8[-1, 1, 10, 15]`
  - `Int8[-1, 1, 10, 16]`
  - `Int8[-1, 1, 11, 12]`
  - `Int8[-1, 1, 11, 13]`
  - `Int8[-1, 1, 11, 14]`
  - `Int8[-1, 1, 11, 15]`
  - `Int8[-1, 1, 11, 16]`
  - `Int8[-1, 1, 12, 13]`

### `pauli_3qubit_ground_state`

3-qubit periodic Heisenberg, order 2

- constraints: 1
- orphan category: No orphans
- orphan reason counts:
  - none

### `fermionic_4site_hubbard`

4-site periodic Hubbard, half-filled N_up=N_dn=2, order 2

- constraints: 7347
- orphan category: No orphans
- orphan reason counts:
  - none

### `bosonic_2mode_truncated`

2-mode Bose-Hubbard with N=2 moment-equality truncation, order 2

- constraints: 28
- orphan category: No orphans
- orphan reason counts:
  - none

### `cs_failure_E10`

deterministic 10-variable NC polyball stressor; historical E10 coefficients are unavailable

- constraints: 406
- orphan category: No orphans
- orphan reason counts:
  - none

## Gate decision
At least one audited case has many orphans. Keep explicit `free_keys`; do not pretend pivot coverage is universal.
