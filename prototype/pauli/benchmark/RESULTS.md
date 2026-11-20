# Recorded periodic Heisenberg results

The order-2 radius-4 candidate reached 16×16 with three fresh numerical passes. The energy at that size is **2.360% below the published QMC estimate**. This is a size result, not a reproduction of Wang's Table 10 settings or a rigorous bound on the error. All energies are per spin.

| Square | Accepted | SDP objective | QMC | Gap below QMC | Median fresh wall (s) | Peak worker RSS (GiB) |
|---|---:|---:|---:|---:|---:|---:|
| 8×8 | 3/3 | −0.6869715602 | −0.673487 | 2.002% | 80.65 | 3.34 |
| 10×10 | 3/3 | −0.6863629325 | −0.671549 | 2.206% | 174.49 | 5.99 |
| 12×12 | 3/3 | −0.6860647689 | −0.670685 | 2.293% | 427.72 | 13.29 |
| 14×14 | 3/3 | −0.6858875905 | −0.670222 | 2.337% | 907.44 | 31.67 |
| 16×16 | 3/3 | −0.6857884978 | −0.669976 | 2.360% | 2207.27 | 62.02 |

QMC values are from [Wang et al., Table 10](https://arxiv.org/html/2604.01555v2#S6.SS2). Gap = (QMC − SDP) / |QMC|. One thread per worker, Julia 1.12.7, MOSEK 11.2.8, 100 GiB memory limit, no swap. The 8×8 result used the SOS submission; the 10×10 SOS form stalled three times and was not accepted. An equivalent unscaled explicit moment submission passed at 10×10 and above. The later detector recovers both translations without hints; at the 10×10 overlap the two layouts agree under reflection, not byte-for-byte. Do not pool timings across these variants.

The accepted 10×10–16×16 runs were all `OPTIMAL` with primal and dual `FEASIBLE_POINT`. **None** of those twelve runs met the separate diagnostic requiring every reconstructed residual to be ≤1e−8: the dual stationarity maxima by size were 1.75e−8, 1.39e−8, 2.63e−8 and 2.57e−8. At 16×16 the first 1800 s attempt timed out without an objective. At 3600 s per worker, fresh wall times were 2198.56–2215.51 s and the cgroup peaked at 61.82 GiB. Stock and candidate row pools differ; there is no square speedup claim.

## Chain certificates are a separate result

On the periodic 200-spin chain, the independently replayed six-word energy-window certificate gave **−0.4432434483609903/spin**. Its error against the finite-ring Bethe reference (used for comparison only) was 7.56821e−5/spin, above the 1e−5 target. This was a separate research replay, not an output of the square-lattice benchmark or a general certification claim for this prototype.

## Identify the replay archives

The workspace backup retains raw records and inputs; only results and content hashes are published here. These hashes identify the original archives, not their expanded contents.

| Record | SHA-256 |
|---|---|
| 8×8 and rejected 10×10 SOS ladder | `2eb36b2c4ece16d1c946c2823bb8c694cbcf35d5d402a8034151aa7e2b9952bc` |
| 10×10–14×14 passes and first 16×16 timeout | `954d86c1c2d80ce6df09dee121fc06ff4cb0bb78e80dfe8178f73a30212fea87` |
| Three 16×16 explicit-moment passes | `ae527eb7094d1e0e21c6dc6570923ef609ceb10fdce3c74c9987ad3e17076186` |
| 200-spin six-word replay | `cb2ebadf705f9cbf202587347272e2b092fc7139ef6460cab9e974230085ecbc` |
| 10×10 pair-box audit inputs (no solve) | `0a5282eb7043f44fe55407dd1d84719b39d822dbf3ad43bbbbd2cc2f992d6479` |
