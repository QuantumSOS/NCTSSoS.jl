# Bell memory growth is after JuMP lowering

No lowering code was changed. The stage measurements do not support attributing the 20 GB failure to `Next.to_jump`. The n=4 profile places the large native work inside `MosekTools.optimize! → Mosek.optimize → MSK_optimizetrm`. At n=5, lowering and optimizer attachment both return before the process crosses the memory cap during `JuMP.optimize!`.

## Measured stages

Dense I_nn22, order 3, with the same `perf/next/probe.jl` inputs as the scale probe. Julia 1.12.7, JuMP 1.32.1, MosekTools 0.15.11, Mosek 11.2.1, Linux x86_64, shared Intel Xeon orb. Julia and BLAS each use one thread; n=4 Mosek uses one thread, n=5 uses 16. Base d97ca7f with the native-writer work; `to_jump` and the JuMP write path are unchanged from that base. Bytes are decimal. Julia allocations do not include Mosek's native heap.

| Stage | n=4 time s | n=4 Julia bytes | n=5 time s | n=5 Julia bytes |
|---|---:|---:|---:|---:|
| `Next.to_jump` | 0.02154 [0.01884,0.03029] | 40,238,296 | 0.81716 | 171,395,592 |
| Lowering plus optimizer attachment | 0.03338 [0.03195,0.03822] | 48,388,088 | — | — |
| Attachment alone | — | — | 3.44060 | 426,847,736 |

n=4 entries are medians and ranges of three warmed calls. The two n=5 calls are single cold stage measurements, not a timing comparison with n=4. Both include a model with the complete dense PSD constraint; no symmetry or sparsity relaxation was substituted. The standard matrices have respectively (23,653, 6,156) and (95,266, 27,555) rows/columns, with 23,652 and 95,265 nonzeros. Their PSD sides are 217 and 436.

The n=4 process peak was 1,251,102,720 bytes before `JuMP.optimize!` and 12,475,506,688 bytes after it. The profiled optimize call took 125.415 s, allocated 59,362,392 Julia bytes, and returned `TIME_LIMIT` at its 120-second solver limit. Its time is diagnostic, not a solve-time benchmark: it includes compilation and sampling overhead. The native log reports 280 million factor nonzeros and 4.41×10¹² factorization flops. Samples include `ipfac_symbolicsetup`, `ipfac_setupfactor`, and native BLAS kernels under `MSK_optimizetrm`. See [the full stage and sampling log](mosek_lowering.txt).

The n=5 watchdog killed the process at 20,037,947,392 bytes after 52.182 s. Its last completed stages were:

```text
BUILD 1.38459842 240506072 SIZES 436 rows, 27556 moments, 27555 real variables, 1 blocks (largest 436)
FORM (95266, 27555) NNZ 95265
WRITE 0.817161008 171395592 ATTACH 3.440595504 426847736
PROCESS bell jump elapsed=52.182 peak_bytes=20037947392 returncode=-9 reason=memory_limit
```

The previous scale probe recorded no optimizer-startup text. That was not evidence that optimization had not started: the new flushed stage markers locate the failure after attachment, during `JuMP.optimize!`. The n=4 C-stack profile identifies the native solver work in that call. We did not profile the n=5 C allocator, so these results do not identify an individual native allocation responsible for its peak.

## Reproduction

Use the disposable solver environment from the writer tests, with MosekTools and a verified license. Run numerical jobs serially.

```sh
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /usr/bin/time -v julia --project=PATH_TO_ENV --startup-file=no perf/next/mosek_lowering.jl
```

The n=5 stage case is `perf/next/mosek_acceptance.jl bell jump`, with `MOSEK_THREADS=16`, `MOSEK_TOLERANCE=1e-8` and `MOSEK_MAX_TIME=120`. It requires an external RSS watchdog: the solver time limit alone does not prevent exhausting memory during setup. The measurement checked `/proc/PID/status` every 100 ms and killed at 20,000,000,000 bytes; the observed peak can overshoot between samples.

The n=4 profile and n=5 stage run took 154.93 + 52.18 = 207.11 s, within the 3–6 minute estimate. This finding changes the attribution, not the requested relaxation. The native SOS writer avoids the moment-form setup in this comparison; its Bell bound and memory are reported separately.
