# Scale measurements

Measured on library revision [7891a73](https://github.com/QuantumSOS/NCTSSoS.jl/commit/7891a73399fb9fabb096d7e3fb98656928db0992), before continuous-symmetry PR #424. These cases use no continuous generators. The follow-up model/display tests use [9b5a127](https://github.com/QuantumSOS/NCTSSoS.jl/commit/9b5a127). No optimization is included.

Julia 1.12.7, JuMP 1.32.1, SCS 2.7.1, Mosek 11.2.1, MosekTools 0.15.11, COSMO 0.8.11. Linux x86_64, Intel Xeon 2.60 GHz, 16 visible CPUs in a shared orb; no CPU pinning. One Julia/BLAS/OpenMP/Mosek thread, serial workers. [Raw samples and full package versions](results.json); [sampling profiles](profiles.txt); [method and relaxation differences](README.md).

Build times are median [minimum, maximum] of three fresh warm builds, in seconds. Allocation and RSS units are decimal GB. RSS is the whole-process maximum, including compilation and lowering. Lowering reports median seconds and allocated GB; sizes are rows/columns/nonzeros for Next, scalar JuMP constraint dimensions/variables/affine coefficients for legacy. These are different representations. `—` means unmeasured, not zero. The 600 s / 20 GB process cap includes all warmups, repeats, and lowering.

| Case | Build s [min,max] | Build GB | Peak GB | Rows | Moments | Real vars | Blocks | Max side | Lower s / GB | Lower sizes | Completion |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| chain-next-16-4-contiguous | 0.223278 [0.209953,0.224369] | 0.0963491 | 2.08111 | 1921 | 1035 | 1034 | 54 | 16 | 0.00175114 / 0.00445257 | 2683/1034/54971 | complete |
| chain-legacy-16-4-contiguous | 0.171023 [0.165675,0.193191] | 0.201842 | 2.22528 | 541 | 1035 | 1035 | 20 | 30 | 0.0426613 / 0.0349112 | 7983/1035/85412 | complete |
| chain-next-32-4-contiguous | 0.779992 [0.726042,0.850309] | 0.305521 | 2.15937 | 3841 | 2755 | 2754 | 94 | 16 | 0.00747792 / 0.0243837 | 5299/2754/238463 | complete |
| chain-legacy-32-4-contiguous | 0.701042 [0.610623,0.730002] | 0.641147 | 2.35444 | 1021 | 2755 | 2755 | 36 | 30 | 0.202732 / 0.136422 | 15423/2755/365204 | complete |
| chain-next-40-4-contiguous | 1.04907 [1.03289,1.29012] | 0.482306 | 2.26232 | 4801 | 3615 | 3614 | 114 | 16 | 0.0108633 / 0.0417121 | 6607/3614/379439 | complete |
| chain-legacy-40-4-contiguous | 1.15477 [1.00442,1.2043] | 0.925111 | 2.32509 | 1261 | 3615 | 3615 | 44 | 30 | 0.326049 / 0.203201 | 19143/3615/576406 | complete |
| chain-next-64-4-contiguous | 3.35231 [2.93385,3.79666] | 1.06557 | 2.35556 | 7681 | 6195 | 6194 | 174 | 16 | 0.0339637 / 0.12167 | 10531/6194/1005775 | complete |
| chain-legacy-64-4-contiguous | 2.91054 [2.68086,3.10013] | 2.04111 | 2.5259 | 1981 | 6195 | 6195 | 68 | 30 | 1.03456 / 0.437393 | 30303/6195/1525244 | complete |
| chain-next-100-4-contiguous | 12.074 [11.8659,12.5341] | 2.81503 | 4.18898 | 12001 | 10065 | 10064 | 264 | 16 | 0.115306 / 0.335303 | 16417/10064/2540536 | complete |
| chain-legacy-100-4-contiguous | 9.85152 [9.48218,10.394] | 4.82658 | 3.45097 | 3061 | 10065 | 10065 | 104 | 30 | 3.64524 / 1.35611 | 47043/10065/3812230 | complete |
| torus-next-6-3-dense | 0.0351112 [0.0348268,0.0385169] | 0.0379648 | 2.08566 | 694 | 24 | 23 | 30 | 17 | 0.000285073 / 0.000476184 | 981/23/5190 | complete |
| torus-legacy-6-3-dense | 0.339865 [0.31707,0.367335] | 0.214545 | 4.59196 | 694 | 4096 | 8192 | 1 | 694 | 3.98207 / 2.15602 | 1926546/8192/1926546 | complete |
| torus-next-6-3-ts | 0.358587 [0.353022,0.37515] | 0.375513 | 2.18153 | 694 | 24 | 23 | 30 | 17 | 0.000225168 / 0.000483864 | 981/23/5427 | complete |
| torus-legacy-6-3-ts | 0.232046 [0.224497,0.248141] | 0.149613 | 2.4684 | 7673 | 1022 | 2044 | 291 | 98 | 1.27716 / 1.52418 | 1363742/2044/1363742 | complete |
| torus-next-9-3-dense | 0.64678 [0.613124,0.658303] | 0.582326 | 2.20012 | 2620 | 75 | 74 | 44 | 25 | 0.00110857 / 0.00356451 | 2323/74/28035 | complete |
| torus-legacy-9-3-dense | 15.3569 [8.94841,16.332] | 3.12162 | 20.026 | 2620 | 104680 | — | 1 | 2620 | — / — | —/—/— | memory_limit at cold_lower |
| torus-next-9-3-ts | 15.0691 [14.5521,15.7102] | 9.38311 | 2.46057 | 2620 | 75 | 74 | 44 | 52 | 0.00200915 / 0.00738701 | 4180/74/61719 | complete |
| torus-legacy-9-3-ts | 12.4667 [7.62261,13.138] | 3.00289 | 6.99465 | 48590 | 26101 | 52202 | 986 | 407 | 80.518 / 31.3804 | 28082498/52202/28082498 | complete |
| torus-next-36-1-dense | 2.30055 [2.13937,2.76064] | 1.68621 | 2.1939 | 109 | 10 | 9 | 11 | 1 | 2.0853e-05 / 1.0992e-05 | 11/9/78 | complete |
| torus-next-36-1-regions | 0.648819 [0.547356,0.653638] | 0.429096 | 2.18035 | 1080 | 4 | 3 | 7 | 2 | 1.9763e-05 / 2.32e-06 | 9/3/8 | complete |
| torus-next-36-1-regions_ts | — | — | 2.47833 | 973 | 4 | 3 | 6 | 2 | — / — | —/—/— | time_limit at warm_build_1 |
| torus-next-36-1-ts | — | — | 2.49981 | — | — | — | — | — | — / — | —/—/— | time_limit at cold_build |
| torus-next-36-2-dense | 6.04284 [5.9993,6.5317] | 2.3499 | 2.69364 | 5779 | 1075 | 1074 | 99 | 19 | 0.00308539 / 0.00890169 | 3642/1074/108427 | complete |
| torus-next-36-2-regions | 0.651565 [0.636708,0.688187] | 0.442711 | 2.17776 | 4968 | 10 | 9 | 21 | 4 | 3.1231e-05 / 2.9016e-05 | 82/9/170 | complete |
| torus-next-36-2-regions_ts | 0.736446 [0.709632,0.777437] | 0.53953 | 2.20058 | 4968 | 10 | 9 | 21 | 4 | 3.3972e-05 / 2.9016e-05 | 82/9/170 | complete |
| torus-next-36-3-dense | — | — | 20.0321 | — | — | — | — | — | — / — | —/—/— | memory_limit at cold_build |
| torus-next-36-3-regions | 0.768684 [0.759687,0.849771] | 0.493311 | 2.18441 | 14688 | 10 | 9 | 29 | 14 | 8.4557e-05 / 0.000194064 | 580/9/1389 | complete |
| torus-next-36-3-regions_ts | 0.937156 [0.882206,1.06593] | 0.560557 | 2.18719 | 14688 | 10 | 9 | 29 | 14 | 8.0405e-05 / 0.000194064 | 580/9/1389 | complete |
| torus-next-36-4-regions | 1.00571 [0.889548,1.15503] | 0.606867 | 2.19863 | 29268 | 10 | 9 | 30 | 30 | 0.00038391 / 0.000535344 | 2260/9/5960 | complete |
| torus-next-36-4-regions_ts | 1.92936 [1.69037,2.05848] | 0.607579 | 2.1903 | 29268 | 10 | 9 | 30 | 30 | 0.000326035 / 0.000535344 | 2260/9/5960 | complete |
| torus-next-36-5-regions | 1.16891 [1.07872,1.25173] | 0.688185 | 2.17349 | 38016 | 10 | 9 | 30 | 40 | 0.000520067 / 0.00106517 | 3906/9/10343 | complete |
| torus-next-36-5-regions_ts | 2.76842 [2.73158,2.84757] | 0.63394 | 2.18647 | 38016 | 10 | 9 | 30 | 40 | 0.000500122 / 0.00106517 | 3906/9/10343 | complete |
| torus-next-36-6-regions | 1.14854 [1.08313,1.25993] | 0.688185 | 2.17985 | 38016 | 10 | 9 | 30 | 40 | 0.000600699 / 0.00106517 | 3906/9/10343 | complete |
| torus-next-36-6-regions_ts | 2.69232 [2.44367,2.73706] | 0.633939 | 2.19203 | 38016 | 10 | 9 | 30 | 40 | 0.000478917 / 0.00106517 | 3906/9/10343 | complete |
| hydrogen-next-4-1-pqg | 0.01449 [0.014284,0.0145325] | 0.0173911 | 2.19503 | 17 | 703 | 1368 | 6 | 65 | 0.000414426 / 0.00118502 | 6922/1368/8980 | complete |
| hydrogen-legacy-4-1-pqg | 0.004492 [0.00410361,0.00571486] | 0.00494132 | 2.20334 | 129 | 1105 | 2210 | 11 | 33 | 0.010673 / 0.0144819 | 9469/2210/14512 | complete |
| hydrogen-next-4-1-t12 | 0.0200078 [0.0182551,0.0266097] | 0.0225445 | 2.18218 | 17 | 703 | 1368 | 8 | 232 | 0.00255203 / 0.00627601 | 63882/1368/42980 | complete |
| hydrogen-legacy-4-1-t12 | 0.0211621 [0.0210268,0.0213621] | 0.037446 | 2.18649 | 417 | 1105 | 2210 | 19 | 92 | 0.0598511 / 0.0904258 | 86525/2210/73040 | complete |
| hydrogen-next-6-1-pqg | 0.0436346 [0.0402354,0.0477218] | 0.0601375 | 2.19431 | 25 | 3160 | 6240 | 6 | 145 | 0.00166398 / 0.00516115 | 32194/6240/40302 | complete |
| hydrogen-legacy-6-1-pqg | 0.0190999 [0.0187177,0.0221395] | 0.0199089 | 2.1797 | 289 | 4351 | 8702 | 11 | 73 | 0.0516661 / 0.0680802 | 46645/8702/64794 | complete |
| hydrogen-next-6-1-t12 | 0.123381 [0.123328,0.156199] | 0.121075 | 2.17444 | 25 | 3160 | 6240 | 8 | 804 | 0.023942 / 0.063533 | 727010/6240/334542 | complete |
| hydrogen-legacy-6-1-t12 | 0.197741 [0.174938,0.199023] | 0.394278 | 2.60792 | 1313 | 4351 | 8702 | 19 | 312 | 0.872027 / 0.900176 | 958197/8702/545386 | complete |
| hydrogen-next-8-1-pqg | 0.115095 [0.113804,0.128885] | 0.147546 | 2.24272 | 33 | 9453 | 18768 | 6 | 257 | 0.00466912 / 0.0163009 | 99146/18768/119720 | complete |
| hydrogen-legacy-8-1-pqg | 0.0677789 [0.0623595,0.0842143] | 0.0581927 | 2.196 | 513 | 11825 | 23650 | 11 | 129 | 0.186916 / 0.199022 | 146413/23650/190868 | complete |
| hydrogen-next-8-1-t12 | 0.443003 [0.424861,0.475263] | 0.43318 | 2.30052 | 33 | 9453 | 18768 | 8 | 1936 | 0.108632 / 0.253215 | 4160842/18768/1439496 | complete |
| hydrogen-legacy-8-1-t12 | 1.79841 [1.21294,1.81786] | 2.17502 | 6.40231 | 3009 | 11825 | 23650 | 19 | 744 | 7.02825 / 4.78395 | 5402605/23650/2296276 | complete |
| bell-legacy-3-2-dense | 0.000945373 [0.000797841,0.00103078] | 0.00100664 | 2.13262 | 28 | 154 | 154 | 1 | 28 | 0.000731087 / 0.000804824 | 407/154/407 | complete |
| bell-next-5-3-dense | 0.00906493 [0.00873106,0.00912725] | 0.00354074 | 2.07906 | 436 | 27556 | 27555 | 1 | 436 | 0.00387083 / 0.011506 | 95266/27555/95265 | complete |
| bell-legacy-5-3-dense | 0.231092 [0.219255,0.269325] | 0.17909 | 2.21472 | 436 | 27556 | 27556 | 1 | 436 | 0.173614 / 0.174413 | 95267/27556/95267 | complete |
| bell-next-5-5-cs_ts | 71.2258 [70.5935,71.9558] | 14.2201 | 13.1151 | 11236 | 11235556 | 11235555 | 1 | 11236 | 10.963 / 10.8017 | 63129466/11235555/63129465 | complete |
| bell-legacy-5-5-cs_ts | 0.75544 [0.706457,0.801234] | 0.378987 | 2.40592 | 948 | 23838 | 23838 | 5 | 628 | 0.511593 / 0.366165 | 217667/23838/217667 | complete |
| bell-next-6-3-dense | 0.0344437 [0.0290747,0.0359828] | 0.0178853 | 2.09027 | 769 | 91969 | 91968 | 1 | 769 | 0.0139274 / 0.0530111 | 296065/91968/296064 | complete |
| bell-legacy-6-3-dense | 1.00095 [0.895269,1.12764] | 0.570314 | 2.61887 | 769 | 91969 | 91969 | 1 | 769 | 0.830994 / 0.550618 | 296066/91969/296066 | complete |

## Solve checks

| Case | Status / stop | Bound | Solve call s | Whole process s | Peak GB |
|---|---|---:|---:|---:|---:|
| hydrogen-next-4-1-pqg-solve | OPTIMAL | -2.168620419132367 | 1.03482 | 46.1293 | 2.05388 |
| bell-next-3-2-dense-solve | OPTIMAL | 0.25093973126127517 | 8.92676 | 40.3112 | 2.12351 |
| bell-legacy-3-2-dense-solve | OPTIMAL | 0.25093973173010703 | 4.92423 | 46.0824 | 2.25124 |
| bell-next-5-3-dense-solve | memory_limit | — | — | 55.6158 | 20.0879 |
| bell-next-5-5-cs_ts-solve | memory_limit | — | — | 161.497 | 20.0307 |
| bell-next-6-3-dense-solve | memory_limit | — | — | 56.112 | 20.0147 |

## Largest completed points and limits

- **Chain:** all requested sizes through N=100 finished on both routes. Next
  used 12.074 s and 2.815 GB allocated per warm build; legacy used 9.852 s and
  4.827 GB. Next was not faster at this point, and its whole-process peak was
  higher (4.189 versus 3.451 GB). Next standard-form lowering took 0.115 s;
  legacy native JuMP lowering took 3.645 s. Both have 10,065 moments but emit
  different cone decompositions. N=100 is the requested endpoint, not a limit.
- **Recorded roster:** all 3×2 and 3×3 order-3 builds finished. Legacy 3×3 dense
  lowering exceeded 20 GB; its build alone completed. The separate successful
  build profile counted 209,360 native-writer variables for this row. Next 3×3
  dense took 0.647 s, versus 15.069 s with default term splitting. Legacy MMD
  term splitting took 12.467 s, then 80.518 s and 31.380 GB cumulative allocation
  for its native lowering. The split algorithms and symmetry reductions differ.
- **Scale point, not a QMBCertify reproduction:** global 6×6 order 2 finished
  (6.043 s, 2.350 GB allocated, 2.694 GB peak); global order 3 crossed 20 GB in
  the cold build after 222.6 s. Global `TermSplit()` order 1 hit 600 s in the
  cold build; higher global split orders were not attempted. Radius-one regions
  finished through orders 5 and 6, with and without splitting. These five-site
  Pauli regions saturate at order 5: order 6 has the same row and cone sizes.
  This is not a global order-5 relaxation. Without splitting, order 5 took
  1.169 s and 0.688 GB allocated. Regional splitting at order 1 was unusually
  expensive: one cold build completed in 498.383 s with 359.464 GB cumulative
  allocations, then the process reached 600 s during the first warm repeat.
  Orders 2–6 completed; order-1 failure does not imply failure at higher orders.
- **Hydrogen:** every H4/H6/H8 PQG and T1/T2prime build and lowering finished.
  H8/T1/T2prime took 0.443 s and 0.433 GB allocated in Next, versus 1.798 s and
  2.175 GB in the legacy structured route. Next lowering took 0.109 s versus
  7.028 s for native legacy JuMP lowering. H8 is the requested endpoint. The
  sector equations, spin blocking, and cone representations differ; this is not
  a matched-solver speedup claim.
- **Bell:** all roster builds finished. The largest, n=5/order-5 with correlative
  and term splitting, emits one side-11,236 cone in Next, with 11,235,556 moments,
  71.226 s warm build time, 14.220 GB allocated and 13.115 GB peak. Legacy emits
  five cones (largest side 628), with 23,838 moments. Next's correlative grouping
  acts on the two party sites, not the measurement letters. Thus this sparse
  comparison is **not the same relaxation**. All three Next roster solve attempts
  exceeded 20 GB inside `JuMPBackend`; no optimizer-startup output was recorded.
  The precise allocation site was not profiled. No roster row therefore qualified
  for a same-relaxation paired solve; legacy roster solves were not attempted.
  The smaller dense n=3/order-2 control returned 0.25093973126127517 (Next) and
  0.25093973173010703 (legacy), difference −4.6883e−10. No literature
  Bell value is claimed.

The H4/PQG total energy is −2.168620419132367 Eh, 0.0022329704976025 Eh below
the committed FCI value −2.1663874486347643 Eh. This is a relaxation lower bound,
not FCI equality or a rigorous certificate. SCS returned `OPTIMAL` at tolerance
1e−5; primal residual 9.43135e−6, dual residual 2.68511e−5, gap 1.72347e−6,
425 iterations. No other hydrogen solve was run.

## Build profiles identify costs, not optimizations

Separate profiles use the largest completed requested chain, hydrogen and Bell
builds; the torus profile uses the highest completed **global** 6×6 order on Next
and the largest recorded 3×3 order on legacy. Regional orders are incomparable
with global orders. Sampling interval is 2 ms, after compilation and three warm
builds, with enough fresh builds to sample at least about 3 s (maximum 100 builds).
The raw flat profiles include inclusive and self counts. The paths below are
nested, so their costs must not be added.

| Profile | Three large cost paths / functions | Ownership |
|---|---|---|
| Next chain N=100 | `_fourier_kernel!` → `_orbit_pair_kernel!`; `orbit_coordinates!` and its dictionary work; `act(::CliffordMap, ...)` within orbit traversal | `pipeline/lattice.jl`, `pipeline/moments.jl`, `core/actions.jl` |
| Legacy chain N=100 | `_translation_momentum_block`; `_mirror_real_blocks`; `_build_moment_linear_data` | `optimization/pauli_chains.jl`, `optimization/moment.jl` |
| Next global 6×6/order 2 | `_fourier_kernel!`; `shared_projections` → `CharacterTable`; `compose(::CliffordMap, ...)` within conjugacy-class construction | `pipeline/abelian.jl`, `pipeline/finite.jl`, `core/actions.jl` |
| Legacy 3×3/order 3 dense | `_build_moment_linear_data`; `_build_constraint_matrix` / `_shrink_terms!`; `_polynomial_total_basis` | `optimization/moment.jl`, `types/polynomial.jl` |
| Next H8/T1/T2prime | fermionic `_mul!`; `_close_orbit_form!`; `_sector_equalities!`, reached through `_fermionic_strengthen!` / `_append_family_sectors!` | `core/fermi.jl`, `pipeline/orbits.jl`, `pipeline/sectors.jl` |
| Legacy H8/T1/T2prime | `_v2rdm_three_index_entry_form`; `normal_order_terms`; `_assert_moment_linear_data_invariants` | `optimization/v2rdm_structured.jl`, `optimization/moment_linear.jl` |
| Next Bell n=5/order 5 | `_block_product!`; `_split_family` / `_support_graph`; `_close_orbit_form!` | `pipeline/families.jl`, `pipeline/selection.jl`, `pipeline/orbits.jl` |
| Legacy Bell n=5/order 5 | `_build_moment_linear_data`; `_build_constraint_matrix`; `_polynomial_total_basis` | `optimization/moment.jl` |

For scale: the Next chain profile has 5,308 inclusive samples in `_fourier_kernel!`
out of 5,839 build-stack samples. The 6×6 profile has 1,297 in `_fourier_kernel!`
and 1,071 in `shared_projections`, out of 2,673. The Bell profile has 24,697 at
the `_block_product!` call path and 12,965 at the splitting call site, out of
36,498; these paths overlap. These are sampled CPU costs, not allocation profiles.
The slow 6×6 split order-1 case was not profiled, so no root cause is asserted.

Large build batches took about 14, 34, 24 and 18 minutes, within their respective
10–25, 20–40, 15–30 and 15–30 minute estimates. The final capped solve batch took
about 5 minutes (estimate 5–30), followed by 17 minutes of profiles (estimate
10–20). Individual process elapsed times and stop stages are in `results.json`.
