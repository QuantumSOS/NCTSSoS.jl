# DEV-382 Track B: local invariant cones at unchanged moment order

A scalar moment is the real expectation of a Hermitian Pauli word. A
second-order cone (SOC) constrains a vector's Euclidean norm. This pilot adds
SOCs on scalar moments already present in an ordinary Pauli moment relaxation;
it adds neither basis words nor moment variables. It does not implement a
Jordan-algebra engine, weighted separation, or a state/two-copy hierarchy.

## Validity does not require a symmetry of the Hamiltonian

Let a cluster contain \(k\in\{3,4\}\) qubits. Label a Pauli word by binary
vectors \(v=(x,z)\), with Y represented by \(x_j=z_j=1\). Define
\(q_0(v)=x\cdot z\pmod2\), the parity of its Y count, and the symplectic
pairing \([r,v]=r_x\cdot z+r_z\cdot x\pmod2\). Complex conjugation in the
computational basis is denoted by \(K\). For a Pauli word \(P_r\), the
antilinear norm-preserving map \(\Theta_r=P_rK\) satisfies
\(\Theta_r^2=(-1)^{q_0(r)}I\). Only odd-Y labels \(r\) are eligible.

Set \(q_r(v)=q_0(v)+[r,v]\pmod2\) and \(J_r=\{v:q_r(v)=0\}\).
For any physical cluster density matrix \(\rho\), the positive, trace-one
matrix
\[
\bar\rho=\tfrac12(\rho+\Theta_r\rho\Theta_r^{-1})
\]
retains precisely its invariant Pauli expectations. Its eigenvalues come in
pairs: if \(\bar\rho\psi=\lambda\psi\), then \(\Theta_r\psi\) is an
orthogonal eigenvector with the same eigenvalue because \(\Theta_r^2=-I\).
Consequently \(\operatorname{tr}\bar\rho^2\le1/2\). Orthogonality of the
Pauli basis gives
\[
1+\sum_{v\in J_r\setminus\{0\}}y_v^2
=2^k\operatorname{tr}\bar\rho^2\le2^{k-1}.
\]
Dropping nonnegative summands proves, for **any available subset**
\(S\subseteq J_r\setminus\{0\}\),
\[
\|(y_v)_{v\in S}\|_2\le\sqrt{2^{k-1}-1}.
\]
The radii are therefore √3 and √7. This argument applies to the reduced
density matrix of every cluster, whether or not the Hamiltonian is invariant
under the chosen map. Purity is used only in this proof, not as an optimization
variable or an additional hierarchy.

## Separation and audit records

`pauli_invariant_soc_nctssos(pop, config)` constructs an unreduced Pauli
moment model. `separate_pauli_invariant_socs!(model, mp, clusters)` also accepts
an existing matching model, allowing composition with anticommutation cones.
`add_pauli_invariant_soc!` supports explicit subsets. The default direct
`build_jump_model` lowering is required; symmetry-identified moments and
state polynomials are deliberately unsupported.
Cluster site numbers are positions in the Pauli registry, not arbitrary
display subscripts passed to `create_pauli_variables`.

Clusters grow as connected subsets of the Hamiltonian's site-interaction
graph, restricted to correlative-sparsity cliques. Three- and four-site
clusters are both generated. Growth is deterministic and capped at 64 subsets
per size, including the intermediate size-two frontier; it is a local
heuristic, not exhaustive large-system separation. Explicit clusters may be
supplied instead.

For available nonidentity words, define \(a_v=(-1)^{q_0(v)}y_v^2\), and set
missing entries to zero. A binary Walsh–Hadamard transform computes signed
sums; swapping its output index's x/z halves yields
\(W_r=\sum_v(-1)^{[r,v]}a_v\). Thus the squared norm for candidate \(r\) is
\[
\tfrac12\left(\sum_v y_v^2+W_r\right).
\]
The transform takes \(O(k4^k)\) time and \(O(4^k)\) storage per cluster.
The separator ranks squared-norm violations, excludes even-Y labels, adds
at most two cones per round in this study, and re-solves for up to three
rounds. It stops early when no new violation exceeds 1e-6. Equivalent selected
sets within a cluster are deduplicated. Every returned record contains the
cluster, `r`, selected Pauli strings in cluster-site order, radius, squared
radius, round, pre-addition violation, and JuMP constraint reference.

## Frozen examples and search provenance

**Supplied:**
\(H=-(X_1Y_2+X_1Y_3+X_2X_3+Y_2Y_3+Z_2Z_3)\).
The degree-one relaxation gives −5. One `r=YII` cone has radius √3 and selects
`IXI, IIX, IXX, IZI, XYI, IZX, IIZ, IXZ, XIY, IZZ, IYY`.
The five objective terms belong to that set, so Cauchy–Schwarz proves the
lower bound −√15 ≈ −3.872983346. The solved relaxation attains it to tolerance.
The exact 8×8 matrix has minimum eigenvalue −3. Exhaustive anticommutation
does **better** here, about −3.828427125; this negative comparison is retained.

**Rich3:** write
\(H_3=\sum_{a,b\in\{X,Y,Z\}}(J^{12}_{ab}a_1b_2+
J^{13}_{ab}a_1b_3+J^{23}_{ab}a_2b_3)\).
Rows and columns below are ordered X, Y, Z. This is the full coefficient
table, including the zero coefficient:

| Pair | X row | Y row | Z row |
|---|---|---|---|
| 1,2 | 0.4, −1.0, 1.0 | −1.0, −0.8, −0.6 | −0.6, −1.0, −0.4 |
| 1,3 | −1.0, 0.8, −0.8 | 0.6, −0.6, −0.2 | −0.6, −0.8, −0.8 |
| 2,3 | 0.8, 0.4, 0.8 | 0.0, 0.8, 0.4 | −1.0, 0.8, −0.4 |

The bounded search used Julia 1.12.7, `Random.seed!(382)`, and
`rand(-5:5, 27) ./ 5` assigned to degree-two words in `get_ncbasis` order.
It was capped at 40 draws each for three and four sites. The **first**
three-site draw qualified, so the search stopped after 31.04 seconds including
compilation. The frozen integer coefficient vector is in the script and CI
fixture; replay does not depend on random-number-generator changes. The search
solved the exhaustive anticommutation model, scored all eligible local cones,
and, when a squared-norm violation exceeded 1e-4, added all eligible invariant
cones and accepted an objective increase above 1e-3. The limited-batch replay
recovered the same combined optimum with just one `r=YYY` cone. Its selected
set is all 27 weight-two words, radius √3. The unchanged degree-one universe
contains 37 moments: identity, nine single-site words, and 27 two-site words.

**Rich4:**
\(H_4=H_3-0.7(X_1Y_4+X_2X_4+Y_2Y_4+Z_2Z_4)+0.23Z_1+0.17X_4\).
It has 67 degree-one-model scalar moments. The fields make the model more
general than a time-reversal-invariant two-body Hamiltonian. Both automatic
three-/four-site separation and four-site-only separation are measured.

## Baselines and certification

All methods use the repository's default no-elimination sparsity configuration,
not manually forced dense matrix blocks. The supplied problem has 23 available
scalar moments; Rich3 has 37 and Rich4 has 67. Comparisons preserve this same
construction except for the explicitly enlarged bases.

Two operators anticommute when \(AB=-BA\). The anticommutation graph has one
vertex for each available scalar Pauli
moment and an edge for each anticommuting pair. **All maximal cliques** of
size at least two receive \(\|y_C\|_2\le1\). This is not a selected-edge or
objective-support-only baseline. The raw records include all clique vertex
indices and their complete Pauli moment list. Rich3 has 363 cones; Rich4 has
1,324. CI independently constructs this graph with explicit Pauli matrices.

There are two separate comparisons:

1. **Fixed degree one:** bare, exhaustive anticommutation, invariant SOC, and
   their combination share exactly the same scalar moment set and variables.
2. **Larger basis:** `extra1`, `extra3`, and `extra5` append the first one,
   three, or five degree-two words in canonical basis order. They introduce
   additional moments and are **not fixed-degree-one comparisons**. Selection
   is deterministic, not tuned to coefficients or solver results. The complete
   selected bases are recorded. More basis words can win; no universal
   superiority over enlargement is claimed.

The solver's primal/dual gap must be below 1e-6 for every reported sample.
In addition, the benchmark repairs the final numerical dual solution using
exact `Rational{BigInt}` arithmetic. A positive-semidefinite (PSD) matrix is
one with nonnegative quadratic forms. Each PSD multiplier is replaced by
\(BB^T\), with rational B; every SOC multiplier is rational and checked by
the exact inequality \(t^2\ge\sum_i v_i^2\), \(t\ge0\). Irrational SOC
radii are replaced by rational upper bounds verified by squaring. Subtracting
these valid multiplier-weighted constraints from the objective leaves an
affine residual. Identity expectation is exactly one; imaginary Pauli
expectations vanish; every other Pauli expectation lies in [−1,1]. Subtracting
the sum of absolute residual coefficients therefore yields a **certified
quantum lower bound**, not merely an optimal-status assertion. The exact
rational bound is saved alongside its floating display. Certification is for
the implemented Float64 Hamiltonian coefficients interpreted exactly.
The audit also records the exact affine constant, every real residual
coefficient in the saved moment order, the identity index, and the exact and
floating sum of absolute nonidentity residual coefficients.
Writing that residual as \(a+\sum_v\delta_v y_v\), its certified value is
\(a+\delta_I-\sum_{v\ne I}|\delta_v|\). Imaginary-moment residuals cost
nothing because physical Pauli expectations are real. The nonidentity
residual sums are 4.18e-11 for the supplied single cone, 3.12e-8 for Rich3
SOC-only, 1.75e-8 for Rich3 combined, and 3.18e-8 for Rich4 combined. The
additional conservative cost of rational cone-radius/multiplier rounding is
also included in each certificate, not omitted from these lower bounds.

Matrix diagonalization supplies an independent physical-energy reference,
not a relaxation oracle. COSMO 0.8.11 reproduces the primal objectives in CI,
but its reported dual objective on the supplied mixed-cone model is loose
(−4.54328433 versus primal −3.87298335), with and without chordal decomposition.
The study therefore uses licensed Mosek and the independent rational repair;
it does not silently substitute COSMO primal objectives for certified bounds.

## Measurement protocol and results

Environment: Linux x86-64 Amp orb, Intel Xeon @ 2.60 GHz, Julia 1.12.7,
Mosek 11.2.7, MosekTools 0.15.10; one solver thread and one BLAS thread.
Mosek primal-feasibility, dual-feasibility, and relative-gap tolerances are
1e-9. A trivial solve verified licensing before measurements.

Each case/method runs in a separate Julia process, warms up once, then builds
and solves three independent models. Construction includes basis generation,
sparsity, symbolic moments, JuMP lowering, anticommutation generation where
applicable, and cluster generation. Separation times include scoring and cone
addition; solve times sum every `optimize!` call, including intermediate
relaxations. Exact diagonalization and rational certification are outside
these timers. Peak resident memory is the **whole-process high-water mark**,
including compilation, all trials, and certificate checking, not incremental
solver memory. Shared-orb noise and compiler-dominated memory preclude a
speedup or solver-memory advantage claim. Raw samples and construction
allocation counts are retained in `perf/results/dev382/*.toml`.

The following times are warmed medians in milliseconds. Parentheses give the
interquartile range (IQR), the 75th percentile minus the 25th, using linear
interpolation over the three samples; these are not confidence intervals.
Memory uses binary mebibytes (MiB). Bounds are rounded displays of the exact
rational certificate. The isolated sweep took 1,105 seconds (18 min 25 s),
within its 15–20 minute estimate. The largest primal/dual gap across all 57
measured final primal/dual pairs was 8.69e-9.

| Case | Method | Certified lower bound | Build ms (IQR) | Separate ms (IQR) | Solve ms (IQR) | Peak RSS MiB |
|---|---|---:|---:|---:|---:|---:|
| supplied | bare | −5.00000000 | 2.654 (0.412) | 0.000 (0.000) | 7.012 (1.086) | 1833.0 |
| supplied | anti | −3.82842722 | 3.020 (0.635) | 0.000 (0.000) | 15.630 (1.107) | 1833.9 |
| supplied | one_soc | −3.87298335 | 1.967 (0.137) | 0.218 (0.033) | 10.001 (1.445) | 1909.5 |
| supplied | extra1 | −4.30097600 | 2.055 (0.146) | 0.000 (0.000) | 9.577 (1.010) | 1925.2 |
| supplied | extra3 | −3.98927571 | 2.617 (0.488) | 0.000 (0.000) | 12.549 (1.004) | 1833.0 |
| supplied | extra5 | −3.00000000 | 2.942 (0.407) | 0.000 (0.000) | 16.095 (1.216) | 1835.2 |
| rich3 | bare | −7.18522231 | 2.276 (0.499) | 0.000 (0.000) | 12.745 (0.294) | 1905.5 |
| rich3 | anti | −6.23562216 | 5.598 (0.100) | 0.000 (0.000) | 76.260 (1.327) | 1905.0 |
| rich3 | soc | **−6.09384688** | 2.455 (0.137) | 0.428 (0.071) | 25.303 (2.309) | 1912.9 |
| rich3 | combined | **−6.07499146** | 4.528 (0.312) | 0.276 (0.007) | 139.726 (7.696) | 1831.5 |
| rich3 | extra1 | −6.63898921 | 2.319 (0.347) | 0.000 (0.000) | 16.540 (1.927) | 1821.5 |
| rich3 | extra3 | **−6.01121292** | 3.141 (0.367) | 0.000 (0.000) | 20.716 (2.124) | 1835.1 |
| rich4 | bare | −9.94307688 | 3.207 (0.529) | 0.000 (0.000) | 23.760 (0.808) | 1834.1 |
| rich4 | anti | −7.84926475 | 11.739 (0.134) | 0.000 (0.000) | 319.320 (4.199) | 1824.2 |
| rich4 | soc | −8.08930147 | 3.249 (0.548) | 0.710 (0.092) | 63.126 (3.245) | 1896.7 |
| rich4 | combined | **−7.81771427** | 12.002 (1.700) | 0.450 (0.061) | 671.336 (105.859) | 1838.7 |
| rich4 | extra1 | −8.94063910 | 2.469 (0.073) | 0.000 (0.000) | 29.108 (2.309) | 1916.2 |
| rich4 | extra3 | −8.25187874 | 3.461 (0.365) | 0.000 (0.000) | 67.793 (2.558) | 1838.1 |
| rich4 | soc4 | −9.37178602 | 2.743 (0.075) | 0.850 (0.026) | 72.800 (3.097) | 1954.9 |

**Primary success:** Rich3's SOC-only certificate improves over exhaustive
anticommutation by 0.141775 and over the one-word enlargement by 0.545142.
Against that enlargement, warmed construction costs differ by +5.9% and
whole-process RSS by +5.0%; separation adds 0.428 ms, and cumulative solving
costs +53.0%. These are observed resource deltas, not speedup claims. The
degree-one SOC model retains 37 scalar moments and 74 real variables;
`extra1` has 40 moments and 80 variables. The extra word is `XXI`.

Adding just the `YYY` cone to exhaustive anticommutation improves the
certificate by 0.160631, with no moment change. Rich3 SOC-only adds `YYY` and
`ZIY` in one batch; their full selected sets are in `rich3-soc.toml`.

**Counterexample:** adding `XXI, XYI, XZI` (`extra3`) wins on Rich3, reaching
−6.01121292. Its construction median is 28.0% above SOC-only, its solve time
is lower, and its compiler-dominated RSS is lower. The pilot establishes
additional fixed-order strength and one close-cost comparison, **not** a
general resource advantage or dominance of the best basis-enlargement choice.

Rich4's combined certificate improves over exhaustive anticommutation by
0.0315505. The standalone four-site-only search improves the bare bound but
is much weaker than using three-site subclusters. The independent matrix
ground energies are −3 (supplied), −4.83107523520 (Rich3), and −5.95031888322
(Rich4), so substantial relaxation gaps remain.

## Reproduction and coverage

Use a separate environment so optional benchmark dependencies do not change
the library's dependency list:

```sh
julia --project -e 'using Pkg; Pkg.activate("/tmp/dev382-env"); Pkg.develop(path=pwd()); Pkg.add(["MosekTools", "COSMO"])'
JULIA_LOAD_PATH='/tmp/dev382-env:@:@stdlib' julia --project perf/pauli_invariant_soc.jl rich3 combined result.toml
JULIA_LOAD_PATH='/tmp/dev382-env:@:@stdlib' julia --project -e 'include("test/TestUtils.jl"); include("test/relaxations/pauli_invariant_soc.jl")'
```

A Mosek license is required for the benchmark, not the COSMO CI tests.
The regression belongs in **CI**, while the heavier resource comparison and
derivation belong in this benchmark report rather than duplicate executable
website examples. Tests cover the supplied analytic result, exact spectrum,
unchanged moments, arbitrary selected subsets, invalid labels and limits,
three-/four-qubit scores against direct matrix conjugation, and composition
with exhaustive anticommutation constraints. No remote deployment or merge
is part of this pilot.
