# SU(2)-Invariant Moment Quotient Design

## Status

Approved architecture for closing the remaining Phase 3 solver-formulation
gate in `plan/qmbcertify_parity.md`.

## Problem

The translation/SU(2) fast path already computes the intended reduced PSD
blocks.  The remaining scale failure occurs after that reduction because the
block entries are expanded back into canonical Pauli-word moment coordinates
before JuMP lowering.

For the guarded L=30 SU(2)+RDM+state-opt profile, the current formulations
have both failed:

- Native Hermitian SOS dual: 122,385 scalar coefficient equations, about
  120 GiB dense-Schur storage proxy, and a 30-minute timeout at about 81 GiB
  peak RSS.
- Primal moment variables: 122,386 moment variables and only 821 explicit
  scalar equalities, but an N=12 solve still timed out after 20 minutes at
  about 26 GiB peak RSS.

Changing primal to dual moves the same excess coordinate space across the
conic boundary.  Solver tuning cannot remove it.

The current tree also provides positive evidence for a quotient.  The N=10,
order-4 analytic translation target has 121 orbit coordinates but only six
SU(2) singlet channels.  The SU(2) transforms, support-orbit grouping,
coefficient-domain labels, Wigner-Eckart blocks, and singlet-channel selectors
already exist.  The missing mechanism is to use the singlet channels as the
moment coordinate space rather than merely emitting non-singlet zero rows.

## Goal

Represent every supported SU(2)-invariant Pauli moment form in a minimal,
deterministic set of canonical pivot moments spanning the singlet channels,
then lower the already-reduced PSD blocks in those coordinates.

The quotient must preserve the mathematical relaxation.  It is not a new
relaxation, a numerical rank approximation, or a solver-specific heuristic.

## Scope

This design covers the active 1D QMBCertify parity plan through Phase 4:

- periodic Pauli chains;
- translation and optional reflection reduction;
- base SU(2) Wigner-Eckart moment blocks;
- extending full, U(1), and SU(2) contiguous RDM blocks;
- linear and PSD state-optimality add-ons;
- moment and singlet-channel scalar equalities;
- native Hermitian SOS dualization and certificate provenance.

The generic algebra path remains unchanged.

The following are separate projects and are not part of this design:

- Phase 5 two-dimensional lattice symmetry;
- interval/rational certificate reconstruction;
- approximate QR/SVD compression of arbitrary moment problems;
- new non-Pauli continuous-group reducers.

## Mathematical model

Consider one translation support orbit with a complete ordered family of
canonical Pauli-axis words.  Let `y` be its canonical moment vector and let
`S` contain the orthonormal spin-zero rows returned by the existing SU(2)
transform code.  SU(2) invariance is exactly

```text
y = S* z
```

where `S*` is the adjoint and `z` contains the singlet-channel coordinates.
All non-singlet rows vanish.

`MomentLinearData` currently requires coordinate keys to remain canonical
moment keys with representative monomials.  Introducing abstract channel keys
would make `key_to_monomial`, adjoint closure, extraction, and certificate
provenance lie about what a coordinate means.  Instead, choose a deterministic
set of canonical pivot moments `p = P y`, one per singlet channel.  Since the
number of pivots equals the known singlet multiplicity and `P S*` is square,

```text
z = (P S*)^-1 p
y = S* (P S*)^-1 p.
```

The reconstruction matrix

```text
R = S* (P S*)^-1
```

maps every canonical moment in that support orbit to the retained canonical
pivot moments.  Any linear form `a'y` is rewritten as `a'R p`.

The singlet multiplicity comes from representation labels, not from a
floating-point rank threshold.  Numerical linear algebra is used only to pick
and validate a coordinate chart of that known dimension.  It must never decide
how many channels exist.

## Pivot selection

For each support orbit:

1. Order canonical axis words with the existing `key_lt`/basis ordering.
2. Read the known singlet-channel count from the SU(2) labels.
3. Select exactly that many pivot columns using deterministic column-pivoted
   QR of `S*`.
4. Preserve the selected canonical keys as quotient coordinate keys.
5. Form `R` by solving against the selected square submatrix; do not form a
   generic pseudoinverse.
6. Validate all of the following before accepting the map:
   - selected pivot count equals the known channel count;
   - selected submatrix is nonsingular;
   - `P R` is the identity to `quotient_atol`;
   - `(I - S*S) R` is zero to `quotient_atol`;
   - reconstructing `S*` through the pivots is within `quotient_atol`;
   - the condition estimate is finite and below a documented hard ceiling.

Failure is an error with the support orbit, channel count, residuals, and
condition estimate.  There is no silent rank truncation.

## Scalable singlet transforms

The existing `pauli_su2_word_transform` constructs and diagonalizes the full
`3^s` by `3^s` spin-1 tensor representation.  It is a valuable diagnostic for
small support, but it is not the quotient implementation: support 8 has 6,561
axis words, exceeds the current 4,096-dimension guard, and would make a dense
full-transform eigendecomposition the new bottleneck.

Add a cached singlet-only Clebsch-Gordan recurrence for fixed active support:

- couple one local spin-1 letter at a time;
- retain deterministic coupling-path labels and square-root-rational
  coefficient provenance;
- materialize only the final `J=0, m=0` rows over Cartesian Pauli words;
- store the result as sparse column forms or a thin matrix with
  `singlet_count` rows and `3^s` columns;
- reuse the same coefficient matrix for every geometric support orbit of size
  `s`;
- never materialize a dense `3^s` by `3^s` transform in the quotient path.

The expected fixed-support singlet multiplicities for support sizes 0 through
8 are `[1, 0, 1, 1, 3, 6, 15, 36, 91]`.  For a full eight-site local operator
space, identities on inactive sites give 1,430 scalar channels, matching the
sum of squared Schur-Weyl multiplicities of the existing k=8 SU(2) RDM blocks.

For support sizes 0 through 4, the recurrence must match the existing dense
transform up to row phase/unitary rotations within repeated singlet spaces.
For support sizes 5 through 8, validation uses the known multiplicity,
orthonormality, total-`J_z=0`, and annihilation by the total raising operator.
The cached support-8 rows and reconstruction matrix must remain thin; the test
suite rejects any implementation that allocates a full 6,561-square matrix.

## Components

### `src/optimization/pauli_su2_quotient.jl`

Add a focused internal module file rather than extending the already-large
`pauli_chains.jl` further.  It owns:

- `PauliSU2MomentQuotientDescriptor`: immutable provenance for site count,
  support orbit, source keys, pivot keys, singlet labels, reconstruction
  coefficients, coefficient domains, residuals, and condition estimate;
- support-orbit collection from `MomentLinearData.key_to_monomial`;
- deterministic pivot selection and reconstruction-map validation;
- cached singlet-only Clebsch-Gordan rows through active support 8;
- sparse canonical-key-to-pivot forms;
- form, PSD-block, zero-row, and objective rewriting;
- construction of a new invariant-checked `MomentLinearData` through
  `MomentLinearBuilder`.

The file is included after `pauli_chains.jl`, because it reuses the existing
Pauli SU(2) transform and translation-support helpers.  Calls from
`pauli_chains.jl` resolve after module initialization without duplicating the
representation code.

### `MomentLinearData` integration

The quotient is a pure transformation:

```text
raw MomentLinearData
  -> validate SU(2)-invariant profile
  -> build support-orbit reconstruction maps
  -> rewrite objective, PSD entries, and scalar rows
  -> drop rows that become structurally zero
  -> finalize a new MomentLinearData
  -> run existing invariants
  -> lower with existing JuMP/SOS code
```

No mutable post-finalization edits are permitted.  Rebuilding through
`MomentLinearBuilder` preserves key ownership, adjoint closure, pivot
discovery, block metadata, and invariant checks.

### Metadata

Each rewritten block keeps its existing translation/reflection/RDM/SU(2)
origin and gains a quotient descriptor.  The descriptor records both
`:complex_algebraic_float64` coefficients and the existing
`:complex_sqrt_rational` exact-domain tag.

Zero-row provenance is preserved for rows that remain nonzero.  Wigner-copy,
Wigner-offblock, and explicit singlet-channel rows that become identities are
counted by origin before being removed.  Certificate diagnostics must report
that count rather than making the rows disappear without explanation.

### Public integration

Add `su2_moment_quotient::Bool` to the supported Pauli-chain fast-path options.
During development it is explicitly selectable so the old and new models can
be compared in the same test process.  After the equivalence gates pass, it is
enabled by default for the supported direct-linear base-SU(2) path.

The option is valid only when the recognizer has proved the required Pauli,
translation, and SU(2)-invariance conditions.  Supplying it on an unsupported
profile raises `ArgumentError`; it does not activate the generic path and does
not weaken constraints.

`SolverConfig`, `cs_nctssos`, the low-level translation constructor, the
profile harness, and the solver-probe harness must agree on the option.

## Failure behavior

The quotient fails closed when any of these conditions holds:

- the algebra is not Pauli;
- SU(2) invariance was not established;
- a canonical key cannot be assigned to a translation support orbit;
- canonicalization makes the support family inconsistent;
- the expected singlet multiplicity disagrees with transform labels;
- the selected coordinate chart is singular or ill-conditioned;
- reconstruction or orthogonality residuals exceed tolerance;
- rewritten forms violate adjoint/self-adjoint invariants;
- an objective or constraint requires a moment outside the quotient map.

Errors identify the feature and support orbit.  The generic fallback is chosen
only when the fast-path recognizer rejects the profile before quotient
construction.  Once the supported SU(2) fast path is selected, quotient
failure remains visible; default activation never turns a construction error
into a silent fallback.

## Solver form

Use the reduced coordinate model with native Hermitian PSD blocks and the
existing SOS dualizer.  This follows the natural sparse-LMI form:

- keep the already-computed small PSD blocks separate;
- do not real-lift Hermitian blocks unless the optimizer lacks native support;
- do not introduce PSD slacks aligned by thousands of binding equalities;
- do not keep canonical variables whose only role is to be projected out;
- do not replace exact symmetry structure with first-order solver tolerances.

The existing primal moment-variable path remains a diagnostic comparison, not
the default Phase 3 scale path.

## Metrics and harness evidence

Extend `TranslationInvariantReport` and `translation_report_metrics` with:

- raw linear moment count;
- quotient moment count;
- quotient reduction ratio;
- support-orbit count;
- singlet-channel count by support size;
- maximum reconstruction residual;
- maximum pivot identity residual;
- maximum coordinate-chart condition estimate;
- Wigner/singlet zero rows eliminated by quotient;
- remaining scalar equality count;
- actual SOS coefficient equation count;
- updated dense-Schur storage proxy;
- quotient enabled/supported status and blocker provenance.

The profile and solver-probe harnesses print these fields.  Persisted fixture
rows include the quotient option and counts so a future run cannot be mistaken
for the old 122k-coordinate formulation.

## Test strategy

Every behavioral change follows red-green-refactor through `easy-ssh`.

### Unit tests

- support sizes 0 through 8: expected singlet multiplicities and deterministic
  pivot keys;
- support sizes 0 through 4: singlet-only recurrence agrees with the existing
  dense transform up to the expected internal singlet-basis freedom;
- support sizes 5 through 8: rows are orthonormal, have total `J_z=0`, and are
  annihilated by total `J_+` without allocating a full square transform;
- full k=8 local operator accounting reports 1,430 scalar channels;
- reconstruction fixes every selected pivot exactly within tolerance;
- reconstructed canonical vectors lie in the singlet subspace;
- non-singlet vectors are rejected or projected only when the caller has
  established SU(2) invariance;
- permutation of input key order produces the same ordered pivots and map;
- singular/ill-conditioned synthetic charts fail with useful errors;
- owned keys survive scratch-buffer mutation;
- adjoint maps and self-adjoint scalar rows remain valid.

### Model equivalence tests

At small valid chains with `N > 2d`, compare quotient disabled/enabled for:

- base translation/SU(2);
- reflection-fixed and conjugate momentum sectors;
- closed and extending full/U(1)/SU(2) RDM blocks;
- linear state optimality;
- PSD state optimality;
- moment equalities;
- singlet-channel equalities;
- combined RDM+LSO+PSO profile.

Checks include block labels and sizes, objective forms, surviving zero-row
origins, solver objective within `1e-8` for small stable cases, and SOS
certificate residual at the existing solver tolerance.

### Rejection tests

- non-Pauli algebra;
- missing SU(2) assumption;
- incomplete/ambiguous support orbit;
- non-invariant objective;
- unsupported public option combinations;
- residual and condition-limit failures.

### Performance ladder

Large runs are never queued merely because code exists.  Each rung requires
fresh `uptime`/`free -h`, one Julia thread, one MOSEK thread, a hard timeout,
and explicit wall/RSS estimates.

1. N=10 target-only: quotient target counts agree with analytic singlet counts.
2. N=8 or N=10 constructed no-solver case: raw/quotient forms and residuals.
3. N=8 solve: objective and certificate equivalence to the current native
   Hermitian SU(2) reference.
4. N=12 no-solve and solve: the quotient model must avoid the recorded
   20-minute primal timeout behavior.
5. L=20 no-solve, then solve only if the model-size gate is safe.
6. L=30 no-solve, then solve only if L=20 scaling and the live pressure gate
   justify it.

## Acceptance criteria

The Phase 3 formulation gap is closed only when all of the following hold:

- quotient and unquotiented small-N models are mathematically equivalent to
  the stated tolerances;
- all focused relaxation tests pass;
- package-level tests pass with only the existing intentional broken checks;
- N=12 quotient solve finishes with accepted MOSEK status and validated
  primal/dual/certificate residuals;
- guarded L=30 SU(2)+RDM+LSO+PSO solve finishes with accepted status;
- L=30 objective matches the reviewed symmetry-equivalent reference within
  `1e-6`;
- L=30 runtime is at most twice the reviewed QMBCertify runtime, or at least a
  twofold improvement over the recorded Phase 2 NCTSSoS runtime, as required
  by the active plan;
- peak RSS stays within the run-specific safety estimate and causes no swap;
- report/fixture provenance proves that the quotient formulation, not a
  smaller or weakened profile, produced the result;
- the generic fallback and all non-Pauli algebras remain unchanged.

If L=30 still fails, the result is not declared complete.  The new metrics
must identify the next structural bottleneck before any further formulation
change or large run.

## Implementation order

1. Add quotient-map unit tests and the immutable descriptor.
2. Implement deterministic pivot selection and reconstruction validation.
3. Rewrite individual `LinearMomentForm` values.
4. Rebuild complete `MomentLinearData` and preserve metadata.
5. Integrate the low-level translation/SU(2) direct builder.
6. Add public configuration and fail-closed recognizer checks.
7. Add report, harness, and fixture metrics.
8. Run the small-N equivalence ladder.
9. Run the load-gated N=12, L=20, and L=30 performance ladder.
10. Update the active parity plan with measured evidence and perform the full
    completion audit.
