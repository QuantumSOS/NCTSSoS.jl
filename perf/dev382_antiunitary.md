# DEV-382 Track A: Pauli antiunitary reduction

An antiunitary is a conjugate-linear norm-preserving map. This implementation
uses Θ = PᵣK, where K is conjugation in the computational basis and Pᵣ is a
Hermitian Pauli word. It is opt-in and does not change existing solver defaults.

## Contract and audit trail

A binary Pauli label is v = (x,z). Define q₀(v) = x·z modulo 2 and
ω(r,v) = rₓ·z + r_z·x modulo 2. The detector performs Gaussian elimination
over the binary field on ω(r,v) = q₀(v), including every nonzero support term
of the Hamiltonian, inequality constraints, equality constraints, and one-sided
moment equality constraints. Each input polynomial must be Hermitian. Thus
every input polynomial is fixed individually; permutation of constraints is
not accepted. An inconsistent system returns `nothing`, and attempting the
reduction then throws rather than silently changing the problem.

The detector sets free variables to zero and returns one solution. It does not
enumerate solutions or search for a requested square sign. An explicit
`PauliAntiunitary(x,z)` can select another solution. Both Θ² = +I and Θ² = −I
are supported and tested, including automatic detection of a forced −I case.
The −I case requires no invariant pure vector: averaging a state with its
antiunitary image preserves positivity, normalization, and the objective.
This is a moment relaxation, not a purity hierarchy.

Define qᵣ(v) = q₀(v) + ω(r,v) modulo 2. Odd moments, those with qᵣ(v) = 1,
are eliminated by exact substitution yᵥ = 0. For each original block, its
Hermitian Pauli row labels determine Dᵢᵢ = (−1)^qᵣ(vᵢ) and Uᵢᵢ = 1 or i
for even or odd qᵣ respectively. Invariance gives M = D conjugate(M) D.
The implementation emits U M U† directly as one real symmetric positive
semidefinite (PSD) block of side m, rather than the usual side-2m real lift.
It does not split parity sets, drop cross-parity entries, or remove odd rows.
The original clique and term-sparsity row labels remain attached to each block.

The returned `original` is the unreduced symbolic relaxation; `phases[b]`
records U for block b. Keeping both symbolic problems costs memory but makes
certificate auditing straightforward. The reduction currently constructs the
ordinary symbolic problem first; it avoids the doubled **solver** matrix, not
all intermediate symbolic storage. The reported moment count includes only
surviving moments appearing in moment matrices, not localizer-only moments.

## Usage and certificates

```julia
using NCTSSoS, JuMP
registry, (x,y,z) = create_pauli_variables(1:2)
pop = polyopt(1.0x[1] + 2z[1] + y[1]*y[2], registry)
config = SolverConfig(optimizer=nothing, order=1)
sparsity = compute_sparsity(pop, config)
reduction = antiunitary_moment_relaxation(pop, sparsity)
model, extract_moments = build_jump_model(reduction.moment_problem)
# Attach an optimizer, then optimize!(model).
# Missing odd physical moments have value zero by construction.
```

A Gram matrix is a PSD matrix multiplying a polynomial moment/localizing
matrix in a sum-of-squares certificate. Given a solved real Gram matrix Z,
`recover_antiunitary_gram(reduction, b, Z)` returns Q = U† Z U. Its convention
is tr(Q M) = Σᵢⱼ Qⱼᵢ Mᵢⱼ in the **original** Pauli algebra. For individually
fixed localizers, the odd coefficients cancel in this expression exactly.
Equality multipliers must also be included: the returned real scalar equality
forms are the real/imaginary parts after invariant averaging. Averaging is
valid because the input constraints are fixed and Pauli row bases are closed
under this diagonal antiunitary action. One-sided equalities express g|ψ⟩ = 0;
they must not be confused with two-sided operator equality localizers.

The regression test extracts Gram variables from the existing
`NCTSSoS.sos_dualize` model, checks their eigenvalues, recovers Q, and checks
every coefficient of the original polynomial identity, including dense and
sparse constrained problems with one-sided equalities. Tests also compare
primal and dual optima and actual JuMP PSD-cone dimensions.

Validation on this checkout: the final 64-check suite passed with Clarabel
(decomposition disabled) and Mosek. COSMO passed all 3,841 relaxation-suite
checks, including those 64 checks. Aqua passed 11 checks and ExplicitImports
completed successfully. The standalone expanded 61-check certificate suite
also passed with both COSMO and Mosek before the final three rejection guards
were added. The full repository test suite was not run.

COSMO's default chordal decomposition returned an incomplete matrix when
reading a primal PSD constraint's dual in the initial one-block experiment:
the recovered matrix had minimum eigenvalue −0.2704 although the coefficient
residual was small. A small residual alone is therefore not enough. Direct
SOS Gram extraction passes the PSD and coefficient checks with COSMO 0.8.11
and Mosek. Clarabel 0.11.1 required `chordal_decomposition_enable=false` to
avoid its `OrderedSet` PSD-completion error. These are solver-extraction
limitations, not reasons to omit PSD validation.

## Reproducible small benchmark

Run `perf/pauli_antiunitary.jl MODE 4 3` in an environment containing this
checkout, JuMP, and MosekTools. Modes are `complex`, `antiunitary`,
`complex-ts`, `antiunitary-ts`, and `specialized`. Use a fresh process for
each mode. The harness warms the exact workload once and takes three samples;
the raw samples and cold costs are in `results/dev382_antiunitary_mosek.txt`.
Measured modes ran sequentially, with no concurrent tests or example builds.

Environment: Linux x86-64 orb, Intel Xeon @ 2.60 GHz (`icelake-server`),
Julia 1.12.7, JuMP 1.32.1, Mosek.jl 11.2.1, MosekTools 0.15.10. Julia,
BLAS, and Mosek each use one thread. Mosek feasibility and relative-gap
tolerances are 1e-9; the license was verified with a trivial optimal solve.

Problem: periodic four-site Heisenberg Hamiltonian Σᵢ σᵢ·σᵢ₊₁/4, order two,
the same 49-word contiguous Pauli basis for every mode. Independent dense
spin-matrix diagonalization gives ground energy −2. `-ts` uses one `MF()`
term-sparsity step; this changes the relaxation strength. The specialized
baseline is `pauli_translation_invariant_moment_relaxation` with its existing
translation, sign, reflection, and conjugation reductions, not just this new
antiunitary reduction.

Construction includes sparsity, symbolic assembly, SOS lowering, and optimizer
attachment. Solver time is Mosek's reported time; full `optimize!` wall time
is also recorded in the raw data. Allocated MiB counts construction allocations,
not peak live memory. Peak RSS is the process high-water resident memory over
loading, compilation, warmup, solving, and certificate checks; it is **not**
isolated solver memory. Timings below are medians, with min–max ranges.

| Mode | Real cone sides | Build ms (range) | Allocated MiB | Peak RSS GiB | Solver ms (range) | Bound | Max coefficient residual |
|---|---|---:|---:|---:|---:|---:|---:|
| Complex baseline | 98 | 6.81 (6.75–7.54) | 7.589 | 1.762 | 62.14 (62.13–64.58) | −2.000000000558 | 3.16e-12 |
| Antiunitary | 49 | 5.56 (5.50–5.64) | 2.932 | 1.758 | 12.53 (12.44–12.56) | −2.000000000247 | 2.20e-12 |
| Complex, term-sparse | 16×8, 6×6 | 2.58 (2.55–2.86) | 1.234 | 1.771 | 6.02 (5.99–6.05) | −3.000000001332 | 3.33e-12 |
| Antiunitary, term-sparse | 16×4, 6×3 | 2.24 (2.22–3.09) | 0.641 | 1.765 | 3.36 (3.34–3.37) | −3.000000001147 | 5.96e-12 |
| Specialized real baseline | 1×4, 8×3, 3×2, 3×1 | 1.86 (1.83–1.93) | 0.619 | 1.725 | 1.84 (1.81–1.90) | −2.000000000055 | 6.57e-12 |

Each antiunitary cone side is halved with its row count unchanged: the dense
cone has 1,225 instead of 4,851 scalar triangular coordinates (74.7% fewer).
The dense solve is about 5× faster in this small experiment; the existing
specialized real construction remains substantially faster. Construction
allocations fall, but whole-process peak RSS is essentially unchanged and
does not establish a memory advantage. The sparse bound remains −3 rather
than −2: reduction preserves the selected relaxation, not dense strength.
All measured certificates have minimum eigenvalue above zero and primal–dual
gap below 2.2e-9. The specialized certificate is checked in its existing
translation-orbit coordinates; antiunitary certificates are pulled back to
the original Pauli polynomials. These are floating-point residual checks,
not interval-certified rigorous bounds or a large-system performance claim.

## Scope

This API accepts ordinary Pauli polynomial problems and existing Pauli-word
sparsity results. It does not compose with another symmetry-adapted basis,
enumerate antiunitaries, certify a pure state, implement a Jordan engine,
or add a state-polynomial/two-copy hierarchy. It inherits existing clique
filtering of custom word bases. Non-Hermitian equalities and constraints that
are invariant only as a set are deliberately rejected. Tests belong in CI;
the small reproducible benchmark and this audit note are developer material,
not an additional always-on performance test or duplicated tutorial workload.

## Documentation setup blocker

`make examples` stalled before executing examples: its log stopped changing
for about 20 minutes while Julia's resident memory grew to 24.6 GiB. A bounded
retry used `JULIA_PKG_PRECOMPILE_AUTO=0`, two precompile tasks, one BLAS/GC
thread, and `--heap-size-hint=4G`. After 6 minutes it still had emitted only
the command line and resident memory had risen to about 7 GiB, so it was
stopped under the delivery instruction to avoid renewed memory pressure.
Command/resource evidence is in `results/dev382_documentation_setup.txt`.
No generated documentation or example stamp was changed. Regeneration and
the subsequent `make servedocs` preview remain the sole delivery blocker;
the source-hash stamp is intentionally stale rather than fabricated.
