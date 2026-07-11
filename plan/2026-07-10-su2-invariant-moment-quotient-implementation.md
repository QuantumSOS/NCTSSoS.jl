# SU(2)-Invariant Moment Quotient Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace canonical non-singlet Pauli moment coordinates with a deterministic SU(2)-invariant quotient so the full L=30 translation/SU(2)+RDM+LSO+PSO relaxation meets the active solver runtime gate without weakening the model.

**Architecture:** Generate thin singlet-only Clebsch-Gordan rows through active support 8, choose canonical pivot moments as a coordinate chart of each singlet subspace, and rewrite finalized `MomentLinearData` through a fresh `MomentLinearBuilder`. The translation constructor activates the quotient only for a validated direct-linear SU(2) profile; existing JuMP/SOS lowering and the generic fallback remain unchanged.

**Tech Stack:** Julia 1.12, `LinearAlgebra`, existing NCTSSoS Pauli/SU(2) helpers, `MomentLinearBuilder`, JuMP, MOSEK for gated performance solves, Clarabel/COSMO for focused regression tests.

## Global Constraints

- All tests, profiles, examples, model builds, and solver runs execute through `easy-ssh`; never run Julia locally.
- Before each computational run, report a wall/RSS estimate; after it, report actual wall/RSS and whether the estimate matched.
- Every production behavior starts with a failing real-code test and a verified RED result.
- Preserve the existing generic fallback and all non-Pauli algebra behavior.
- The quotient is exact up to validated transform arithmetic; no numerical rank truncation or approximate feasible-set compression.
- Once a supported SU(2) fast path selects the quotient, failures are visible and never silently fall back.
- Large runs require fresh load/RSS telemetry, one Julia thread, one MOSEK thread, explicit estimates, and hard timeouts.
- Preserve all pre-existing dirty-worktree changes. Commit only isolated new files or explicitly staged hunks; never absorb unrelated existing modifications.
- Phase 5 2D symmetry and interval/rational certification remain separate projects.

---

## File Structure

- Create `src/optimization/pauli_su2_quotient.jl` — singlet recurrence, quotient maps, form rewriting, provenance wrappers, and `MomentLinearData` rebuilding.
- Create `test/relaxations/pauli_su2_quotient.jl` — focused unit, rejection, equivalence, and certificate tests.
- Modify `src/NCTSSoS.jl` — include the new source after `pauli_chains.jl`.
- Modify `src/optimization/pauli_chains.jl` — direct-builder keyword, quotient invocation, report fields/metrics, solve support, and public wrapper forwarding.
- Modify `src/optimization/interface.jl` — fast-path keyword validation/forwarding and user-facing option tests.
- Modify `perf/pauli_translation_profile.jl` — target/constructed quotient metrics.
- Modify `perf/pauli_translation_solver_probe.jl` — quotient selection, model-size fields, and fixture emission.
- Modify `test/relaxations/runtests.jl` — include the new focused test file.
- Modify `test/relaxations/interface.jl` — `SolverConfig`/`cs_nctssos` routing and rejection tests.
- Modify `test/expectations_loader.jl` — emitted fixture contract and deleted-run safety policy for quotient rows.
- Modify `test/data/expectations/heisenberg_qmbcertify_rdm.toml` only after real measured rows exist.
- Modify `plan/qmbcertify_parity.md` only after verification, recording measured evidence rather than predictions.

---

### Task 0: Establish the current-tree baseline

**Files:**
- Read only: current worktree and remote test output.

**Interfaces:**
- Produces: authoritative pre-change Git status and focused Pauli-chain test baseline.

- [x] **Step 1: Record worktree state without modifying it**

Run local read-only Git metadata commands:

```bash
git status --short --branch
git diff --stat
git rev-parse HEAD
```

Save the output in the execution notes. Do not clean, stash, reset, or stage pre-existing changes.

- [x] **Step 2: Check remote load and estimate the baseline**

Run `easy-ssh run "uptime; nproc; free -h"`. Use the prior focused-suite history to estimate 8–11 minutes and under 6 GiB RSS for `test/relaxations/pauli_chains.jl`.

- [x] **Step 3: Run the focused baseline remotely**

```bash
easy-ssh submit "timeout 900s env JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. --startup-file=no test/relaxations/pauli_chains.jl"
easy-ssh monitor
```

Expected: zero failures. If it fails, invoke `superpowers:systematic-debugging`, prove whether the failure predates quotient work, and do not proceed with an unexplained red baseline.

---

### Task 1: Thin singlet-only Clebsch-Gordan rows

**Files:**
- Create: `src/optimization/pauli_su2_quotient.jl`
- Create: `test/relaxations/pauli_su2_quotient.jl`
- Modify: `src/NCTSSoS.jl:132-134`
- Modify: `test/relaxations/runtests.jl`

**Interfaces:**
- Consumes: `_pauli_su2_word_local_spherical_transform()`, `_SparseTransformRows`, `_sparse_transform_rows`, `pauli_su2_word_blocks(s)`.
- Produces: `_pauli_su2_spin1_cg`, `_pauli_su2_singlet_spin_paths`, `_pauli_su2_word_singlet_rows`, and `_pauli_su2_word_singlet_multiplicity`.

- [x] **Step 1: Add the focused test file and failing multiplicity tests**

```julia
using Test
using LinearAlgebra
using NCTSSoS

@testset "SU(2)-invariant moment quotient" begin
    @testset "thin singlet rows" begin
        expected = [1, 0, 1, 1, 3, 6, 15, 36, 91]
        for support_size in 0:8
            rows = NCTSSoS._pauli_su2_word_singlet_rows(support_size)
            @test size(rows) == (expected[support_size + 1], 3^support_size)
            dense = NCTSSoS._dense_sparse_transform_rows(rows)
            @test dense * dense' ≈ I atol = 1e-11
        end
    end
end
```

Add `include("pauli_su2_quotient.jl")` to `test/relaxations/runtests.jl` and `include("optimization/pauli_su2_quotient.jl")` immediately after `pauli_chains.jl` in `src/NCTSSoS.jl`.

Create `src/optimization/pauli_su2_quotient.jl` with this initial content before adding the include, so package loading succeeds and the behavioral test—not a missing file—causes RED:

```julia
# SU(2)-invariant Pauli moment-coordinate quotient.
```

- [x] **Step 2: Run the focused test and verify RED**

Run remotely:

```bash
easy-ssh submit "timeout 180s env JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. --startup-file=no test/relaxations/pauli_su2_quotient.jl"
easy-ssh monitor
```

Expected: FAIL with `UndefVarError: _pauli_su2_word_singlet_rows not defined`.

- [x] **Step 3: Implement exact spin-1 CG coefficients and path enumeration**

Add to `src/optimization/pauli_su2_quotient.jl`:

```julia
const _PAULI_SU2_SINGLET_COUNTS = (1, 0, 1, 1, 3, 6, 15, 36, 91)

function _pauli_su2_factorial(n::Int)
    n >= 0 || throw(DomainError(n, "factorial argument must be non-negative"))
    return factorial(big(n))
end

function _pauli_su2_spin1_cg(
    parent_spin2::Int,
    parent_m2::Int,
    local_m2::Int,
    child_spin2::Int,
)
    all(iseven, (parent_spin2, parent_m2, local_m2, child_spin2)) ||
        throw(ArgumentError("spin-1 coupling requires even doubled quantum numbers"))
    local_m2 in (-2, 0, 2) || throw(DomainError(local_m2, "local m2 must be -2, 0, or 2"))

    j1 = parent_spin2 ÷ 2
    m1 = parent_m2 ÷ 2
    j2 = 1
    m2 = local_m2 ÷ 2
    j = child_spin2 ÷ 2
    m = m1 + m2
    abs(m1) <= j1 && abs(m2) <= j2 && abs(m) <= j || return 0.0
    abs(j1 - j2) <= j <= j1 + j2 || return 0.0

    triangle = (
        BigInt(2j + 1) * _pauli_su2_factorial(j + j1 - j2) *
        _pauli_su2_factorial(j - j1 + j2) *
        _pauli_su2_factorial(j1 + j2 - j)
    ) // _pauli_su2_factorial(j1 + j2 + j + 1)
    magnetic = (
        _pauli_su2_factorial(j + m) * _pauli_su2_factorial(j - m) *
        _pauli_su2_factorial(j1 - m1) * _pauli_su2_factorial(j1 + m1) *
        _pauli_su2_factorial(j2 - m2) * _pauli_su2_factorial(j2 + m2)
    ) // one(BigInt)

    kmin = maximum((0, j2 - j - m1, j1 + m2 - j))
    kmax = minimum((j1 + j2 - j, j1 - m1, j2 + m2))
    kmin <= kmax || return 0.0
    series = zero(Rational{BigInt})
    for k in kmin:kmax
        denom = _pauli_su2_factorial(k) *
            _pauli_su2_factorial(j1 + j2 - j - k) *
            _pauli_su2_factorial(j1 - m1 - k) *
            _pauli_su2_factorial(j2 + m2 - k) *
            _pauli_su2_factorial(j - j2 + m1 + k) *
            _pauli_su2_factorial(j - j1 - m2 + k)
        series += (isodd(k) ? -one(BigInt) : one(BigInt)) // denom
    end
    value = sqrt(Float64(triangle * magnetic)) * Float64(series)
    return iszero(value) ? 0.0 : value
end

function _pauli_su2_singlet_spin_paths(support_size::Int)
    support_size >= 0 || throw(DomainError(support_size, "support size must be non-negative"))
    paths = Vector{Vector{Int}}()
    function visit!(path::Vector{Int}, remaining::Int)
        current = last(path)
        if iszero(remaining)
            iszero(current) && push!(paths, copy(path))
            return
        end
        for child in abs(current - 2):2:(current + 2)
            child <= 2 * (remaining - 1) || continue
            push!(path, child)
            visit!(path, remaining - 1)
            pop!(path)
        end
    end
    visit!(Int[0], support_size)
    return paths
end
```

- [x] **Step 4: Implement thin Cartesian singlet rows**

```julia
function _pauli_su2_word_singlet_rows(support_size::Integer; atol::Real=1e-13)
    s = Int(support_size)
    s >= 0 || throw(DomainError(support_size, "support size must be non-negative"))
    s <= 8 || throw(ArgumentError("moment quotient supports active support through 8; got $s"))
    local_transform = _pauli_su2_word_local_spherical_transform()
    rows = Vector{Vector{Tuple{Int,ComplexF64}}}()
    for path in _pauli_su2_singlet_spin_paths(s)
        row = Tuple{Int,ComplexF64}[]
        for code0 in 0:(3^s - 1)
            amplitudes = Dict{Int,ComplexF64}(0 => 1.0 + 0.0im)
            code = code0
            for site in 1:s
                cartesian = mod(code, 3) + 1
                code = div(code, 3)
                parent_spin2 = path[site]
                child_spin2 = path[site + 1]
                next = Dict{Int,ComplexF64}()
                for (parent_m2, amplitude) in amplitudes
                    for (local_idx, local_m2) in enumerate((-2, 0, 2))
                        child_m2 = parent_m2 + local_m2
                        abs(child_m2) <= child_spin2 || continue
                        cg = _pauli_su2_spin1_cg(
                            parent_spin2,
                            parent_m2,
                            local_m2,
                            child_spin2,
                        )
                        iszero(cg) && continue
                        value = amplitude * cg * local_transform[local_idx, cartesian]
                        next[child_m2] = get(next, child_m2, 0.0 + 0.0im) + value
                    end
                end
                amplitudes = next
            end
            value = get(amplitudes, 0, 0.0 + 0.0im)
            abs(value) <= atol || push!(row, (code0 + 1, value))
        end
        if !isempty(row)
            phase = conj(first(row).second) / abs(first(row).second)
            row = [(column, coefficient * phase) for (column, coefficient) in row]
        end
        push!(rows, row)
    end
    length(rows) == _PAULI_SU2_SINGLET_COUNTS[s + 1] || throw(ArgumentError(
        "singlet path count mismatch for support $s",
    ))
    return _sparse_transform_rows(rows, 3^s)
end
```

Add a small `_dense_sparse_transform_rows` helper for validation and pivot selection; it materializes only `singlet_count × 3^s`, never a square transform.

- [x] **Step 5: Add raising-operator and dense-small-support cross-checks**

For supports 0:4, compare the row-space projectors of thin and existing dense singlet rows. For 5:8, construct the sparse total raising action from Cartesian rows via the local spherical transform and assert the maximum raised coefficient is at most `1e-10`. Assert the support-8 thin matrix has `91 * 6561` or fewer stored coefficients and that no helper returns a `6561 × 6561` matrix.

- [x] **Step 6: Run focused tests and verify GREEN**

Run the same focused command. Expected: all thin-row tests pass; wall under 180 seconds and peak RSS under 4 GiB.

- [x] **Step 7: Record an isolated checkpoint**

Commit the two new files and include lines only if they can be staged without pre-existing hunks. Otherwise record `git diff --check` plus the exact changed-file list and leave overlapping files unstaged.

---

### Task 2: Deterministic canonical pivot quotient maps

**Files:**
- Modify: `src/optimization/pauli_su2_quotient.jl`
- Modify: `test/relaxations/pauli_su2_quotient.jl`

**Interfaces:**
- Consumes: `MomentLinearData.key_to_monomial`, `_pauli_translation_support_orbit_word_patterns`, `_pauli_su2_translation_orbit_basis_lookup`, `_pauli_su2_translation_orbit_support_columns`, `key_lt`.
- Produces: `PauliSU2MomentOrbitQuotient`, `PauliSU2MomentQuotientDescriptor`, `_pauli_su2_moment_quotient_descriptor`.

- [x] **Step 1: Write failing pivot/reconstruction tests**

Build a support-complete translation-orbit basis for N=10, support 0:4. Assert:

```julia
descriptor = NCTSSoS._pauli_su2_moment_quotient_descriptor(linear, 10)
@test descriptor.raw_moment_count == length(linear.moments)
@test descriptor.quotient_moment_count == sum(length, getfield.(descriptor.orbits, :pivot_keys))
@test descriptor.quotient_moment_count < descriptor.raw_moment_count
@test descriptor.max_pivot_residual <= 1e-11
@test descriptor.max_reconstruction_residual <= 1e-10
@test isfinite(descriptor.max_condition)
```

Permute the input moment order and assert identical pivot keys and reconstruction coefficients. Add rejection tests for a missing axis pattern, a singular synthetic singlet matrix, and condition limit below the measured chart condition.

- [x] **Step 2: Run focused tests and verify RED**

Expected: `UndefVarError: _pauli_su2_moment_quotient_descriptor not defined`.

- [x] **Step 3: Add quotient data types**

```julia
struct PauliSU2MomentOrbitQuotient{K,C}
    support_orbit::Tuple{Vararg{Int}}
    source_keys::Vector{K}
    pivot_keys::Vector{K}
    singlet_labels::Vector{Any}
    reconstruction::Matrix{C}
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol
    pivot_residual::Float64
    invariant_residual::Float64
    reconstruction_residual::Float64
    condition::Float64
end

struct PauliSU2MomentQuotientDescriptor{K,C}
    n_sites::Int
    orbits::Vector{PauliSU2MomentOrbitQuotient{K,C}}
    source_to_orbit_row::Dict{K,Tuple{Int,Int}}
    raw_moment_count::Int
    quotient_moment_count::Int
    support_orbit_count::Int
    singlet_channel_support_counts::Vector{Pair{Int,Int}}
    max_pivot_residual::Float64
    max_invariant_residual::Float64
    max_reconstruction_residual::Float64
    max_condition::Float64
end
```

- [x] **Step 4: Implement deterministic pivot selection**

```julia
function _pauli_su2_select_pivot_columns(S::Matrix{ComplexF64}; condition_limit::Real)
    channel_count, source_count = size(S)
    iszero(channel_count) && return Int[], zeros(ComplexF64, source_count, 0), 1.0
    factor = qr(S, ColumnNorm())
    pivots = sort!(Int[factor.p[i] for i in 1:channel_count])
    Sstar = adjoint(S)
    chart = Sstar[pivots, :]
    chart_condition = cond(chart)
    isfinite(chart_condition) && chart_condition <= condition_limit || throw(ArgumentError(
        "SU(2) moment quotient coordinate chart condition $chart_condition exceeds $condition_limit",
    ))
    reconstruction = adjoint(adjoint(chart) \ adjoint(Sstar))
    return pivots, Matrix{ComplexF64}(reconstruction), Float64(chart_condition)
end
```

The implementation must not sort pivot indices until after selection if doing so would mismatch chart columns; when sorted for deterministic key order, apply the same permutation to `chart` and the reconstruction coordinates.

- [x] **Step 5: Build per-support-orbit maps**

Group `linear.key_to_monomial` by canonical translation support orbit. Use existing support-completeness validation and column lookup. Reuse `_pauli_su2_word_singlet_rows(s)` by active support size, choose pivots, map source keys to rows of `R`, and compute:

```julia
pivot_residual = maximum(abs, reconstruction[pivot_rows, :] - I)
projector = adjoint(S) * S
invariant_residual = maximum(abs, (I - projector) * reconstruction)
reconstruction_residual = maximum(abs, reconstruction * chart - adjoint(S))
```

Convert singlet rows to real coefficients only when the imaginary residual is below tolerance; otherwise keep complex data or fail when the source `MomentLinearData` coefficient type is real.

- [x] **Step 6: Run focused tests and verify GREEN**

Expected: deterministic pivot tests, support-completeness rejection, and residual checks pass.

---

### Task 3: Rewrite complete `MomentLinearData`

**Files:**
- Modify: `src/optimization/pauli_su2_quotient.jl`
- Modify: `test/relaxations/pauli_su2_quotient.jl`

**Interfaces:**
- Consumes: `MomentLinearBuilder`, `register_moment!`, `add_objective_terms!`, `add_psd_block!`, `_add_zero_constraint_trusted!`, `finalize!`.
- Produces: `_pauli_su2_rewrite_form`, `_pauli_su2_quotient_linear_data`, `PauliSU2QuotientBlockOrigin`, `PauliSU2QuotientTransform`.

- [x] **Step 1: Write failing form and linear-data equivalence tests**

For a small direct-linear N=6/order-2 SU(2) model, build the descriptor and quotient. Assert:

- every rewritten form evaluates identically on random singlet-coordinate samples;
- selected pivot moments reconstruct exactly;
- block count, cones, block sizes, labels, and logical row labels are unchanged;
- quotient moment count is smaller;
- Wigner/singlet zero rows that rewrite to zero are counted and removed;
- remaining zero rows retain concrete origin types and self-adjointness;
- `assert_moment_linear_data_invariants(quotient)` passes;
- mutating caller-owned key buffers after quotient construction changes nothing.

- [x] **Step 2: Run focused tests and verify RED**

Expected: `_pauli_su2_quotient_linear_data` is undefined.

- [x] **Step 3: Add provenance wrappers**

```julia
struct PauliSU2QuotientTransform
    family::Symbol
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol
    base::Any
    descriptor::Any
end

struct PauliSU2QuotientBlockOrigin <: BlockOrigin
    label::Any
    logical_row_labels::Vector{Any}
    transform::PauliSU2QuotientTransform
    source_origin::BlockOrigin
end
```

Construct the wrapper by forwarding `label`, `logical_row_labels`, and the base `transform` through `hasproperty`. This keeps `translation_linear_provenance` and SOS diagnostics working without special cases.

- [x] **Step 4: Implement sparse form rewriting**

```julia
function _pauli_su2_rewrite_form(
    form::LinearMomentForm{K,C},
    descriptor::PauliSU2MomentQuotientDescriptor{K,QC},
    ::Type{C};
    atol::Real,
) where {K,C,QC}
    terms = Pair{K,C}[]
    for (source_key, source_coef) in form
        orbit_idx, source_row = descriptor.source_to_orbit_row[source_key]
        orbit = descriptor.orbits[orbit_idx]
        for pivot_col in eachindex(orbit.pivot_keys)
            coefficient = source_coef * orbit.reconstruction[source_row, pivot_col]
            abs(coefficient) <= atol && continue
            converted = convert(C, coefficient)
            push!(terms, orbit.pivot_keys[pivot_col] => converted)
        end
    end
    return _linear_moment_form_from_owned_pairs!(terms)
end
```

For real `C`, reject an imaginary coefficient above tolerance and use `real(coefficient)` otherwise.

- [x] **Step 5: Rebuild immutable linear data**

Create a fresh `MomentLinearBuilder(K, C, M)`. Register every pivot key with its original monomial. Rewrite and add the objective, each PSD block with wrapped metadata, and each nonzero scalar row with its original origin/trust flag. Track eliminated zero rows by `_translation_zero_origin_histogram_key`. Finalize and return:

```julia
(
    linear=finalize!(builder; stage_times_ns, stage_prefix=:su2_moment_quotient),
    descriptor=descriptor,
    eliminated_zero_row_count=eliminated_count,
    eliminated_zero_feature_histogram=histogram,
)
```

- [x] **Step 6: Verify RED-GREEN and mutation safety**

Run the focused file. Expected: all equivalence/invariant tests pass and no mocks are used.

---

### Task 4: Integrate the direct translation/SU(2) constructor and reports

**Files:**
- Modify: `src/optimization/pauli_chains.jl`
- Modify: `src/optimization/pauli_su2_quotient.jl`
- Modify: `test/relaxations/pauli_su2_quotient.jl`
- Modify: `test/relaxations/pauli_chains.jl`

**Interfaces:**
- Consumes: `_pauli_translation_base_linear_relaxation`, `TranslationInvariantReport`, `translation_report_metrics`, `translation_solve_support`.
- Produces: `su2_moment_quotient` keyword and report metrics.

- [x] **Step 1: Write failing integration and rejection tests**

Add N=6/order-2 direct-linear cases with quotient disabled/enabled for base SU(2), reflection, SU(2) RDM extend, LSO, PSO, moment equalities, and singlet equalities. Assert equal block shapes and small-solve objectives, lower moment counts, quotient metrics, and fail-closed errors for:

- `su2_moment_quotient=true` without `base_su2_extend_rdm=true`;
- quotient without `su2_symmetry=true`;
- quotient on symbolic or QMBCertify-base routes;
- quotient with a non-invariant objective when checks are enabled.

- [x] **Step 2: Run focused test and verify RED**

Expected: unsupported keyword or missing report fields.

- [x] **Step 3: Extend `TranslationInvariantReport` with defaulted quotient fields**

Append fields:

```julia
su2_moment_quotient::Bool
su2_moment_raw_count::Int
su2_moment_quotient_count::Int
su2_moment_quotient_reduction_ratio::Float64
su2_moment_support_orbit_count::Int
su2_moment_singlet_channel_support_counts::Vector{Pair{Int,Int}}
su2_moment_max_pivot_residual::Float64
su2_moment_max_invariant_residual::Float64
su2_moment_max_reconstruction_residual::Float64
su2_moment_max_condition::Float64
su2_moment_eliminated_zero_row_count::Int
su2_moment_eliminated_zero_feature_histogram::Vector{Pair{Any,Int}}
```

Update every `TranslationInvariantReport(...)` constructor. Non-quotient paths use `false`, identical raw/quotient counts, ratio `1.0`, zero residuals, and empty histograms.

- [x] **Step 4: Invoke the quotient after direct builder finalization**

Add `su2_moment_quotient::Bool=false`, `su2_moment_quotient_atol::Real=1e-11`, and `su2_moment_quotient_condition_limit::Real=1e10` to `_pauli_translation_base_linear_relaxation`. Validate the option combination. Immediately after raw `finalize!`:

```julia
raw_linear = finalize!(builder)
quotient_result = su2_moment_quotient ?
    _pauli_su2_quotient_linear_data(
        raw_linear,
        n;
        atol=su2_moment_quotient_atol,
        condition_limit=su2_moment_quotient_condition_limit,
        stage_times_ns=construction_stage_times_ns,
    ) : nothing
linear = isnothing(quotient_result) ? raw_linear : quotient_result.linear
```

Populate report fields from `raw_linear` and the descriptor.

- [x] **Step 5: Extend metrics and solve support**

Expose all quotient fields in `translation_report_metrics`. Compute SOS coefficient equations from the actual post-quotient `linear_moment_count`. `translation_solve_support` accepts the quotient only when its residual fields are finite/below tolerance and the quotient count is positive; malformed reports return blocker `:su2_moment_quotient`.

- [x] **Step 6: Run focused and Pauli-chain tests**

Run the new focused file first, then `test/relaxations/pauli_chains.jl`. Expected: quotient integration green; existing tests unchanged except intentional updated report fields.

---

### Task 5: Public routing, profile harness, solver probe, and fixtures

**Files:**
- Modify: `src/optimization/pauli_chains.jl`
- Modify: `src/optimization/interface.jl`
- Modify: `perf/pauli_translation_profile.jl`
- Modify: `perf/pauli_translation_solver_probe.jl`
- Modify: `test/relaxations/interface.jl`
- Modify: `test/expectations_loader.jl`
- Modify: `test/data/expectations/heisenberg_qmbcertify_rdm.toml`

**Interfaces:**
- Produces: public keyword forwarding, `NCTS_TRANSLATION_SU2_MOMENT_QUOTIENT`, `NCTS_SOLVER_PROBE_SU2_MOMENT_QUOTIENT`, and persisted quotient evidence.

- [x] **Step 1: Write failing public routing tests**

Through both `pauli_translation_invariant_nctssos(pop, config; ...)` and `cs_nctssos(pop, config; ...)`, assert the quotient-enabled direct SU(2) route returns `MomentLinearData`, reports a smaller coordinate count, and preserves objective/certificate residual. Add rejection tests for unsupported option combinations.

- [x] **Step 2: Verify RED**

Run `test/relaxations/interface.jl`; expected failure is an unsupported/unforwarded `su2_moment_quotient` keyword.

- [x] **Step 3: Forward the public option**

Add keywords to `pauli_translation_invariant_nctssos` and pass them only to the direct-linear route. The `SolverConfig` bridge already forwards keyword arguments; update its supported-key validation and error text, not the generic `SolverConfig` struct.

- [x] **Step 4: Add harness flags and printed metrics**

Profile harness:

```julia
su2_moment_quotient = _parse_bool_env(
    "NCTS_TRANSLATION_SU2_MOMENT_QUOTIENT",
    base_su2_extend_rdm,
)
```

Solver probe:

```julia
probe_su2_moment_quotient = _parse_bool_env(
    "NCTS_SOLVER_PROBE_SU2_MOMENT_QUOTIENT",
    base_su2_extend_rdm,
)
```

Print raw/quotient counts, ratio, support/channel counts, residuals, condition, eliminated rows, actual scalar equation count, dense-Schur proxy, and the exact flag value. Update usage comments with estimates and telemetry requirements.

- [x] **Step 5: Add fixture-emitter and loader contracts**

The synthetic loader test uses internally consistent concrete values:

```toml
su2_moment_quotient = true
su2_moment_raw_count = 100
su2_moment_quotient_count = 10
su2_moment_max_reconstruction_residual = 1.0e-12
su2_moment_max_condition = 5.0
```

The synthetic row lives only in the test body. Do not insert predicted numeric values into the scientific TOML fixture. Add loader tests requiring positive reduction, finite residuals/condition, consistent actual equation counts, and execution-state metadata. Promote only real emitted rows after runs complete.

- [x] **Step 6: Run interface and loader tests**

Run `test/relaxations/interface.jl` and `test/expectations_loader.jl` through `easy-ssh` with one Julia thread and hard timeouts. Expected: all routing and fixture contracts pass.

---

### Task 6: Small-N equivalence and certificate gate

**Files:**
- Modify: `test/relaxations/pauli_su2_quotient.jl`
- Modify: `test/relaxations/pauli_chains.jl`
- Modify: `plan/qmbcertify_parity.md` only after evidence exists.

- [x] **Step 1: Run the complete small-N matrix**

Run quotient disabled/enabled at valid `N > 2d` for base, reflection, closed/extended RDM, LSO, PSO, moment equality, singlet equality, and combined profiles. Use native Hermitian SOS where supported and Clarabel/COSMO for cheap stable checks.

- [x] **Step 2: Validate solver quality**

For every solved case record termination status, primal feasibility, native Hermitian minimum eigenvalue, objective delta, dual feasibility/gap when available, and `sos_dual_certificate_residual`. Required objective delta: `<= 1e-8` for stable small cases; certificate residual: existing suite tolerance or tighter.

- [x] **Step 3: Run touched relaxation suites**

Run:

```bash
test/relaxations/pauli_su2_quotient.jl
test/relaxations/pauli_chains.jl
test/relaxations/moment_linear.jl
test/relaxations/lowering.jl
test/relaxations/sos.jl
test/relaxations/interface.jl
test/expectations_loader.jl
```

Each uses `easy-ssh submit` plus `monitor`, one Julia thread, and a hard timeout based on recorded suite history.

- [x] **Step 4: Record evidence in the active plan**

Add measured commands, wall time, peak RSS, test counts, objective deltas, residuals, and quotient counts. Do not write completion language yet.

---

### Task 7: Load-gated performance ladder

**Files:**
- Modify: `test/data/expectations/heisenberg_qmbcertify_rdm.toml` with real emitted rows.
- Modify: `plan/qmbcertify_parity.md` with measured evidence.

- [x] **Step 1: N=10 target-only and constructed no-solver**

Verify analytic quotient counts, then construct the combined SU(2)+k=8 RDM+LSO7+PSO3 profile. Record raw/quotient counts, build stages, scalar equations, dense-Schur proxy, wall, and RSS.

- [x] **Step 2: N=8/N=10 small MOSEK solve**

Use the docs MOSEK environment, native Hermitian SOS, one thread, explicit tolerances, and a 10-minute guard. Require accepted status, objective parity, and certificate residual.

- [x] **Step 3: N=12 no-solve and solve**

No-solve must show a safe model-size gate. The solve must finish rather than repeat the recorded 20-minute primal timeout. Record full quality metrics.

- [x] **Step 4: L=20 no-solve and solve**

Launch only after fresh telemetry and an estimate derived from N=12. If model-size or pressure gates reject the run, fix/downsize the formulation rather than override the gate.

- [x] **Step 5: L=30 no-solve and solve**

Launch only if L=20 scaling predicts the run stays within safety limits. Acceptance requires:

- accepted MOSEK status;
- objective within `1e-6` of the reviewed symmetry-equivalent reference;
- runtime `<= 2 ×` reviewed QMBCertify or `>= 2 ×` faster than recorded Phase 2 NCTSSoS;
- no swap and RSS within estimate;
- emitted fixture proves quotient-enabled full profile.

- [x] **Step 6: Promote real fixture rows and rerun loader**

Copy the harness-emitted TOML exactly, review it, add it to the expectation fixture, and rerun the full loader.

---

### Task 8: Full verification and completion audit

**Files:**
- Modify: `plan/qmbcertify_parity.md`

- [x] **Step 1: Run package hygiene checks relevant to touched code**

Run formatting check, ExplicitImports/Aqua/JET quality suites already wired by the package, and `git diff --check`. Do not run local Julia.

- [x] **Step 2: Run full package tests**

After fresh telemetry, estimate 24–35 minutes and under 7 GiB from the prior baseline, then run `Pkg.test()` through `easy-ssh` with a hard timeout. Require zero failures and only the existing intentional broken checks.

- [x] **Step 3: Requirement-by-requirement audit**

For every design acceptance criterion, cite authoritative evidence: source line, focused test, full test, harness row, solver status, objective delta, residual, runtime, RSS, and fixture contract. Any missing or indirect item remains unfinished.

- [x] **Step 4: Update the active plan**

Replace stale “remaining formulation gap” language with measured completion evidence only if every gate passed. Keep Phase 5 and interval certification explicitly separate.

#### Completion audit

- **Exact quotient construction:** thin Clebsch-Gordan singlet rows, exact sparse-integer
  projected-rank computation, deterministic conditioned pivot charts, and the quotient
  descriptor are implemented in `src/optimization/pauli_su2_quotient.jl`. Direct/custom
  inputs require complete coordinate support; only the generated translation-relaxation
  path may use an explicitly recorded generated projection. Descriptor and per-orbit
  provenance record that distinction. Complete sparse rewriting covers every
  `MomentLinearData` component, including rewriting labeled zero rows before deciding
  whether they are redundant. The focused regression suite passed all 2,804 checks in
  85.1 s, including exact ranks through support eight, missing-pattern rejection,
  generated-mask provenance, singular-chart rejection, row permutation, mislabeled
  residual rows, and four-thread cache coverage.
- **No equality-family fallback:** quotient routing is integrated in
  `src/optimization/pauli_chains.jl`; cone-only Wigner-Eckart construction is controlled
  by `emit_invariance_rows`, and generated-coordinate projection is enabled only when
  the caller did not supply a custom basis. Focused tests cover enabled/disabled,
  reflection, RDM, LSO, PSO, moment-equality, singlet-equality, and combined profiles.
- **Solver equivalence and conditioning:** small native-Hermitian solves meet the
  `1e-8` objective-delta gate and the suite certificate tolerance. Capped max-absolute
  coefficient scaling is exercised at `test/relaxations/sos.jl:137-143`; the default
  remains unscaled.
- **Large target:** fixture
  `NCTSSOS_A2_SU2_MOMENT_QUOTIENT_L30_source_like_sos_dual_mosek_2026_07_11` is stored
  in `test/data/expectations/heisenberg_qmbcertify_rdm.toml` and is labeled large-N
  NCTSSoS backend evidence rather than a reviewed QMBCertify parity run. It records strict
  `OPTIMAL`, a `1.4046251926234205e-8` relative gap, `3.3181295178152936e-7`
  certificate residual, `5.96901113e-7` same-profile objective delta, 821.14 s remote
  wall, 12,801,622,016-byte peak RSS, and no swap.
- **Post-review large-model recheck:** the final exact-rank implementation rebuilt and
  dualized the L=30 model without solving: raw moments `94_129`, quotient moments `3_529`,
  185 cones with maximum block 28, `12_001` dual variables, and `4_528` scalar
  equalities. Construction took 420.5785 s, the remote run took 8:48.20, peak RSS was
  13,071,404 KiB, and no swap occurred. Those counts exactly match the accepted solved
  fixture, so the accepted strict-`OPTIMAL` solution remains the authoritative solve.
- **Independent review:** two fresh-context reviews were completed. The first found
  incomplete-support, rewrite-order, adversarial-test, cache-synchronization, and fixture
  labeling defects; all were reproduced and fixed. The final review reported no
  actionable findings.
- **Package-wide regression gate:** the final remote `Pkg.test()` run passed 17,378
  tests with only the two existing intentional broken checks in 26:08.0 test time
  (26:31.28 wrapper wall), at 6,287,336 KiB maximum RSS and no swap. Aqua,
  ExplicitImports, doctests, all touched suites, loader validation, and `git diff
  --check` are green.
- **Scope boundary:** Phase 3 is complete. Phase 4 interval certification and Phase 5
  two-dimensional lattice work remain explicitly separate projects and are not part of
  this implementation claim.

- [x] **Step 5: Finish the development branch**

Invoke `superpowers:finishing-a-development-branch`, re-verify, and present integration options. Do not merge, push, or open a PR without the user’s choice.
