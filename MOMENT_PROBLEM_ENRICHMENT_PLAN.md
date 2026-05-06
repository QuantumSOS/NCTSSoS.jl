# MomentProblem enrichment plan: linear-form cache + clean lowering

**Status**: proposal, not yet started.
**Supersedes**: `LOWERING_REFACTOR_PLAN.md` (steps 1–4 already shipped in `src/optimization/lowering.jl`) and `MOMENT_SOS_PIPELINE_ANALYSIS.md` (analysis material, key conclusions folded in here).

---

## TL;DR

The lowering layer (`src/optimization/lowering.jl`) currently re-derives linear-form structure (pivots, phases, orphans, identity, key universe) on every build by walking polynomial matrices at runtime. This wastes work, requires seven open contracts (C1–C7 in the old plan), has subtle bugs (Hermitian pivot injectivity, `/ phase` type-widening, missing provenance for Phase 2), and forces every consumer (lowering, SOS dual, demos) to re-import algebra-specific code (`simplify`, `symmetric_canon`).

**Fix**: in `moment_relax`, after all symbolic mutations are finished, compute and cache a `MomentLinearData` view onto `MomentProblem`. Lowering and SOS dual then read this cache and become mechanical.

This is a precondition for Phase 2's $\mathcal{A}\mathcal{A}^*$ diagnostic getting block-level provenance and for Phase 3 production at H₄/Nk=2.

---

## What is in scope

- Add a `MomentLinearData{K,C,M}` field to `MomentProblem` carrying linear-form views, pivots, key universe, identity, adjoint map, and per-block provenance.
- Populate it at the end of `moment_relax`, after `_add_parity_constraints!` and `_add_moment_eq_constraints!` have finished mutating `mp.constraints`.
- Run a pivot-coverage audit on real problems (H₂/Nk=2, small Pauli/Fermionic/Bosonic) before declaring the design final.
- Rewrite `src/optimization/lowering.jl` to consume the cache. Drop runtime `discover_pivots`, `_all_moment_keys`, `_pivot_candidate`, `_pivot_entry_lookup`. Keep the public `build_jump_model` signature.
- Fix the four oracle/jinguo-identified bugs:
  1. Hermitian upper-triangle pivot injectivity (add `adjoint::Bool`).
  2. `conj(phase) * X` instead of `X / phase`.
  3. `JuMP.add_to_expression!` instead of `sum(...)` for deterministic accumulation.
  4. Typed `BlockOrigin` instead of `Tuple{Symbol, Int, Int}`.
- Add `ScalarLinearConstraint{K,C}` so zero/binding/normalization constraints carry provenance for Phase 2 and the SOS dual.
- Migrate `sos_dualize` to consume the cache. The factor-of-2 Hermitian dual scaling is preserved by an explicit helper, not auto-collapsed.
- Add an HAI smoke test (`test/relaxations/h2_nk2_bridge_smoke.jl`) checking the cone-count predicate from the old plan.

## What is out of scope

- Solving H₄/Nk=2.
- Symmetry block-diagonalization.
- BPSDP tuning (warm start, ALM schedule, rank).
- Facial reduction / Mazziotti-CG integration.
- `:psd_blocks, representation = :real` (real-lift PSD-first). Add only if a non-Hermitian-capable solver demands it.
- `StateMomentProblem` enrichment. Same pattern can apply later but is a separate PR.
- `SOSDualForm` symbolic dual + reusable SOS lowering. Tracked separately; not blocked on this plan.

---

## Background: what the codebase has today

### Current `MomentProblem` (`src/optimization/moment.jl:46`)

```julia
struct MomentProblem{A<:AlgebraType, T<:Integer, M<:NormalMonomial{A,T}, P<:Polynomial{A,T}}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}   # :PSD | :HPSD | :Zero
    total_basis::Vector{M}
    n_unique_moment_matrix_elements::Int
end
```

`moment_relax` builds polynomial matrices, then *mutates* `mp.constraints` via:
- `_add_parity_constraints!(mp)` (fermionic parity superselection)
- `_add_moment_eq_constraints!(mp, pop, ...)` (one-sided localizing rows)

Any cache must be populated **after** these mutations.

### Current `lowering.jl` (~587 LOC)

Already implements:
- `Pivot{key::Any, constraint_idx::Int, block_idx::Int, i::Int, j::Int, phase::ComplexF64, cone::Symbol}` (no `adjoint` flag)
- `AffineResolver`, `PivotResolver`, `BlockMomentResolver`
- `discover_pivots(mp; orphan_policy)` — runtime O(entries) pivot search
- `build_jump_model(mp; formulation, representation, orphan_policy)` with three orphan policies (`:error`, `:free_variables`, `:aux_psd_free`)
- Resolver uses `block[i,j] / phase` (the type-widening bug)
- Resolver uses `expr += coef * resolver(α)` (non-deterministic accumulation)

### Current `sos.jl` (~470 LOC)

`sos_dualize(mp::MomentProblem)` and `_sos_dualize_hermitian` walk polynomial matrices term-by-term, build coefficient-matching equations per canonical key, and emit lifted real `2n × 2n` PSD duals for Hermitian primal cones. The factor-of-2 scaling is intentional and tested.

### The H₂/Nk=2 pathology

From `output/phase2/h2_nk2_bpsdp_bridge_inspect.json`:

```
symbolic complex moments        9393
current direct real variables   18786
JuMP variables before bridge    18786
BPSDP blocks after bridge       406724
BPSDP scalar 1x1 blocks         406706
BPSDP non-scalar PSD blocks     18
BPSDP A shape                   333405 x 479158
BPSDP A nnz                     1279400
```

Free-moment variable lowering creates 18,786 real moment variables. After MOI bridges, 406,706 scalar 1×1 PSD cones. BPSDP cannot solve this. The fix is `formulation=:psd_blocks, representation=:complex`.

---

## The new types

```julia
# ── linear forms ────────────────────────────────────────────────────

"""
    LinearMomentForm{K,C}

Sorted, deduplicated, zero-pruned linear combination of canonical moments.
`terms[i] = (key, coefficient)` with keys in ascending `key_lt` order.
Iteration order is deterministic across runs and Julia versions.
"""
struct LinearMomentForm{K,C}
    terms::Vector{Pair{K,C}}
end

# ── pivots ──────────────────────────────────────────────────────────

"""
    Pivot{C}

Records that physical canonical moment `α` is represented by entry
`X[block][row,col]` with phase `phase` (a unit complex number with
`abs(phase) == 1`). The relation is

    ⟨α⟩ = conj(phase) * X[block][row,col]                    if !adjoint
    ⟨α†⟩ = conj(phase) * X[block][row,col]                   if adjoint

For Hermitian blocks under upper-triangle storage, a moment and its
adjoint may share a position with the conjugate phase; `adjoint`
disambiguates which one this `Pivot` declares.
"""
struct Pivot{C}
    block::Int
    row::Int
    col::Int
    phase::C
    adjoint::Bool
end

# ── block provenance ────────────────────────────────────────────────

abstract type BlockOrigin end

struct MomentMatrixOrigin <: BlockOrigin
    clique::Int
    ts_block::Int
end

struct LocalizingOrigin <: BlockOrigin
    clique::Int
    cons_idx::Int          # index into PolyOpt.eq/ineq concatenated list
    ts_block::Int
end

struct GlobalOrigin <: BlockOrigin
    cons_idx::Int
end

struct AuxOrigin <: BlockOrigin
    reason::Symbol         # :orphan_packing
    key_group::Vector      # canonical keys this aux block serves
end

struct BlockMeta{M}
    cone::Symbol           # :PSD | :HPSD | :AuxHPSD
    origin::BlockOrigin
    row_labels::Vector{M}  # row_labels[i] = the basis monomial used at row i
end

struct PSDBlockLin{K,C,M}
    size::Int
    entries::Matrix{LinearMomentForm{K,C}}
    meta::BlockMeta{M}
end

# ── scalar constraints (zeros / bindings / normalization) ───────────

abstract type ConstraintOrigin end

struct ZeroMatrixOrigin <: ConstraintOrigin
    constraint_idx::Int    # index into mp.constraints
    row::Int
    col::Int
    part::Symbol           # :real | :imag | :scalar
end

struct ParityOrigin <: ConstraintOrigin
    cons_idx::Int          # which fermionic parity row
end

struct MomentEqOrigin <: ConstraintOrigin
    cons_idx::Int          # which moment_eq_constraint
    row::Int
end

struct IdentityOrigin <: ConstraintOrigin end

struct ScalarLinearConstraint{K,C}
    form::LinearMomentForm{K,C}
    kind::Symbol                # :zero | :identity_norm
    origin::ConstraintOrigin
end

# ── linear data cache on MomentProblem ──────────────────────────────

"""
    MomentLinearData{K,C,M}

Derived linear-form view of a `MomentProblem`. Populated by `moment_relax`
after all symbolic constraint mutations are complete. Treated as a cache:
once attached, `mp.constraints`, `mp.objective`, and `mp.total_basis` must
not be mutated.
"""
struct MomentLinearData{K,C,M}
    moments::Vector{K}                                     # sorted, unique
    moment_index::Dict{K,Int}                              # α → position in moments
    identity::K                                            # canonical key for empty word
    key_to_monomial::Dict{K,M}                             # α → representative monomial
    adjoint_key::Dict{K,K}                                 # α → ᾱ (identity for real algebras)

    psd_blocks_lin::Vector{PSDBlockLin{K,C,M}}             # mirror of :PSD/:HPSD
    psd_block_constraint_idx::Vector{Int}                  # block b ↔ mp.constraints[idx]
    zero_constraints::Vector{ScalarLinearConstraint{K,C}}  # already split for complex
    objective_lin::LinearMomentForm{K,C}

    pivots::Dict{K, Pivot{C}}                              # primary pivot per key
    pivot_at::Dict{Tuple{Int,Int,Int}, Vector{K}}          # (b,i,j) → keys (≤ 2)
    free_keys::Vector{K}                                   # keys with no pivot (orphans)
end

# ── enriched MomentProblem ──────────────────────────────────────────

struct MomentProblem{A<:AlgebraType, T<:Integer, C, M<:NormalMonomial{A,T}, P<:Polynomial{A,T}, K}
    # legacy fields, untouched
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}
    total_basis::Vector{M}
    n_unique_moment_matrix_elements::Int

    # new cache field
    linear::MomentLinearData{K,C,M}
end
```

**Type-parameter contract**: `K` is the canonical-key type used by
`symmetric_canon(expval(::M))`. In current code `K = Vector{T}`. Code that
needs to be `K`-generic must use `key_lt`/`key_isequal` traits, not
`Base.<` on raw vectors.

**Migration shim for `MomentProblem{A,T,M,P}` callers**: provide an inner
constructor `MomentProblem{A,T,M,P}(...)` that infers `C` from `P` and `K`
from a representative `total_basis` entry, builds `MomentLinearData`, and
returns `MomentProblem{A,T,C,M,P,K}`. Existing field accesses
(`mp.objective`, `mp.constraints`, etc.) continue to work.

---

## Constructor invariants

`moment_relax` must enforce these on the returned `MomentProblem`. Violation is a programmer error and raises during construction, not at lowering time.

| # | Invariant | Enforcement |
|---|---|---|
| I1 | `Set(linear.moments) == keys(linear.key_to_monomial)` | Constructor assertion |
| I2 | `linear.moments[linear.moment_index[α]] == α` for all α | Constructor assertion |
| I3 | `symmetric_canon(expval(linear.key_to_monomial[α])) == α` | Constructor assertion |
| I4 | `linear.identity == symmetric_canon(expval(one(M)))` | Constructor assertion |
| I5 | Every key in any `LinearMomentForm.terms` is in `linear.moments` | Constructor assertion |
| I6 | `LinearMomentForm.terms` is sorted by `key_lt`, no duplicate keys, no zero coefficients | Smart constructor only path |
| I7 | `Set(keys(pivots)) ∪ Set(free_keys) == Set(moments)` | Constructor assertion |
| I8 | `Set(keys(pivots)) ∩ Set(free_keys) == ∅` | Constructor assertion |
| I9 | `pivot_at[(b,i,j)]` lists exactly the keys whose pivot points to `(b,i,j)` | Constructor assertion |
| I10 | For each `(b,i,j)`, `length(pivot_at[(b,i,j)]) ≤ 2`; if `==2`, the two keys are adjoints | Constructor assertion |
| I11 | `pivots[α].phase` is a unit complex (`abs(phase) == 1`); for real algebras, `phase == 1` | Constructor assertion |
| I12 | For complex algebras, `adjoint_key[adjoint_key[α]] == α` | Constructor assertion |
| I13 | `zero_constraints` contains only real linear forms over the moment universe (Hermitian/anti-Hermitian split done) | Constructor by construction |
| I14 | `psd_blocks_lin[b].size == size(mp.constraints[psd_block_constraint_idx[b]][2], 1)` | Constructor assertion |

Smart constructor for `LinearMomentForm`:

```julia
function LinearMomentForm{K,C}(pairs) where {K,C}
    terms = Vector{Pair{K,C}}()
    for (k, c) in pairs
        iszero(c) && continue
        push!(terms, k => convert(C, c))
    end
    sort!(terms; by = x -> x.first, lt = key_lt)
    # Combine duplicate keys
    out = Vector{Pair{K,C}}()
    for (k, c) in terms
        if !isempty(out) && key_isequal(out[end].first, k)
            out[end] = out[end].first => (out[end].second + c)
            iszero(out[end].second) && pop!(out)
        else
            push!(out, k => c)
        end
    end
    return LinearMomentForm{K,C}(out)
end
```

---

## Algorithm: populating `MomentLinearData` in `moment_relax`

**Position in `moment_relax`**: after `_add_parity_constraints!` and `_add_moment_eq_constraints!`, before returning `mp`.

```
1.  Walk every entry of mp.constraints and mp.objective. For each polynomial term,
    compute symmetric_canon(expval(mono)) and accumulate into a temporary
    Dict{K, NormalMonomial{A,T}}: representative monomial per canonical key.
    This builds key_to_monomial.

2.  moments := sort(collect(keys(key_to_monomial)); lt = key_lt)
    moment_index := Dict(α => i for (i, α) in enumerate(moments))
    identity := symmetric_canon(expval(one(M)))

3.  adjoint_key:
    for each α in moments:
        m := key_to_monomial[α]
        α_adj := symmetric_canon(expval(adjoint(m)))
        adjoint_key[α] := α_adj
    For real algebras (where _is_complex_problem(A) == false), adjoint_key is
    the identity map.

4.  Linearize objective:
    objective_lin := linearize(mp.objective, moment_index, key_to_monomial)
    where linearize walks poly.terms, applies symmetric_canon, and produces
    a LinearMomentForm via the smart constructor.

5.  Linearize PSD/HPSD blocks. For each (idx, (cone, mat)) in pairs(mp.constraints)
    where cone in (:PSD, :HPSD):
        block_idx := length(psd_blocks_lin) + 1
        entries[i,j] := linearize(mat[i,j], ...) for all i, j
        meta := BlockMeta(cone, infer_origin(idx, mp), row_labels)
        push!(psd_blocks_lin, PSDBlockLin(size(mat,1), entries, meta))
        push!(psd_block_constraint_idx, idx)

    infer_origin uses the order of construction in moment_relax: clique k's
    moment matrix is the first PSD block from clique k; localizing matrices
    follow in cons-index order; global constraints come last; aux blocks
    (added in step 8) are AuxOrigin.

    row_labels[i] := basis_monomial used to build row i of the block. The
    builder must thread this from _build_constraint_matrix's local_basis.

6.  Linearize zero constraints. For each (idx, (:Zero, mat)) in pairs(mp.constraints):
        For each (i, j) in upper triangle of mat:
            iszero(mat[i,j]) && continue
            f := linearize(mat[i,j], ...)
            assert imag-component is provably zero (already split in _zero_constraint_components)
            push!(zero_constraints, ScalarLinearConstraint(
                form = f,
                kind = :zero,
                origin = ZeroMatrixOrigin(idx, i, j, :scalar)
            ))

    For complex problems, _zero_constraint_components has already split
    non-Hermitian zero matrices into Hermitian and anti-Hermitian-as-Hermitian
    parts. Both parts produce real-valued linear forms by construction.

7.  Discover primary pivots. For each (b, block) in enumerate(psd_blocks_lin):
        For (i, j) in iteration order (rows first, then columns; upper triangle
        only for HPSD blocks):
            f := block.entries[i, j]
            length(f.terms) == 1 || continue
            (key, coef) := f.terms[1]
            phase := unit_phase(coef)        # one of ±1, ±i; else nothing
            phase === nothing && continue

            # adjoint disambiguation for complex Hermitian blocks
            adj := false
            if cone == :HPSD && haskey(pivots, key) && !haskey(pivots, adjoint_key[key])
                # primary already taken; this entry pivots the adjoint
                key := adjoint_key[key]
                phase := conj(phase)
                adj := true
            end
            haskey(pivots, key) && continue

            pivots[key] := Pivot(b, i, j, phase, adj)

8.  Determine orphans: free_keys := [α for α in moments if !haskey(pivots, α)]

    For each orphan α, decide:
      - real algebras and small problems: leave in free_keys; lowering
        builds free scalars (this is what :free_variables does today).
      - complex algebras and BPSDP target: still leave in free_keys; the
        :aux_psd_free policy is preserved for backward compatibility but is
        not the recommended path (oracle finding: 1×1 aux block ≠ free
        complex moment, and BPSDP behaves better with explicit free scalars).

    The pivot-coverage audit (see § Pivot-coverage audit) determines whether
    free_keys is empty in practice for our target workloads. If it is, the
    BPSDP path needs no orphan handling at lowering time.

9.  Build pivot_at:
    pivot_at := Dict{Tuple{Int,Int,Int}, Vector{K}}()
    for (α, p) in pivots
        push!(get!(pivot_at, (p.block, p.row, p.col), K[]), α)

10. Validate invariants I1..I14. Construct MomentLinearData and attach to
    MomentProblem.
```

`linearize` helper (the one piece of algebra-aware code that runs once per entry):

```julia
function linearize(poly::Polynomial{A,T,C}, moment_index, key_to_monomial) where {A,T,C}
    pairs = Pair{Vector{T}, C}[]
    for (coef, mono) in poly.terms
        for (sub_coef, word) in _simplified_to_pairs(simplify(A, mono.word))
            key = symmetric_canon(expval_word(word))
            push!(pairs, key => coef * sub_coef)
        end
    end
    return LinearMomentForm{Vector{T}, C}(pairs)
end
```

---

## Algorithm: simplified `lowering.jl`

After enrichment, `lowering.jl` shrinks dramatically. Public API stays the same; internals consume `mp.linear`.

```julia
function build_jump_model(mp::MomentProblem;
        formulation::Symbol = :moment_variables,
        representation::Symbol = :real,
        orphan_policy::Symbol = :error)

    # Validate orphan policy against mp.linear.free_keys
    if !isempty(mp.linear.free_keys) && orphan_policy == :error
        throw(ArgumentError("$(length(mp.linear.free_keys)) orphan moment(s); see mp.linear.free_keys"))
    end

    if formulation == :psd_blocks
        _is_complex_problem_K(mp) || throw(ArgumentError(":psd_blocks only for complex algebras"))
        representation == :complex || throw(ArgumentError(":psd_blocks supports :complex only"))
        return _build_psd_blocks_model(mp, orphan_policy)
    end

    return _is_complex_problem_K(mp) ?
        _build_complex_moment_variable_model(mp) :
        _build_real_moment_variable_model(mp)
end

function _build_psd_blocks_model(mp::MomentProblem, orphan_policy::Symbol)
    L = mp.linear
    C = real(eltype(coefficients(mp.objective)))
    model = GenericModel{C}()

    # 1. Declare PSD/HPSD block variables
    X = [_declare_block!(model, b.meta.cone, b.size) for b in L.psd_blocks_lin]

    # 2. Declare free variables for orphans (if any and policy allows)
    free = _declare_free!(model, L.free_keys, orphan_policy, C)

    # 3. Resolver: returns physical ⟨α⟩ as JuMP scalar
    resolver = _make_pivot_resolver(L, X, free)

    # 4. Identity normalization
    @constraint(model, resolver(L.identity) == one(C))

    # 5. Bind PSD/HPSD entries (skip pivot-defining positions)
    for (b, block) in enumerate(L.psd_blocks_lin)
        upper = block.meta.cone == :HPSD
        for i in 1:block.size, j in (upper ? i : 1):block.size
            haskey(L.pivot_at, (b, i, j)) && continue   # binding is tautological
            @constraint(model, X[b][i, j] == _eval_form(block.entries[i, j], resolver, C))
        end
    end

    # 6. Zero constraints
    for zc in L.zero_constraints
        @constraint(model, _eval_form(zc.form, resolver, C) == 0)
    end

    # 7. Objective
    @objective(model, Min, real(_eval_form(L.objective_lin, resolver, C)))

    # 8. Extraction closure
    extract_monomap = () -> Dict(α => value(resolver(α)) for α in L.moments)

    return model, extract_monomap
end

function _make_pivot_resolver(L::MomentLinearData, X, free)
    function resolver(α)
        p = get(L.pivots, α, nothing)
        if p !== nothing
            # ⟨α⟩ = conj(phase) * X[block][row,col]
            return conj(p.phase) * X[p.block][p.row, p.col]
        end
        return free[α]
    end
    return resolver
end

function _eval_form(form::LinearMomentForm{K,C}, resolver, ::Type{T}) where {K,C,T}
    expr = zero(GenericAffExpr{T, VariableRef})
    for (α, c) in form.terms
        add_to_expression!(expr, c, resolver(α))
    end
    return expr
end
```

The lowering body is roughly **80–120 LOC**, down from 587. The deleted material:
- `_all_moment_keys`, `_collect_polynomial_keys!`, `_collect_moment_key!`
- `_pivot_candidate`, `_strict_unit_phase`
- `_discover_pivots_unchecked`, `discover_pivots`, `orphan_keys`
- `_pivot_entry_lookup`
- `substitute(::Polynomial, ::Resolver)`

These are replaced by direct field reads on `mp.linear`.

---

## Algorithm: simplified `sos.jl`

`_sos_dualize_real` and `_sos_dualize_hermitian` currently walk polynomial matrices. With the cache, they walk linear forms instead. The factor-of-2 trap stays — it's a property of the Hermitian inner product, not of the representation.

Pseudocode:

```julia
function _sos_dualize_complex(mp::MomentProblem)
    L = mp.linear
    RC = real(eltype(coefficients(mp.objective)))
    dual_model = GenericModel{RC}()

    # Per-block dual variables: lifted real 2n×2n PSD for HPSD, n×n symmetric free for Zero
    G = [_declare_dual!(dual_model, block) for block in L.psd_blocks_lin]
    Z = [_declare_zero_dual!(dual_model, _zero_block_size(zc)) for zc in L.zero_constraints]

    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    # One coefficient-matching equation per (canonical key, real/imag component)
    eqs_re = [zero(...) for _ in L.moments]
    eqs_im = [zero(...) for _ in L.moments]

    # Objective contribution
    for (α, c) in L.objective_lin.terms
        idx = L.moment_index[α]
        eqs_re[idx] += real(c)
        eqs_im[idx] += imag(c)
    end
    eqs_re[L.moment_index[L.identity]] -= b

    # PSD-block contribution: for each block entry, accumulate Tr(G_b * E_{ij}) per key
    for (b, block) in enumerate(L.psd_blocks_lin)
        for i in 1:block.size, j in 1:block.size
            for (α, c) in block.entries[i,j].terms
                idx = L.moment_index[α]
                # see _sos_dualize_hermitian for the factor-of-2 derivation
                _accumulate_dual_contribution!(eqs_re, eqs_im, idx, c, G[b], i, j, block.meta.cone)
            end
        end
    end

    # Zero-constraint contribution: Z_k * f
    for (k, zc) in enumerate(L.zero_constraints)
        for (α, c) in zc.form.terms
            idx = L.moment_index[α]
            eqs_re[idx] += real(c) * Z[k]
            eqs_im[idx] += imag(c) * Z[k]
        end
    end

    @constraint(dual_model, eqs_re .== 0)
    @constraint(dual_model, eqs_im .== 0)

    return SOSProblem(dual_model, length(L.moments))
end
```

Key point: **the factor-of-2 helper `_accumulate_dual_contribution!` is the only place the Hermitian inner-product convention lives.** The existing factor-of-2 regression test in `test/relaxations/interface.jl` must still pass.

---

## Pivot-coverage audit (precondition)

Before declaring this design final, run this audit and record results in `output/phase2/pivot_coverage_audit.md`:

```julia
# probes/pivot_coverage_audit.jl
for case in [
    h2_chain_nk2,                          # primary BPSDP target
    h4_chain_nk2_proxy_small,              # downscaled H₄
    pauli_3qubit_ground_state,             # small Pauli
    fermionic_4site_hubbard,               # small Fermionic
    bosonic_2mode_truncated,               # small Bosonic
    cs_failure_E10,                        # known correlated-sparsity edge case
]
    mp = build_moment_problem(case)
    pivots = discover_pivots(mp)           # use existing runtime discovery
    orphans = orphan_keys(mp, pivots)
    record(case, length(orphans), categorize_orphans(orphans, mp))
end
```

The audit answers two questions:

1. **Does I7 (every moment has a pivot or is in `free_keys`) hold trivially?** If orphans only ever appear for known reasons (e.g., one-sided moment-equality rows that don't enter any PSD block), the orphan handling stays as `free_keys` + free scalar variables in lowering.

2. **Are there algebra/sparsity combinations where orphans dominate?** If yes, the `:aux_psd_free` policy stays but is documented as "compatibility, not recommended."

Result categories the audit must distinguish:
- **No orphans** → cleanest case; `free_keys` empty in production.
- **Few orphans (< 0.1% of moments)** → `:free_variables` is the right default.
- **Many orphans** → algebra-specific issue; need investigation, possibly a moment-relax fix.

This audit is the gate before the implementation can ship — without it we'd be designing for hypothetical cases.

---

## Migration order

Each step is an independently-mergeable PR with its own tests.

### PR 1 — Add types, no behavior change

- Add `LinearMomentForm`, `Pivot`, `BlockMeta`, `BlockOrigin` hierarchy, `PSDBlockLin`, `ScalarLinearConstraint`, `ConstraintOrigin` hierarchy, `MomentLinearData` to `src/optimization/moment_linear.jl` (new file).
- Don't touch `MomentProblem` yet.
- Unit tests for `LinearMomentForm` smart constructor (sorting, dedup, zero-pruning) and key ordering (`key_lt`).

### PR 2 — Pivot-coverage audit

- Add `probes/pivot_coverage_audit.jl`.
- Run on HAI; commit results to `output/phase2/pivot_coverage_audit.md`.
- This determines whether subsequent PRs can drop `:aux_psd_free` from default suggestions.

### PR 3 — Enrich `MomentProblem`

- Add `linear::MomentLinearData{K,C,M}` field to `MomentProblem`.
- Update `moment_relax` to populate the cache as the last step.
- Provide a `MomentProblem{A,T,M,P}(...)` migration shim so old call sites still work.
- All existing tests pass with no edits — the cache is populated but unused.

### PR 4 — Rewrite `lowering.jl`

- Replace `discover_pivots`-driven path with cache reads.
- Fix the four bugs:
  - Use `conj(p.phase) * X[...]` not `X[...] / p.phase`.
  - Use `add_to_expression!` not `expr += coef * resolver(α)`.
  - Hermitian upper-triangle pivots use `pivot.adjoint`.
  - Identity normalization uses `L.identity`, not `findfirst`.
- Drop deleted helpers (see § Algorithm: simplified `lowering.jl`).
- Existing lowering tests must all pass with no edits. **This is the gate** — semantic equivalence with the current lowering on every COSMO-stable test.

### PR 5 — H₂/Nk=2 bridge smoke test

- Add `test/relaxations/h2_nk2_bridge_smoke.jl` (predicate from the old plan; preserved verbatim below).
- Run on HAI per server policy.
- Acceptance: < 100 scalar 1×1 PSD cones, vs current 406,706.

### PR 6 — Migrate `sos_dualize` to cache

- Rewrite `_sos_dualize_real` and `_sos_dualize_hermitian` to consume `mp.linear`.
- The factor-of-2 regression test in `test/relaxations/interface.jl` is the gate.
- Document `_accumulate_dual_contribution!` as the single home for Hermitian inner-product scaling.

### PR 7 — Remove dead code

- Delete `discover_pivots`, `_all_moment_keys`, `_pivot_candidate`, `_strict_unit_phase`, `_pivot_entry_lookup`.
- Update docstrings in `lowering.jl` and `moment.jl`.

---

## Acceptance test (the H₂/Nk=2 smoke gate)

```julia
# test/relaxations/h2_nk2_bridge_smoke.jl
using NCTSSoS, BPSDP, MathOptInterface
const MOI = MathOptInterface

@testset "H2/Nk=2 :psd_blocks delivers BPSDP-favorable cones" begin
    mp = build_h2_nk2_moment_problem()        # from demos/h2_periodic_nk2_moment_sos.jl
    model, _extract = build_jump_model(mp;
        formulation    = :psd_blocks,
        representation = :complex,
    )
    set_optimizer(model, BPSDP.Optimizer)
    MOI.Utilities.attach_optimizer(model)

    backend = MOI.get(unsafe_backend(model), MOI.ListOfConstraintTypesPresent())
    counts  = Dict(t => MOI.get(unsafe_backend(model), MOI.NumberOfConstraints{t...}())
                   for t in backend)

    n_hermitian = sum(v for ((F,S),v) in counts
                      if S <: MOI.HermitianPositiveSemidefiniteConeTriangle; init=0)
    n_psd_1x1   = sum(v for ((F,S),v) in counts
                      if S <: MOI.PositiveSemidefiniteConeTriangle &&
                         dimension_is_1x1(F, S); init=0)

    n_symbolic_blocks_hermitian = count(b -> b.meta.cone == :HPSD, mp.linear.psd_blocks_lin)

    @test n_hermitian == n_symbolic_blocks_hermitian
    @test n_psd_1x1 <= length(mp.linear.free_keys)
    @test n_psd_1x1 <  100                    # absolute ceiling vs current 406,706
end
```

Three assertions:
1. **Cone count match.** Number of `HermitianPositiveSemidefiniteConeTriangle` instances equals the number of `:HPSD` blocks. No silent realification by MOI bridges.
2. **No 1×1 clutter beyond orphans.** Under `:error` policy the second is zero.
3. **Hard ceiling.** < 100 vs the current 406,706.

If (1) fails, MOI is silently realifying `HermitianPSDCone()` — investigate which bridge.
If (3) fails, a bridge is misbehaving and the refactor is incomplete.

---

## Open contracts (decide before coding)

| ID | Contract | Default |
|---|---|---|
| **D1** | `K = Vector{T}` (matches today's `symmetric_canon` output) | Use `Vector{T}`; document the mutability risk; revisit only if profiling shows Dict overhead |
| **D2** | `key_lt` is `Base.lexless` on the underlying `Vector{T}` | Lexicographic; matches `_sorted_symmetric_basis` |
| **D3** | Pivot iteration order: outer over blocks, inner row-major upper-triangle for HPSD, full matrix for PSD | Deterministic by data layout |
| **D4** | When two keys α and α† pivot at the same `(b,i,j)` of an HPSD block, the lex-smaller key gets `adjoint=false` | Lex order on canonical keys |
| **D5** | Phase coefficient set: `{±1, ±i}` exactly | Same as today |
| **D6** | Identity normalization is emitted unconditionally as a regular constraint | Same as old C6 |
| **D7** | Hermitian zero matrices are upper-triangle bound only; the lower triangle is implied by the cone | Same as today's `_add_zero_bindings!` |
| **D8** | `mp.constraints`, `mp.objective`, `mp.total_basis` are immutable after `moment_relax` returns | Constructor invariant; assertion at `linear`-cache populate time |
| **D9** | `expval_word(word)` is the canonical-key projection used by `linearize` | Match the existing `symmetric_canon(expval(...))` chain |

---

## Known traps

Carried forward from the analysis material:

1. **`Symmetric(M) in PSDCone()` ≠ `M in PSDCone()`**. JuMP wraps the latter via `SquareBridge` into triangular cone + off-diagonal equalities. For PSD blocks declared with `@variable(model, [1:n,1:n] in PSDCone())` this is moot, but if a future builder constructs a `Symmetric` wrapper from affine expressions, beware.

2. **Real-lift PSD-first needs four bindings per Hermitian entry.** `:psd_blocks, representation = :real` is out of scope for this plan; if added later, every entry of a real-lift block must bind all four positions: `R[i,j], R[i,n+j], R[n+i,j], R[n+i,n+j]`.

3. **Pivot phases include `±1, ±i`; phase must be tracked.** Encoded in `Pivot.phase`. Resolver multiplies by `conj(phase)`.

4. **Native `HermitianPSDCone()` survival is bridge-dependent.** The smoke test asserts it; if it fails with a particular MOI version, fix MOI or supply `representation = :real` for that backend.

5. **`symmetric_canon` does not collapse `w ↔ w†` for complex algebras.** This means `α` and `adjoint_key[α]` are distinct keys in `moments`. The cache stores both and `pivot_at` maps both to the same `(b,i,j)` for HPSD blocks under upper-triangle storage.

6. **`Vector{T}` as `Dict` key is mutable.** Don't mutate keys post-construction. The cache treats canonical keys as values; constructing them via `symmetric_canon` returns fresh vectors each time, but downstream code must not call `push!`/`pop!` on retrieved keys.

7. **`zero_constraints` are real linear forms, not complex.** The Hermitian/anti-Hermitian split in `_zero_constraint_components` produces matrices whose entries become real linear forms after canonicalization. The constructor asserts this by checking `imag(c) ≈ 0` for every coefficient.

---

## File-level change list

| File | Change |
|---|---|
| `src/optimization/moment_linear.jl` | **New**. Type definitions: `LinearMomentForm`, `Pivot`, `BlockOrigin` hierarchy, `BlockMeta`, `PSDBlockLin`, `ConstraintOrigin` hierarchy, `ScalarLinearConstraint`, `MomentLinearData`. Smart constructors and invariant assertions. |
| `src/optimization/moment.jl` | Add `linear::MomentLinearData{K,C,M}` field to `MomentProblem`. Update `moment_relax` to populate the cache after `_add_parity_constraints!` and `_add_moment_eq_constraints!`. Add migration shim constructor. Add `linearize` helper. |
| `src/optimization/lowering.jl` | Strip ~400 LOC of runtime discovery. Rewrite `_build_*` builders to consume `mp.linear`. Apply the four bug fixes. Public `build_jump_model` signature unchanged. |
| `src/optimization/sos.jl` | Rewrite `_sos_dualize_real` and `_sos_dualize_hermitian` to consume `mp.linear`. Factor Hermitian inner-product scaling into `_accumulate_dual_contribution!`. |
| `src/optimization/interface.jl` | Update result-struct types if they expose `MomentProblem` parameters. |
| `probes/pivot_coverage_audit.jl` | **New**. Run audit, write `output/phase2/pivot_coverage_audit.md`. |
| `test/relaxations/lowering.jl` | New unit tests for `LinearMomentForm` constructor, pivot computation, adjoint disambiguation. |
| `test/relaxations/h2_nk2_bridge_smoke.jl` | **New**. HAI gate. |
| `output/phase2/pivot_coverage_audit.md` | **New**. Audit results. |
| `AGENTS.md` | Update references from `LOWERING_REFACTOR_PLAN.md` and `MOMENT_SOS_PIPELINE_ANALYSIS.md` to this file. |

---

## Estimated change footprint

- **New code**: `moment_linear.jl` (~250 LOC types + invariants), `linearize` helper in `moment.jl` (~50 LOC), unit tests (~200 LOC), pivot-coverage probe (~100 LOC), smoke test (~80 LOC).
- **Removed code**: ~400 LOC of runtime discovery in `lowering.jl`.
- **Modified code**: `moment.jl` (~60 LOC for cache populate + assertions), `lowering.jl` (~150 LOC rewritten as cache reads), `sos.jl` (~200 LOC rewritten as cache reads).

**Net**: roughly −100 LOC, but the *real* win is that the seven open contracts collapse to constructor invariants and four bugs get fixed.

---

## What this enables downstream

| Phase | Without this | With this |
|---|---|---|
| Phase 2 diagnostic on H₂/Nk=2 | Reads anonymous polynomial matrices; no per-block provenance | `mp.linear.psd_blocks_lin[b].meta` has block origin and row labels; symmetry sector inference attaches a separate `Dict{Int, SymmetryTag}` keyed by block index |
| Phase 3 production at H₄/Nk=2 | Cannot build clean JuMP shape; bridge clutter swamps BPSDP | Build is clean. Remaining work: symmetry block-diagonalization, Phase-2-recommended solve strategy, BPSDP tuning |
| COSMO compatibility | n/a | Unchanged. `:moment_variables, :real` is still default |
| SOS dual factor-of-2 | Test-protected but logic spread across `_sos_dualize_hermitian` | Confined to `_accumulate_dual_contribution!`; same test still gates |

---

## Bottom line

This plan moves linear-form discovery from runtime (in `lowering.jl`) into construction-time (in `moment_relax`), caches it on `MomentProblem` via `MomentLinearData`, and fixes four bugs (Hermitian pivot injectivity, divide-vs-conj, non-deterministic accumulation, untyped origin). It is the smallest cut that makes:

- `lowering.jl` readable in a single sitting,
- `sos.jl` factor-of-2 trap a one-helper concern,
- Phase 2 provenance a first-class feature, and
- H₂/Nk=2 (and H₄/Nk=2) buildable in a JuMP shape BPSDP can consume.

Implementation is gated on the pivot-coverage audit (PR 2) and the H₂/Nk=2 bridge smoke test (PR 5). Ship in the order PR 1 → PR 2 → PR 3 → PR 4 → PR 5 → PR 6 → PR 7.
