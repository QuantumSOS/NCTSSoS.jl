# Paper $^2D$ block structure vs. NCTSSoS term-sparsity blocks — H4 periodic

Source paper: Sajjan, Schlimgen, Warren, Mazziotti, *arXiv:2504.02861* (Apr 2025).
Target code: `NCTSSoS.jl` as of this repo snapshot.
Demo data: `demos/h4_periodic_term_sparsity_size.jl` (Nk=2, [4,8] active space).

---

## 1. Two different objects to block-diagonalize

| | Paper (periodic V2RDM) | NCTSSoS TS pipeline |
|---|---|---|
| Matrix being blocked | The **2-RDM** $^2D$ and its partners $^2Q, ^2G$ (each of size $r^2 N_k \times r^2 N_k$). | The **moment matrix** $M_L[b_i^\dagger b_j]$ over the full order-$d$ monomial basis, plus every localizing matrix. |
| What goes in a row/column | A pair operator: $a_{p k_1} a_{q k_2}$ (for $^2D$), $a^\dagger_{p k_1} a_{q k_2}$ (for $^2G$), $a^\dagger_{p k_1} a^\dagger_{q k_2}$ (for $^2Q$). Exactly degree-2 fermionic words. | Every normal-ordered word of total degree $\le d$, regardless of particle number, spin, or momentum. For H4 Nk=2, order 2: 2,081 basis words. |
| Block count | $N_k$ sectors, indexed by total pair-momentum $K = k_1 + k_2 \pmod{G}$ (Eq. 16, p. 3). | One per maximal clique of the chordal extension of the term-sparsity graph. |
| Block size | $\approx r^2 N_k$ per sector. For Nk=2 [4,8] that is 240 (K=0) or 256 (K=1). | Data-dependent; for Nk=2 [4,8] at order 2 with MMD: **313 blocks, max 168×168** (see demo). |
| Non-zero element count (paper metric) | Drops from $r^4 N_k^4$ to $r^4 N_k^3$ — factor $N_k$ savings. | Drops from $\binom{2n}{\le 4} \approx 6.8\times 10^5$ unique moments to $\approx 2.7\times 10^4$ (demo number). |

## 2. How each construction decides who connects to whom

**Paper.** Entries $^2D^{p_1 k_1, p_2 k_2}_{p_3 k_3, p_4 k_4}$ are zero unless $k_1+k_2 \equiv k_3+k_4 \pmod G$. This is a *group-theoretic selection rule* from the abelian translation symmetry — every irreducible representation of $\mathbb{Z}_{N_k}$ is 1-dimensional and labeled by $K$. The 2-RDM decomposes as

$$^2D \;=\; \bigoplus_{K=0}^{N_k-1} \, ^2D^{(K)}.$$

The blocks live in the rep label, not in the chosen orbital basis. The paper extends the same selection rule to $^2Q, ^2G$ (all three carry a pair momentum) and the 1-RDM (momentum-diagonal).

**NCTSSoS.** Nodes are the basis monomials $b_i$. An edge $b_i \sim b_j$ exists iff, for some monomial $g$ in the constraint polynomial (or the identity, for the moment matrix), the normal-ordered product $b_i^\dagger g\, b_j$ lies in the current *activated support* — a set seeded from the objective/constraints and then iterated (`iterate_term_sparse_supp` in `src/optimization/sparsity.jl`). The edge test is **literal monomial equality** after normal ordering, not a symmetry projection. Chordal extension + maximal-clique decomposition then hands out the PSD blocks.

The crucial consequence, observed empirically in `h4_periodic_term_sparsity_size.jl`:

> All PSD blocks are K-pure: TS does recover momentum symmetry, but it *shards each K sector into many smaller PSD blocks*.

Why? Because inside (say) $K=0$ there are many sub-words that only appear together in a specific bilinear; anything that never co-occurs in an activated-support entry stays disconnected. The chordal extension only bridges them if there is already a chordal path — and for H4 at order 2 there isn't enough connectivity to reassemble the single dense sector block that the group-theoretic projection gives.

## 3. What the two approaches share, and what they do not

Shared:
- Both exploit the same underlying symmetry (translation / momentum conservation) in the sense that **no block mixes $K=0$ with $K=1$**.
- Both produce block-diagonal positive-semidefinite constraints that replace one dense one.

Different in kind:
- Paper blocks are **irreducible under an abelian group action** — they are as coarse as the symmetry allows. The paper's $^2D^{(K)}$ block is *one* matrix.
- TS blocks are **chordal cliques on a graph defined by concrete monomial supports** — they are as coarse as the specific monomial algebra allows. TS gives dozens of small $K=0$ blocks whose direct sum is contained inside the paper's single $K=0$ block.
- Paper only certifies 2-positivity on pair operators ($^2D$, $^2Q$, $^2G$); NCTSSoS at order 2 certifies PSD on the full degree-$\le 2$ moment matrix. The NCTSSoS relaxation is **strictly stronger** (tighter lower bound) but also **strictly bigger**.

A one-line summary: *TS refines the symmetry blocks when it can, but it never coarsens a sector into a single group-theoretic block*.

## 4. Current scalability wall for H4 Nk=2 [4,8]

From `test/data/expectations/h4_periodic_v2rdm.toml` and the benchmark header:

| Route | scalar moments | PSD blocks (max) | constraints | Solve status |
|---|---:|---:|---:|---|
| Order-2 CS-only (no TS) | 679,121 | 1 × **2081** | 508,611 | Mosek OOM (`MSK_RES_ERR_SPACE`) |
| Order-2 CS + TS (MMD), SOS dual route | — | 313 (max 168) | 172,195 vars, 2,006,058 JuMP vars | COSMO did not terminate in 45 min |
| Order-2 CS + TS (MMD), primal moment route | 26,817 | 313 (max 168) | (assembly OK) | Stage-level workable but still dominated by blocks near 168 |

The paper's equivalent (V2RDM[4,8] at Nk=2) sees one PSD sector per K: roughly **240×240 and 256×256**. Two matrices total. Converges in minutes with a boundary-point method.

So the cost gap has two factors compounding:

1. **NCTSSoS's degree-$\le 2$ moment matrix is ~10× wider than the paper's pair-only matrix.** (2,081 basis words vs. $\sim 500$ pair operators, for the same $n=32$ modes.)
2. **Inside each K sector the moment matrix is further fractured by TS** into many small blocks whose sum of squared sizes is larger than the single coarse sector block would be.

## 5. Concrete ways to scale NCTSSoS on this class of problems

Ranked from biggest lever to smallest, and from least to most invasive to the codebase.

### (A) Use the explicit-`moment_basis` API to reproduce the paper's 2-RDM shape

`correlative_sparsity(pop, moment_basis::AbstractVector, elim_algo)` already exists (see `src/optimization/sparsity.jl`). Feed it a *pair-only* basis:

```julia
# pair operators that populate rows/columns of ^2D
pair_basis = NormalMonomial[
    a_up[p] * a_up[q]                        # for ^2D, σ = τ = ↑
    for p in 1:n for q in (p+1):n            # antisymmetric pair
]
append!(pair_basis, a_up[p] * a_dn[q]  for p in 1:n, q in 1:n)  # σ=↑, τ=↓
# ... and similarly for ^2Q (aa†…), ^2G (a†a)
pushfirst!(pair_basis, one(NormalMonomial))  # identity is mandatory
```

Expected impact for Nk=2 [4,8]:

- basis drops from 2,081 to O(500) → moment matrix ≤ 500×500 pre-TS.
- unique moments drop from 679k toward 10k-ish (rough scaling, dominated by degree-$\le 4$ words reachable from pair-pair products).
- This is the paper's scale, modulo the extra 1-RDM rows.

Gotchas:
- You lose order-2 certification strength; the resulting SDP is now a *pair-2-positivity* relaxation, not a full order-2 Lasserre. That is exactly the trade the paper makes.
- You must include the identity and the degree-1 words if you want the 1-RDM constraint.

### (B) Add symmetry-adapted basis partitioning *before* TS

The demo already demonstrates that every TS block is K-pure. Formalize that into the pipeline so that:

1. A user supplies a `quantum_numbers(m::NormalMonomial) -> Tuple`. For H4 that is `(ΔN, K, S_z)`.
2. `compute_sparsity` partitions `clq_mom_mtx_bases` by quantum-number tuple *before* calling `init_activated_supp` and `term_sparsities`.
3. Each partition is then one pre-declared PSD block. TS runs *inside* each partition to sparsify further if desired.

This is the minimal change that turns the existing TS-is-K-pure observation into a guaranteed invariant and also produces the coarse paper-style blocks when TS cannot fragment further.

Rough implementation sketch:

```julia
# in src/optimization/sparsity.jl, CorrelativeSparsity path
function _partition_by_quantum_numbers(basis::Vector{M}, qn::Function) where {M}
    d = Dict{Any, Vector{M}}()
    for b in basis
        push!(get!(d, qn(b), M[]), b)
    end
    return collect(values(d))
end
```

and then wrap `term_sparsities` so that it receives a *vector of partitions* and returns a `TermSparsity` with `block_bases` that is the concatenation of per-partition blocks. The moment matrix assembly (`_build_constraint_matrix` in `src/optimization/moment.jl`) already consumes a per-block basis, so nothing downstream needs to change.

Expected impact for H4 Nk=2 [4,8], order 2:

- At most $2 \times (2N_{\max ΔN}+1) \times (2S_{\max}+1)$ sectors, i.e. on the order of tens.
- Each sector ≤ few hundred basis words, matching the paper's block size.
- TS inside the sector either splits it further or leaves it alone — you keep the paper's guarantee and optionally get more.

### (C) Expose a `symmetry_sectors` kwarg on `polyopt` / `cs_nctssos`

Give users a first-class way to pass in a symmetry action (abelian character table, or a tagging function). Internally:

- For each basis word, compute its irrep label.
- Assert that every constraint polynomial is a sum of words with *consistent* irrep labels (this is the fermionic parity check generalized — and you already do that for $\mathbb Z_2$ fermion parity in `_add_parity_constraints!` in `src/optimization/moment.jl`).
- Block the moment matrix accordingly.

This generalizes the existing fermionic parity machinery from $\mathbb Z_2$ to any finite abelian group. That covers translation, spin-Sz, and particle number in one sweep.

### (D) Avoid the SOS-dual blow-up on problems this large

The benchmark header explicitly notes:

> A first run with `dualize=true` (cs_nctssos's default) built a JuMP SOS dual with 2,006,058 variables and 172,195 constraints, and COSMO did not terminate within 45 minutes. The primal moment route keeps the same 313 PSD blocks, but the assembled Hermitian JuMP model is 171,882 scalar vars representing 26,817 unique symbolic moments.

So for wide problems with many zero-cone constraints, primal-moment assembly is already about **11.7×** smaller at the JuMP-model level. The larger **75×** gap is the structural comparison between dual JuMP vars and the symbolic unique-moment count. Two things to do:

1. Document it and make it selectable in `SolverConfig` (current exposure is only via the demo script assembling a primal model by hand).
2. Investigate whether `sos_dualize` should dualize all zero-cone blocks at all. On this asset, most of the dual blow-up comes from lifted matrix multipliers on the HPSD blocks, but the many scalar zero-cone blocks still add extra lifted multipliers; dualizing only the PSD cones and keeping zeros primal may be the cleaner route.

### (E) Boundary-point solver hook instead of first-order ADMM or IPMs

The paper solves its SDPs with a boundary-point method (D. A. Mazziotti, PRL 106, 083001 (2011); Ref. [47] in the paper). COSMO's first-order ADMM and Mosek's interior-point both struggle on this particular problem shape (many small PSD blocks + a huge equality block). Adding a MOI-compatible boundary-point backend (e.g. via a thin wrapper around an existing C/Fortran implementation, or a pure-Julia prototype) would match the paper's toolchain.

Realistic order of implementation: (A) immediately, (B) next, (C) after the API shape is locked, (D) concurrently with (A) as a solver-config switch, (E) research-grade.

## 6. Short answer to the original question

- The paper's $^2D$ blocks come from *symmetry projection* onto irreps of the translation group. They are the coarsest possible blocks compatible with momentum conservation. One block per total pair momentum $K$.
- NCTSSoS term sparsity instead derives blocks from the *chordal structure of the specific monomial-product graph*. It **respects** the same symmetry (no block crosses $K$) but *further fragments* each symmetry sector into many small blocks, and it operates on a much wider basis (all order-$\le 2$ words, not just pair operators).
- On H4 Nk=2 [4,8] at order 2, TS reduces the dense 2081×2081 block to 313 blocks, max size 168, but that is still well short of the paper's two sector blocks at 240 / 256.
- The practical path to scale NCTSSoS on this problem family is **(A)** restrict the moment basis to paper-style pair operators (big, immediate), combined with **(B)** symmetry-adapted basis partitioning before TS (moderate effort, turns the empirical K-purity into a guaranteed coarser blocking), and switch default assembly to the primal-moment route **(D)**.
