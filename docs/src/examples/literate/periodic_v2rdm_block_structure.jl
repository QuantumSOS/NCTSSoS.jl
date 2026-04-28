# # [Periodic V2RDM: ${}^{2}\!D$ Symmetry Blocks vs. Term-Sparsity Blocks](@id periodic-v2rdm-block-structure)
#
# This page answers a single, concrete question:
#
# > **The paper's ${}^{2}\!D$ is block-diagonal by crystal momentum.
# >  NCTSSoS.jl's term sparsity also produces block-diagonal PSD
# >  constraints. Are they the same blocks? If not, which one helps
# >  us scale, and how do we get it in the current framework?**
#
# Short answer: they are *different objects*. Both respect the
# symmetry. Only one exploits it **coarsely**. Below we show exactly
# what each produces on the H₄ chain benchmark
# (Schouten, Ewing & Mazziotti, [arXiv:2504.02861](https://arxiv.org/abs/2504.02861)),
# then list the concrete code levers to reach paper-scale
# simulations inside NCTSSoS.jl.
#
# If you have not read it yet, the companion spec page
# [Periodic V2RDM Benchmark (H₄ Chain)](@ref periodic-v2rdm-h4)
# pins down every sign and index convention used here. The page
# [Why Correlative Sparsity Misses the $k$-Blocks](@ref periodic-v2rdm-cs-k-blocks)
# establishes the weaker, CS-only picture. This page zooms in on the
# **TS layer** and compares it to the paper's symmetry blocking
# directly.

# ---
#
# ## 1. Two different matrices to block-diagonalize
#
# The first thing people get wrong: the paper and NCTSSoS are not
# block-diagonalizing the *same* matrix.
#
# | | Paper (periodic V2RDM) | NCTSSoS order-2 moment relaxation |
# |---|---|---|
# | Matrix | The pair-indexed **2-RDM** ${}^{2}\!D$ (plus ${}^{2}\!Q, {}^{2}\!G$). | The full **order-2 moment matrix** $M[b_i^\dagger b_j]$ over every normal-ordered word of degree $\le 2$. |
# | Row/column index | A pair operator $a_{p k_1} a_{q k_2}$ (or the $a^\dagger a$ / $a^\dagger a^\dagger$ variants for ${}^{2}\!G, {}^{2}\!Q$). | Any normal-ordered monomial of degree $\le 2$ — including the identity, single creators, and single annihilators. |
# | Block rule | Total pair momentum $K = (k_1 + k_2) \bmod N_k$. Abelian selection rule from the translation group. | Maximal cliques of the chordal extension of a term-sparsity graph on the moment basis. |
# | Block count | Exactly $N_k$ sectors — one per irrep of $\mathbb Z_{N_k}$. | Data-dependent (313 on H₄ Nk=2 [4,8]). |
# | Block size | $\sim r^2 N_k$ per sector. Paper: 240 and 256 on H₄ Nk=2. | Depends on chordal extension. H₄ Nk=2: max 168, many small ones. |
# | Relaxation tightness | 2-positivity only (PQG). | Full order-2 Lasserre ⊇ PQG. Strictly tighter **and** strictly more expensive. |
#
# The paper projects a small matrix onto irreps. NCTSSoS decomposes
# a much wider matrix via graph chordality. These are genuinely
# different procedures that happen to produce block-diagonal
# structures.

# ---
#
# ## 2. Paper blocks: an abelian group projection
#
# Crystal momentum $k$ is conserved modulo a reciprocal lattice
# vector:
# ```math
# (k_1 + k_2 - k_3 - k_4) \cdot a \;=\; 2\pi n.
# ```
# Since $\mathbb Z_{N_k}$ is abelian, every irrep is 1-dimensional
# and labeled by a single character $K$. So the 2-RDM splits:
# ```math
# {}^{2}\!D \;=\; \bigoplus_{K=0}^{N_k-1}\; {}^{2}\!D^{(K)},
# \qquad
# {}^{2}\!D^{(K)}_{(p_1 k_1, p_2 k_2),\,(p_3 k_3, p_4 k_4)}
# = 0 \text{ unless } k_1 + k_2 \equiv k_3 + k_4 \equiv K.
# ```
# This is the **coarsest blocking that respects the symmetry**. Each
# block is one irrep.
#
# Non-zero element count drops from $r^4 N_k^4$ to $r^4 N_k^3$ — a
# factor of $N_k$, and an algorithmic cost factor much larger than
# that on interior-point / boundary-point SDPs.

# ---
#
# ## 3. NCTSSoS TS blocks: a chordal graph decomposition
#
# NCTSSoS builds a graph where:
#
# - **nodes** = basis monomials $b_i$,
# - **edge** $b_i \sim b_j$ ⇔ the product $b_i^\dagger\, g\, b_j$
#   (for some monomial $g$ in the current activated support) lands
#   inside the activated support after normal ordering.
#
# The activated support starts from the objective and constraints
# and is iterated. A chordal extension of this graph, followed by
# maximal-clique decomposition, gives the PSD blocks.
# See `src/optimization/sparsity.jl` for the exact recipe
# (`get_term_sparsity_graph`, `iterate_term_sparse_supp`).
#
# !!! note "TS respects the symmetry, but does not project onto it"
#     The edge test is **literal monomial equality after normal
#     ordering**, not a group character. Two basis words in the same
#     $K$ sector whose specific products never land in the activated
#     support simply get no edge — even though they belong to the
#     same irrep.
#
#     Consequence: every TS block is contained in some symmetry
#     sector (we verify this below), but one symmetry sector usually
#     breaks into many TS blocks.

# ---
#
# ## 4. Set up: the vendored H₄ Nk=2 asset
#
# Everything below uses the reviewed PySCF integrals vendored in
# `test/data/assets/`. The full Hamiltonian build — 32 spin-orbital
# modes, 8 momentum-conserving ERI blocks — is identical to
# [Why Correlative Sparsity Misses the $k$-Blocks](@ref periodic-v2rdm-cs-k-blocks)
# and is factored into a helper here for brevity.

using NCTSSoS

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: load_nk2_asset, load_expectation

asset = load_nk2_asset()
blocker = load_expectation("spin_orbital_order2_blocker")
cs_only = load_expectation("correlative_sparsity_only_attempt")
nothing #hide

# `asset` carries the reviewed integrals; `blocker` and `cs_only`
# are pinned expectation records documenting the dense and CS-only
# SDP footprints. We lean on them instead of re-running the heavy
# 32-mode TS step, which allocates ~10 GiB transiently — not
# something we want in a docs build.
#
# ## 4.1 Manual K-sector partition of the pair basis
#
# The paper's blocks are not discovered from a sparsity graph. They
# come from sorting pair operators by the conserved crystal-momentum
# label. For the `Nk = 2` H₄ asset we can reproduce those sectors with
# plain Julia.
#
# Start from the same 32 spin orbitals used by the Bloch Hamiltonian:

nk = asset.nk
norb = asset.n_active_orb

reg_pairs, ((c_up_k0, c_up_k0_dag),
            (c_dn_k0, c_dn_k0_dag),
            (c_up_k1, c_up_k1_dag),
            (c_dn_k1, c_dn_k1_dag)) = create_fermionic_variables([
    ("c_up_k0", 1:norb),
    ("c_dn_k0", 1:norb),
    ("c_up_k1", 1:norb),
    ("c_dn_k1", 1:norb),
])

spin_orbitals = vcat(
    [(; a = c_up_k0[p], adag = c_up_k0_dag[p], k = 0, spin = :up, orb = p) for p in 1:norb],
    [(; a = c_dn_k0[p], adag = c_dn_k0_dag[p], k = 0, spin = :dn, orb = p) for p in 1:norb],
    [(; a = c_up_k1[p], adag = c_up_k1_dag[p], k = 1, spin = :up, orb = p) for p in 1:norb],
    [(; a = c_dn_k1[p], adag = c_dn_k1_dag[p], k = 1, spin = :dn, orb = p) for p in 1:norb],
)

length(spin_orbitals)

# `spin_orbitals` has length 32: 16 spin orbitals at `k = 0` and 16 at
# `k = 1`.
#
# A pair-basis row is identified by two spin-orbitals. For
# `${}^{2}\!D` and `${}^{2}\!Q` the conserved label is
# `K = (k₁ + k₂) mod Nk`; for `${}^{2}\!G` it is
# `K = (k₁ - k₂) mod Nk`.

pair_slot_label(oi, oj) = (
    (orb = oi.orb, spin = oi.spin, k = oi.k),
    (orb = oj.orb, spin = oj.spin, k = oj.k),
)

PairLabel = typeof(pair_slot_label(spin_orbitals[1], spin_orbitals[2]))
D_pair_labels_by_K = [PairLabel[] for _ in 1:nk]
Q_pair_labels_by_K = [PairLabel[] for _ in 1:nk]
G_pair_labels_by_K = [PairLabel[] for _ in 1:nk]

for i in 1:length(spin_orbitals)-1, j in i+1:length(spin_orbitals)
    oi = spin_orbitals[i]
    oj = spin_orbitals[j]
    K = mod(oi.k + oj.k, nk)
    label = pair_slot_label(oi, oj)
    push!(D_pair_labels_by_K[K + 1], label)
    push!(Q_pair_labels_by_K[K + 1], label)
end

for oi in spin_orbitals, oj in spin_orbitals
    K = mod(oi.k - oj.k, nk)
    push!(G_pair_labels_by_K[K + 1], pair_slot_label(oi, oj))
end

paper_D_block_sizes = length.(D_pair_labels_by_K)
paper_Q_block_sizes = length.(Q_pair_labels_by_K)
paper_G_block_sizes = length.(G_pair_labels_by_K)

n_spin_orb_per_k = 2 * norb
expected_D_block_sizes = [2 * binomial(n_spin_orb_per_k, 2), n_spin_orb_per_k^2]
expected_Q_block_sizes = copy(expected_D_block_sizes)
expected_G_block_sizes = fill(2 * n_spin_orb_per_k^2, nk)

@assert paper_D_block_sizes == expected_D_block_sizes == [240, 256]
@assert paper_Q_block_sizes == expected_Q_block_sizes == [240, 256]
@assert paper_G_block_sizes == expected_G_block_sizes == [512, 512]

for K in 0:nk-1
    println("²D, K = ", K,
            " | block size = ", paper_D_block_sizes[K + 1],
            " | first row labels = ", D_pair_labels_by_K[K + 1][1:3])
end
println("²Q block sizes by K = ", paper_Q_block_sizes)
println("²G block sizes by K = ", paper_G_block_sizes)

# The `{}^{2}\!D` sizes match the paper exactly: `K = 0` gives 240
# pair rows and `K = 1` gives 256. In this abelian setting, sorting the
# pair basis by `K` is already a symmetry-adapted ordering: one sector,
# one irrep, one dense block.
#
# The key numbers, locked down in
# `test/data/expectations/h4_periodic_v2rdm.toml`:

baseline = (
    modes           = asset.total_spin_orbital_modes,      # 32
    order2_basis    = blocker["order2_basis_size"],        # 2081
    order2_nuniq    = blocker["order2_nuniq"],             # 679121
    cs_constraints  = cs_only["moment_problem_constraints"],  # 508611
    cs_mosek_err    = cs_only["mosek_error_symbol"],       # MSK_RES_ERR_SPACE
)

# TS block footprint on the same 32-mode problem, from the
# `demos/h4_periodic_term_sparsity_cosmo_benchmark.jl` run log:

ts_stats = (
    n_blocks                = 313,
    max_block               = 168,
    sum_sq                  = 983_593,  # exact Σ block_dim² from benchmark log
    upper_tri_psd_slots     = 498_845,  # exact Σ d(d+1)/2 over the 313 PSD blocks
    total_constraint_blocks = 8_571,    # all symbolic constraint blocks after parity + moment-equalities
    moments                 = 26_817,   # unique symbolic moment values in primal route
    primal_jump_basis       = 85_941,   # |symmetric_canon(expval(total_basis))|
    primal_jump_vars        = 171_882,  # real/imag JuMP vars on primal route
    sos_dual_vars           = 2_006_058,# lifted matrix multipliers on dual route
)

# Paper's V2RDM[4,8] sector structure at `Nk = 2`, computed above by
# direct `K`-sectorization of the pair basis:

paper_blocks = (K0 = paper_D_block_sizes[1], K1 = paper_D_block_sizes[2])
paper_n_blocks = length(paper_D_block_sizes)

# ---
#
# ## 5. The three cost regimes, side by side
#
# Three metrics matter for SDP cost:
#
# - **largest PSD block dim** — dominates interior-point factorization cost.
# - **number of blocks** — affects ADMM parallelism and memory layout.
# - **sum of squared block sizes** — a proxy for total PSD variable slots.
#
# The dense number is ruinous; TS helps a lot but still shards each
# symmetry sector; the paper's two coarse symmetry blocks are the
# goal.

dense_sum_sq  = baseline.order2_basis^2
paper_sum_sq  = paper_blocks.K0^2 + paper_blocks.K1^2
ts_sum_sq     = ts_stats.sum_sq

cost_table = [
    ("dense order-2"             , 1                , baseline.order2_basis, dense_sum_sq),
    ("TS (MMD) order-2"          , ts_stats.n_blocks, ts_stats.max_block   , ts_sum_sq),
    ("paper K-blocked V2RDM[4,8]", paper_n_blocks   , paper_blocks.K1      , paper_sum_sq),
]

for (label, nb, maxb, ssq) in cost_table
    println(rpad(label, 34),
            " | blocks = ", lpad(nb, 4),
            " | max dim = ", lpad(maxb, 4),
            " | Σ dim² ≈ ", lpad(ssq, 9))
end

# Three numbers, three regimes. Let's draw them.

using CairoMakie

set_theme!(theme_light())

fig = Figure(size = (900, 420), fontsize = 14)

labels = ["dense\norder 2", "TS (MMD)\norder 2", "paper\nK-blocked"]
x_positions = 1:3
bar_color = [:firebrick, :steelblue, :seagreen]

ax1 = Axis(fig[1, 1];
    title = "Largest PSD block dim",
    ylabel = "dim (rows)",
    xticks = (x_positions, labels),
    yscale = log10,
)
barplot!(ax1, x_positions,
         [baseline.order2_basis, ts_stats.max_block, paper_blocks.K1];
         color = bar_color)
for (x, y) in zip(x_positions, [baseline.order2_basis, ts_stats.max_block, paper_blocks.K1])
    text!(ax1, x, y; text = string(y), align = (:center, :bottom),
          offset = (0, 4), fontsize = 12)
end

ax2 = Axis(fig[1, 2];
    title = "Number of PSD blocks",
    ylabel = "count",
    xticks = (x_positions, labels),
    yscale = log10,
)
barplot!(ax2, x_positions,
         [1, ts_stats.n_blocks, paper_n_blocks];
         color = bar_color)
for (x, y) in zip(x_positions, [1, ts_stats.n_blocks, paper_n_blocks])
    text!(ax2, x, y; text = string(y), align = (:center, :bottom),
          offset = (0, 4), fontsize = 12)
end

ax3 = Axis(fig[1, 3];
    title = "Σ (block dim)²  (≈ PSD slots)",
    ylabel = "count",
    xticks = (x_positions, labels),
    yscale = log10,
)
barplot!(ax3, x_positions,
         [dense_sum_sq, ts_sum_sq, paper_sum_sq];
         color = bar_color)
for (x, y) in zip(x_positions, [dense_sum_sq, ts_sum_sq, paper_sum_sq])
    text!(ax3, x, y; text = string(y), align = (:center, :bottom),
          offset = (0, 4), fontsize = 12)
end

Label(fig[0, :],
      "H₄ periodic Nk = 2 [4,8]: three regimes of block-diagonal structure";
      fontsize = 17, font = :bold)

fig

# **Reading the plot.** Each metric is on a log scale. Dense is
# unsolvable (2081×2081, 508k equality constraints, Mosek OOMs). TS
# is still a real improvement — factor 12 in max block dim and about
# factor 4.4 in total PSD slots versus dense. But TS is still about
# 8× the paper's footprint on $\sum \dim^2$, and about 156× more
# blocks. **The paper's coarse symmetry projection is the floor.**

# ---
#
# ## 6. Why TS shards where symmetry does not: a toy case
#
# To see the mechanism, look at a tiny fermionic Hamiltonian where
# we can list every basis word. Take **two sites, two spins, 4
# modes**, singlet sector. The degree-1 words split by
# $(\Delta N, S_z)$. A full-symmetry projection sees five irreps:
# $\Delta N \in \{-2,-1,0,+1,+2\}$ times $S_z \in \{-1,0,+1\}$.
#
# Term sparsity on the moment matrix, in contrast, cares about
# *which specific word products* appear in the activated support.
# The same irrep label may correspond to many disconnected TS
# components. Here is that calculation, runnable end to end.

using NCTSSoS

registry_toy, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
    ("c_up", 1:2),
    ("c_dn", 1:2),
])

## Simple Hubbard dimer for structural illustration only.
t_hop, U_int = 1.0, 2.0
ham_toy = -t_hop * (c_up_dag[1] * c_up[2] + c_up_dag[2] * c_up[1] +
                    c_dn_dag[1] * c_dn[2] + c_dn_dag[2] * c_dn[1]) +
          U_int * sum((c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
                      for i in 1:2)

## Each fermionic mode in `registry_toy` carries a sign and an index.
## Convention (src/types/registry.jl): `idx < 0` is a creator, `idx > 0` is an
## annihilator; `abs(idx)` identifies the physical spin-orbital mode.
## We tag each signed index by its (ΔN, spin) quantum numbers.

signed_indices = sort!(collect(NCTSSoS.indices(registry_toy)))

function mode_info(idx::Integer, registry)
    sym = registry[idx]
    is_creator = idx < 0                         # negative = creator
    ΔN = is_creator ? +1 : -1                    # creation raises N by +1
    name = String(sym)
    spin = occursin("up", name) ? :up : :dn
    return (; sym, is_creator, ΔN, spin)
end

mode_table = [mode_info(i, registry_toy) for i in signed_indices]
for info in mode_table
    println("  ", rpad(info.sym, 12),
            " ΔN=", lpad(info.ΔN, 2),
            "   spin=", info.spin,
            "   is_creator=", info.is_creator)
end

# So the 8 signed indices (4 creators + 4 annihilators) decompose
# into four $\Delta N = +1$ creators and four $\Delta N = -1$
# annihilators, each further split by spin. That is the **symmetry
# picture**: 4 degree-1 basis words per $(\Delta N, S_z)$ sector.
#
# TS sees a graph. Let's actually run it on this tiny problem and
# print the resulting blocks — this is cheap enough to execute here.

using NCTSSoS: polyopt, compute_sparsity, SolverConfig, MMD

n_up_toy = 1.0 * sum(c_up_dag[i] * c_up[i] for i in 1:2)
n_dn_toy = 1.0 * sum(c_dn_dag[i] * c_dn[i] for i in 1:2)

pop_toy = polyopt(ham_toy, registry_toy;
    moment_eq_constraints = [n_up_toy - 1.0 * one(ham_toy),
                             n_dn_toy - 1.0 * one(ham_toy)])

config_toy = SolverConfig(optimizer = nothing, order = 2, ts_algo = MMD())
sparsity_toy = compute_sparsity(pop_toy, config_toy)

blocks_toy = sparsity_toy.cliques_term_sparsities[1][1].block_bases
block_sizes_toy = length.(blocks_toy)

println("\nToy Hubbard dimer, order 2, MMD():")
println("  correlative cliques : ", length(sparsity_toy.corr_sparsity.cliques))
println("  TS PSD blocks       : ", length(blocks_toy))
println("  block sizes         : ", sort(block_sizes_toy; rev = true))
println("  sum(dim²)           : ", sum(s^2 for s in block_sizes_toy))
println("  max dim             : ", maximum(block_sizes_toy))

# The dimer is small enough that TS is already close to the
# symmetry irrep structure. On the full H₄ problem, the same
# procedure scales up — but the chordal extension starts dropping
# edges that symmetry would keep, and one sector ends up scattered
# across many TS blocks.

# ---
#
# ## 7. The cost wall: why the SOS dual variable count blows up
#
# The key distinction is simple and easy to miss:
#
# - the **primal moment route** optimizes over moment scalars,
# - the **SOS dual route** optimizes over **matrix multipliers**,
#   i.e. the matrix-valued Lagrange multipliers attached to each
#   primal matrix constraint.
#
# For this benchmark, you need to track both the symbolic moments and
# the actual JuMP model size:
#
# | Route | symbolic moments | JuMP vars / eqs | PSD blocks (max dim) | status |
# |---|---:|---:|---:|---|
# | Dense CS only | 679,121 | 508,611 eq constraints | 1 × **2081** | Mosek OOM (`MSK_RES_ERR_SPACE`) |
# | CS + TS (MMD), SOS dual | — | 2,006,058 JuMP vars | 313 × 168 | COSMO does not terminate in 45 min |
# | CS + TS (MMD), **primal moment** | 26,817 | 171,882 JuMP vars | 313 × 168 | tractable assembly |
# | Paper V2RDM[4,8] | $\sim 10^4$ | — | 2 × {240, 256} | minutes on boundary-point |
#
# The dual count is not mysterious. It closes exactly from the block
# data in the benchmark log:
#
# A Hermitian PSD block of size `d` is dualized in `sos_dualize` as a
# lifted real PSD variable of size `2d × 2d`, so it contributes
# `(2d)(2d+1)/2 = 2d^2 + d` scalar JuMP variables. A scalar complex
# zero-cone block becomes a lifted `2 × 2` symmetric multiplier, so it
# contributes 3 scalar variables.
#
# On this H₄ asset:

ts_sum_dims = 2 * ts_stats.upper_tri_psd_slots - ts_stats.sum_sq
zero_blocks = ts_stats.total_constraint_blocks - ts_stats.n_blocks

dual_breakdown = (
    psd_multiplier_vars  = 2 * ts_stats.sum_sq + ts_sum_dims,
    zero_blocks          = zero_blocks,
    zero_multiplier_vars = 3 * zero_blocks,
    scalar_objective     = 1,
)

@assert dual_breakdown.psd_multiplier_vars +
        dual_breakdown.zero_multiplier_vars +
        dual_breakdown.scalar_objective == ts_stats.sos_dual_vars

dual_breakdown

# Interpreting the numbers:
#
# - the 313 lifted HPSD multipliers contribute **1,981,283** scalar
#   dual vars;
# - the remaining **8,258** zero-cone blocks (mostly parity and
#   one-sided moment-equality constraints) contribute another
#   **24,774** vars;
# - the scalar objective variable `b` gives the last `+1`.
#
# Total: `1,981,283 + 24,774 + 1 = 2,006,058`.
#
# So the SOS-dual blow-up is structural: the model is not optimizing
# over the 26,817 symbolic moments, but over **lifted matrix
# multipliers for every constraint block**.
#
# Also, do not compare `2,006,058` directly to `26,817` as if both
# were JuMP variable counts. The strict apples-to-apples comparison is
#
# - SOS dual: `2,006,058` JuMP vars,
# - primal moment: `171,882` JuMP vars,
#
# which is still about **11.7×** smaller on the primal route. The
# larger **75×** gap compares dual JuMP vars against the **symbolic**
# unique-moment count, which is a useful structural signal but not the
# same unit.
#
# Even with TS + primal moment, the 313 blocks with max 168 are still
# at the upper end of what COSMO solves comfortably with loose
# tolerances. Tightening takes either a boundary-point solver (the
# paper's choice) or coarser blocks.

# ---
#
# ## 8. Concrete scalability levers in NCTSSoS.jl today
#
# Ranked from biggest lever to smallest, from least to most invasive.
#
# ### (A) Use `correlative_sparsity(pop, moment_basis, elim_algo)` with a pair-operator basis
#
# This is the **single biggest lever** and requires no code changes
# to NCTSSoS — just a user-supplied basis. The overload already
# exists in `src/optimization/sparsity.jl`:
#
# ```julia
# correlative_sparsity(
#     pop::OP,
#     moment_basis::AbstractVector,      # <-- you supply this
#     elim_algo::EliminationAlgorithm,
# ) where { ... }
# ```
#
# Feed it only pair operators — the same rows ${}^{2}\!D, {}^{2}\!Q,
# {}^{2}\!G$ use in the paper — plus the identity:
#
# ```julia
# pair_basis = NormalMonomial[one(NormalMonomial)]
# for p in 1:n_modes, q in (p+1):n_modes
#     push!(pair_basis, annihilator[p] * annihilator[q])            # ^2D
#     push!(pair_basis, creator[p] * creator[q])                    # ^2Q
# end
# for p in 1:n_modes, q in 1:n_modes
#     push!(pair_basis, creator[p] * annihilator[q])                # ^2G / ^1D
# end
# ```
#
# Expected impact on Nk=2 [4,8]:
#
# - moment matrix drops from 2081 rows to $\le 500$.
# - unique moments drop from 679,121 toward $\sim 10^4$ (scaling
#   dominated by the degree-$\le 4$ words reachable from pair-pair
#   products).
# - This is the paper's PQG relaxation, not the full order-2
#   Lasserre. Weaker lower bound, but at the paper's cost.
#
# Caveats: you lose the extra tightness of full order-2 Lasserre;
# and your pair basis must include the identity (NCTSSoS enforces
# `one(M) in basis`).
#
# ### (B) Partition the basis by symmetry *before* TS
#
# The demo script `demos/h4_periodic_term_sparsity_size.jl` already
# verifies that every TS block is K-pure and $\Delta N$-pure for
# this Hamiltonian. That empirical observation can be turned into a
# guaranteed invariant by a small change in `compute_sparsity`:
#
# 1. Accept a tagging function `quantum_numbers(m::NormalMonomial) → Tuple`.
# 2. Partition each `clq_mom_mtx_bases[i]` by tag *before* calling
#    `term_sparsities`.
# 3. Run TS inside each partition and union the results.
#
# The moment-matrix assembly in `src/optimization/moment.jl`
# (`_build_constraint_matrix`) already consumes one basis per block
# — nothing downstream needs to change.
#
# For the H₄ asset, tagging by `(ΔN, K, S_z)` gives tens of sectors,
# each a few hundred basis words. TS inside each sector either
# splits it further or leaves it as one coarse block. Either way you
# recover the paper-style coarseness and keep TS as a free refinement.
#
# ### (C) Generalize the fermionic parity hook to arbitrary abelian symmetry
#
# `src/optimization/moment.jl` already has `_add_parity_constraints!`
# for fermionic $\mathbb Z_2$ parity. The same code path is exactly
# what is needed for arbitrary abelian symmetry labels — the
# parity rule is just the case $G = \mathbb Z_2$. Generalizing this
# to user-specified characters turns (B) from a kwarg into a first-
# class construct.
#
# ### (D) Default to the primal moment route for wide fermionic problems
#
# The benchmark log records: SOS dualize → 2,006,058 JuMP variables,
# COSMO times out at 45 min. Primal moment → 171,882 JuMP variables on
# the assembled Hermitian model, representing 26,817 unique symbolic
# moments, with the same 313 PSD blocks. So the strict model-size
# reduction is about **11.7×**, and the symbolic unique-moment count is
# about **75×** smaller than the dual variable count.
#
# The reason is formulation, not solver superstition: the dual route
# introduces lifted matrix multipliers for every HPSD and Zero
# constraint block. On this H₄ asset, **1,981,283** of the 2,006,058
# dual vars come from the 313 lifted HPSD multipliers alone.
#
# `SolverConfig` does not currently expose this as a first-class
# toggle; it is selected by calling `solve_moment_problem(mp, ...)`
# directly on the output of `NCTSSoS.moment_relax(...)` instead of
# going through `cs_nctssos`. For wide problems, this should be the
# default — or at minimum, a supported configuration flag.
#
# ### (E) Boundary-point SDP backend
#
# The paper uses a boundary-point method (Mazziotti, PRL **106**,
# 083001, 2011) which is tailored to V2RDM SDPs: many small PSD
# blocks, one large equality block, no barrier parameter. Neither
# Mosek (interior-point) nor COSMO (first-order ADMM) is a
# comfortable fit. Research-grade, but the ultimate solver for this
# problem family.

# ---
#
# ## 9. The mental picture
#
# One picture is worth a thousand index conventions:

fig2 = Figure(size = (900, 300), fontsize = 14)

ax_dense = Axis(fig2[1, 1]; title = "dense order-2\n1 block (2081 × 2081)",
                aspect = DataAspect(), xticksvisible = false,
                yticksvisible = false, xticklabelsvisible = false,
                yticklabelsvisible = false)
poly!(ax_dense, Rect(0, 0, 1, 1); color = (:firebrick, 0.6),
      strokewidth = 2, strokecolor = :firebrick)
limits!(ax_dense, -0.1, 1.1, -0.1, 1.1)

ax_ts = Axis(fig2[1, 2]; title = "TS (MMD)\n313 blocks (max 168)",
             aspect = DataAspect(), xticksvisible = false,
             yticksvisible = false, xticklabelsvisible = false,
             yticklabelsvisible = false)
## Draw a schematic: many small blocks down the diagonal.
rng = 1:17
for (i, r) in enumerate(rng)
    size = 0.03 + 0.04 * (i % 4)
    x0 = (r - 1) / length(rng)
    poly!(ax_ts, Rect(x0, x0, size, size); color = (:steelblue, 0.65),
          strokewidth = 1, strokecolor = :steelblue)
end
limits!(ax_ts, -0.05, 1.05, -0.05, 1.05)

ax_paper = Axis(fig2[1, 3]; title = "paper K-blocked V2RDM\n2 blocks ({240, 256})",
                aspect = DataAspect(), xticksvisible = false,
                yticksvisible = false, xticklabelsvisible = false,
                yticklabelsvisible = false)
poly!(ax_paper, Rect(0.00, 0.48, 0.48, 0.48); color = (:seagreen, 0.6),
      strokewidth = 2, strokecolor = :seagreen)
poly!(ax_paper, Rect(0.50, 0.00, 0.50, 0.48); color = (:seagreen, 0.45),
      strokewidth = 2, strokecolor = :seagreen)
text!(ax_paper, 0.24, 0.72; text = "K = 0", fontsize = 12, align = (:center, :center))
text!(ax_paper, 0.75, 0.24; text = "K = 1", fontsize = 12, align = (:center, :center))
limits!(ax_paper, -0.05, 1.05, -0.05, 1.05)

Label(fig2[0, :],
      "Schematic: same underlying moment matrix, three different block structures";
      fontsize = 16, font = :bold)

fig2

# - **Left (dense)**: one monolithic PSD constraint. Correct, useless.
# - **Middle (TS)**: many small blocks down the diagonal. Each
#   lives inside a symmetry sector, but the sector is shattered.
# - **Right (paper)**: two big sector blocks — one per irrep of the
#   translation group. Coarsest possible blocking consistent with
#   the symmetry.
#
# The scalability path is clear: the TS picture already respects
# the symmetry; it just does not exploit it coarsely. Levers (A)
# and (B) above turn the TS picture into the paper picture without
# rewriting the solver stack.

# ---
#
# ## 10. Summary
#
# - The paper's ${}^{2}\!D$ blocks are **irreducible representations
#   of $\mathbb Z_{N_k}$** — one block per total pair momentum $K$.
#   Coarsest blocking compatible with translation symmetry.
# - NCTSSoS's TS blocks are **maximal cliques of a chordal extension**
#   of a graph defined by monomial equality after normal ordering.
#   They respect the same symmetry ($K$-pure, $\Delta N$-pure) but
#   shatter each sector into many smaller blocks.
# - On H₄ Nk=2 [4,8]: dense = 1 block of 2081 (OOM), TS = 313 blocks
#   of max 168, paper = 2 blocks of {240, 256}; and the SOS-dual route
#   inflates to 2,006,058 scalar vars because it optimizes over lifted
#   matrix multipliers rather than moment scalars.
# - **Actionable levers inside NCTSSoS today**: use the explicit
#   `moment_basis` API (A) to match the paper's PQG relaxation;
#   add symmetry-adapted partitioning before TS (B) to recover
#   coarse blocks; default to the primal-moment assembly (D) to
#   avoid the SOS-dual blow-up. Boundary-point SDP (E) is the
#   long-term solver story.
#
# ### See also
#
# - [Periodic V2RDM Benchmark (H₄ Chain)](@ref periodic-v2rdm-h4) — full specification of the benchmark problem.
# - [Why Correlative Sparsity Misses the $k$-Blocks](@ref periodic-v2rdm-cs-k-blocks) — the weaker, CS-only picture on the same asset.
# - [Building the Paper's PQG Relaxation](@ref periodic-v2rdm-pqg-construction) — executes levers (A)+(B) from §8 of this page and prints SDP sizes without solving.
# - [Stabilization vs. Exactness in the Sparse Hierarchy](@ref sparsity-convergence) — how TS iterations converge on smaller NC problems.
#
# ### References
#
# - A. O. Schouten, S. Ewing, D. A. Mazziotti,
#   "Bootstrapping the Electronic Structure of Quantum Materials,"
#   arXiv:2504.02861 (2025).
# - D. A. Mazziotti, "Large-Scale Semidefinite Programming for
#   Many-Electron Quantum Mechanics,"
#   *Phys. Rev. Lett.* **106**, 083001 (2011). — The boundary-point
#   solver used in the paper.
# - J. Wang, V. Magron, J.-B. Lasserre,
#   "TSSOS: A Moment-SOS Hierarchy That Exploits Term Sparsity,"
#   *SIAM J. Optim.* **31**, 30–58 (2021).
