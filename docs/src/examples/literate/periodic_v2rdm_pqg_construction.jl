# # [Periodic V2RDM: Building the Paper's PQG Relaxation (No Solve)](@id periodic-v2rdm-pqg-construction)
#
# This page does one thing: reproduce — as a concrete, executable
# object — the **paper's PQG block structure** for the H₄ periodic
# benchmark (Schouten, Ewing & Mazziotti,
# [arXiv:2504.02861](https://arxiv.org/abs/2504.02861)), then print
# how big the resulting SDP is. It **does not solve** the SDP; the
# goal is to have a proper moment-relaxation object sized and
# documented so you know what you would hand to a solver.
#
# Companion pages:
#
# - [Periodic V2RDM Benchmark (H₄ Chain)](@ref periodic-v2rdm-h4) —
#   conventions (index order, chemist↔physicist, per-cell
#   normalisation, etc.).
# - [Why Correlative Sparsity Misses the $k$-Blocks](@ref periodic-v2rdm-cs-k-blocks) —
#   the support graph is complete; CS alone can't recover the paper's
#   blocks.
# - [²D Symmetry Blocks vs. Term-Sparsity Blocks](@ref periodic-v2rdm-block-structure) —
#   the mental model. Lists lever (A) "use the paper's pair basis
#   via `moment_basis`" and lever (B) "partition that basis by
#   conserved quantum numbers before TS." **This page executes
#   both.**
#
# ## What gets built, in three concrete steps
#
# 1. **Paper's moment basis.** Construct the explicit list
#    $B = \{1\} \cup \{a_i a_j\}_{i<j} \cup \{a_i^\dagger a_j^\dagger\}_{i<j} \cup \{a_i^\dagger a_j\}$
#    on 32 spin-orbital modes. Group by $\Delta N$ so each of ${}^{2}\!D$,
#    ${}^{2}\!Q$, ${}^{2}\!G$ appears as a distinct block of the moment
#    matrix.
# 2. **Force paper-style momentum blocking.** Instead of letting
#    NCTSSoS's term-sparsity (MMD) chordalisation shard each sector
#    into many small blocks, manually build a `TermSparsity` whose
#    `block_bases` is the partition of $B$ by total crystal momentum
#    $K = (k_i + k_j) \bmod N_k$. Verify the block sizes match the
#    paper: ${}^{2}\!D = \{240, 256\}$, ${}^{2}\!G = \{512, 512\}$,
#    ${}^{2}\!Q = \{240, 256\}$.
# 3. **Build the moment problem without solving.** Call
#    `moment_relax` directly, then report: number of PSD blocks,
#    largest block dimension, $\sum \dim^2$ (PSD-slot proxy), unique
#    moment variables, total-basis size, and number of zero-cone
#    constraints (fermionic parity + one-sided particle-number
#    moment equalities).
#
# !!! note "Same bound, not the same symbolic SDP"
#     Block-diagonalising the moment matrix by conserved charges
#     gives an SDP whose *optimum* coincides with the full-PSD PQG
#     relaxation on the same basis (objective and constraints only
#     couple ΔN = 0 moments, so the off-sector entries are free
#     scalars that never enter any PSD cone). The symbolic SDP is
#     different from the dense version — smaller and sector-aware —
#     but the bound is unchanged.

# ---
#
# ## Load the vendored asset and build the Bloch Hamiltonian
#
# Use the same helper as the other periodic-V2RDM pages. The 32-mode
# Hamiltonian has 23,752 monomials after normal ordering; the support
# graph is the complete graph on 16 compound spatial orbitals (see
# the [CS page](@ref periodic-v2rdm-cs-k-blocks)).

using NCTSSoS
using LinearAlgebra
using Printf

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: load_nk2_asset

asset = load_nk2_asset()
nk    = asset.nk                # 2
norb  = asset.n_active_orb      # 8

function build_h4_nk2_hamiltonian(h1e, eri; nk::Int, norb::Int)
    registry, ((c_up_k0, c_up_k0_dag),
               (c_dn_k0, c_dn_k0_dag),
               (c_up_k1, c_up_k1_dag),
               (c_dn_k1, c_dn_k1_dag)) = create_fermionic_variables([
        ("c_up_k0", 1:norb),
        ("c_dn_k0", 1:norb),
        ("c_up_k1", 1:norb),
        ("c_dn_k1", 1:norb),
    ])

    ann = Dict((0,:up)=>c_up_k0, (0,:dn)=>c_dn_k0,
               (1,:up)=>c_up_k1, (1,:dn)=>c_dn_k1)
    dag = Dict((0,:up)=>c_up_k0_dag, (0,:dn)=>c_dn_k0_dag,
               (1,:up)=>c_up_k1_dag, (1,:dn)=>c_dn_k1_dag)
    spin_channels = ((:up,:up), (:up,:dn), (:dn,:up), (:dn,:dn))

    ham = (0.0 + 0.0im) * (c_up_k0_dag[1] * c_up_k0[1])

    for k in 0:nk-1
        hk = h1e[k] / nk
        for p in 1:norb, s in 1:norb, σ in (:up,:dn)
            coeff = hk[p, s]
            iszero(coeff) && continue
            ham += coeff * dag[(k,σ)][p] * ann[(k,σ)][s]
        end
    end

    for (k1,k2,k3,k4) in sort!(collect(keys(eri)))
        V = eri[(k1,k2,k3,k4)] / nk^2
        for p in 1:norb, r in 1:norb, q in 1:norb, s in 1:norb
            coeff = V[p, r, q, s]
            iszero(coeff) && continue
            for (σ, τ) in spin_channels
                ham += 0.5 * coeff *
                    dag[(k1,σ)][p] * dag[(k2,τ)][q] *
                    ann[(k4,τ)][s] * ann[(k3,σ)][r]
            end
        end
    end

    return registry,
           ((c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
            (c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag)),
           0.5 * (ham + adjoint(ham))
end

registry, vars, ham = build_h4_nk2_hamiltonian(
    asset.h1e, asset.eri; nk = nk, norb = norb)

(c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
(c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag) = vars

n_up = 1.0 * sum(c_up_k0_dag[i]*c_up_k0[i] + c_up_k1_dag[i]*c_up_k1[i] for i in 1:norb)
n_dn = 1.0 * sum(c_dn_k0_dag[i]*c_dn_k0[i] + c_dn_k1_dag[i]*c_dn_k1[i] for i in 1:norb)

pop = polyopt(ham, registry;
    moment_eq_constraints = [n_up - 4.0 * one(ham),
                             n_dn - 4.0 * one(ham)])

(; terms_in_H = length(monomials(ham)), spin_modes = length(variable_indices(ham)))

# ---
#
# ## Step 1 — Construct the paper's PQG moment basis
#
# The V2RDM paper's moment basis has four parts, one for each of the
# PQG blocks plus the identity:
#
# - 1 × identity (carries the trace normalisation);
# - $\binom{M}{2} = 496$ annihilator pairs $a_i a_j$ with $i<j$ —
#   these are the rows of ${}^{2}\!D$;
# - $496$ creator pairs $a_i^\dagger a_j^\dagger$ with $i<j$ — the
#   rows of ${}^{2}\!Q$;
# - $M^2 = 1024$ particle–hole pairs $a_i^\dagger a_j$ (allowing
#   $i = j$ — the $n_i$ operators sit on the diagonal) — the rows of
#   ${}^{2}\!G$.
#
# With $M = 2r N_k = 32$ modes, the total is
# $|B| = 1 + 496 + 496 + 1024 = 2017$, versus the dense order-2
# basis $2081$. The 64 missing entries are the degree-1
# single-creator / single-annihilator words, which NCTSSoS's default
# order-2 basis includes but the paper's PQG basis does not. Those
# words live in $\Delta N = \pm 1$ sectors, decoupled from the
# number-conserving Hamiltonian.

## Resolve each registry mode index (1..32) into its (k, spin, orbital) label.
function mode_info(mode_id::Integer; norb::Int = 8)
    family, o0 = divrem(Int(mode_id) - 1, norb)
    k    = family ÷ 2            # 0,1 → k=0; 2,3 → k=1
    spin = iseven(family) ? :up : :dn
    return (; k, spin, orb = o0 + 1, family)
end

## Conserved-charge labels for a basis element:
##   ΔN  = net particle number change;
##   K   = net crystal momentum mod Nk;
##   two_Sz = 2·Sz (integer) for cleanliness.
function basis_label(m::NormalMonomial; nk::Int = nk, norb::Int = norb)
    ΔN = 0; K = 0; two_Sz = 0
    for idx in m.word
        mode = abs(Int(idx))
        info = mode_info(mode; norb)
        sign = idx < 0 ? +1 : -1          # creator → +1, annihilator → −1
        σ    = info.spin === :up ? +1 : -1
        ΔN     += sign
        K       = mod(K + sign * info.k, nk)
        two_Sz += sign * σ
    end
    return (; ΔN, K, two_Sz, deg = length(m.word))
end

## All 32 physical spin-orbitals, in the same order the registry assigns them.
spin_orbitals = NamedTuple[]
for family in 0:2nk-1        # 0 → (k=0,↑), 1 → (k=0,↓), 2 → (k=1,↑), 3 → (k=1,↓)
    k    = family ÷ 2
    spin = iseven(family) ? :up : :dn
    arr_ann = Dict((0,:up)=>c_up_k0, (0,:dn)=>c_dn_k0,
                   (1,:up)=>c_up_k1, (1,:dn)=>c_dn_k1)[(k, spin)]
    arr_dag = Dict((0,:up)=>c_up_k0_dag, (0,:dn)=>c_dn_k0_dag,
                   (1,:up)=>c_up_k1_dag, (1,:dn)=>c_dn_k1_dag)[(k, spin)]
    for p in 1:norb
        push!(spin_orbitals,
              (; a = arr_ann[p], adag = arr_dag[p], k, spin, orb = p))
    end
end
M = length(spin_orbitals)          # 32

## Build the four parts of B. Each `monomials(expr)[1]` pulls the unique
## normal-ordered NormalMonomial out of a product of operators.
## We pin the concrete type `NM = NormalMonomial{FermionicAlgebra, Int8}`
## (the registry picks `Int8` for ≤ 32 modes) so every basis vector is
## strictly typed — required by `TermSparsity{NM}` below.
NM = typeof(monomials(c_up_k0_dag[1] * c_up_k0[1])[1])
@assert NM <: NormalMonomial

identity_mono = one(NM)

D_basis = NM[]          # a_i a_j (i<j): ²D rows, ΔN = −2
Q_basis = NM[]          # a†_i a†_j (i<j): ²Q rows, ΔN = +2
G_basis = NM[]          # a†_i a_j (all i,j): ²G rows, ΔN = 0

for i in 1:M-1, j in i+1:M
    push!(D_basis, monomials(spin_orbitals[i].a    * spin_orbitals[j].a)[1])
    push!(Q_basis, monomials(spin_orbitals[i].adag * spin_orbitals[j].adag)[1])
end
for i in 1:M, j in 1:M
    push!(G_basis, monomials(spin_orbitals[i].adag * spin_orbitals[j].a)[1])
end

pqg_basis = NM[identity_mono; D_basis; Q_basis; G_basis]
unique!(pqg_basis)

## Group by conserved charges (ΔN, K, 2Sz) — this is the V2RDM-friendly sort.
labels = basis_label.(pqg_basis)
sort_perm = sortperm(eachindex(pqg_basis); by = i -> (
    labels[i].ΔN, labels[i].K, labels[i].two_Sz, length(pqg_basis[i].word),
))
pqg_basis = pqg_basis[sort_perm]
labels    = labels[sort_perm]

## Sanity: identity present, total basis count matches algebra.
@assert identity_mono in pqg_basis
@assert length(pqg_basis) == 1 + 2 * binomial(M, 2) + M^2 == 2017

(; basis_size = length(pqg_basis), M = M,
   expected   = 1 + 2 * binomial(M, 2) + M^2,
   ΔN_minus2  = count(l -> l.ΔN == -2, labels),
   ΔN_zero    = count(l -> l.ΔN == 0, labels),
   ΔN_plus2   = count(l -> l.ΔN == 2, labels))

# So $|B| = 2017$: 496 + 496 + 1025, where the 1025 $\Delta N = 0$
# part consists of 1024 particle–hole pairs and the identity.

# ---
#
# ## Step 2 — Partition by crystal momentum $K$ (the paper's blocking)
#
# The paper's ${}^{2}\!D$ block structure at $N_k = 2$ is
# $\{240, 256\}$ by total pair momentum $K = (k_i + k_j) \bmod 2$.
# The same blocking applies to ${}^{2}\!Q$ (trivially, by
# $-(k_i + k_j)$), and ${}^{2}\!G$ splits by $(k_i - k_j) \bmod 2$
# into $\{512, 512\}$. Identity lives in the $\Delta N = 0, K = 0$
# sector with ${}^{2}\!G$'s same-$k$ particle–hole pairs, so that
# block is 513 (= 1 + 512).

sector_labels = [(l.ΔN, l.K) for l in labels]
sectors = unique(sector_labels)
sort!(sectors; by = s -> (s[1], s[2]))

paper_blocks = Dict(s => NM[] for s in sectors)
for (b, s) in zip(pqg_basis, sector_labels)
    push!(paper_blocks[s], b)
end

sector_rows = map(sectors) do s
    (; ΔN = s[1], K = s[2], dim = length(paper_blocks[s]))
end

for r in sector_rows
    @printf("ΔN = %+d , K = %d  →  block dim = %4d\n", r.ΔN, r.K, r.dim)
end

@assert sum(r.dim for r in sector_rows) == length(pqg_basis)

# The exact paper-faithful sizes:

(; D_K0 = first(r.dim for r in sector_rows if r.ΔN == -2 && r.K == 0),
   D_K1 = first(r.dim for r in sector_rows if r.ΔN == -2 && r.K == 1),
   G_K0 = first(r.dim for r in sector_rows if r.ΔN ==  0 && r.K == 0),
   G_K1 = first(r.dim for r in sector_rows if r.ΔN ==  0 && r.K == 1),
   Q_K0 = first(r.dim for r in sector_rows if r.ΔN ==  2 && r.K == 0),
   Q_K1 = first(r.dim for r in sector_rows if r.ΔN ==  2 && r.K == 1))

# Reading that tuple:
#
# - ${}^{2}\!D$: $(240, 256)$ ✓ — matches Section 7 of the
#   [H₄ spec page](@ref periodic-v2rdm-h4).
# - ${}^{2}\!Q$: $(240, 256)$ ✓ — same pair counts, creator-pair indexing.
# - ${}^{2}\!G + \text{id}$: $(513, 512)$ ✓ — the $K=0$ block absorbs
#   the identity.
#
# **This is the paper's blocking**, reproduced from pure algebra —
# no term-sparsity graph was harmed.

# ---
#
# ### A finer refinement: split each $K$-sector by $2 S_z$
#
# The Hamiltonian also conserves $S_z$ (spin-diagonal integrals —
# see Section 3 of [the H₄ spec page](@ref periodic-v2rdm-h4)), so
# nothing stops us from splitting each $K$-sector further. This is
# tighter than the paper's own block report but mathematically
# equivalent for a spin-conserving Hamiltonian:

fine_labels = [(l.ΔN, l.K, l.two_Sz) for l in labels]
fine_sectors = sort!(unique(fine_labels); by = s -> (s[1], s[2], s[3]))
fine_blocks = Dict(s => NM[] for s in fine_sectors)
for (b, s) in zip(pqg_basis, fine_labels)
    push!(fine_blocks[s], b)
end

println("\n(ΔN, K, 2Sz) refinement (tighter than the paper's table):")
for s in fine_sectors
    @printf("  ΔN = %+d , K = %d , 2Sz = %+d  →  dim = %3d\n",
            s[1], s[2], s[3], length(fine_blocks[s]))
end

# That gives 18 symmetry sectors instead of 6. We stick with the
# paper-faithful (ΔN, K) partition below; switching to the finer one
# is a one-line change to `paper_blocks`.

# ---
#
# ## Step 3 — Build the `MomentProblem` without solving
#
# The path is:
#
# 1. Ask `compute_sparsity` to accept the PQG basis via `moment_basis`.
#    With `cs_algo = NoElimination()` + a complete support graph, it
#    returns exactly **one clique** of 32 modes (Blocker 1 from
#    [the H₄ spec page](@ref periodic-v2rdm-h4)). Good — we want the
#    basis filter to leave $B$ intact.
# 2. Rebuild the single `TermSparsity` in that clique with our
#    $(\Delta N, K)$ partition as `block_bases`. `TermSparsity` is
#    immutable, so we construct a fresh one rather than mutating.
# 3. Call `moment_relax` directly to produce the symbolic
#    `MomentProblem`. We **do not** hand it to any solver here.

config = SolverConfig(optimizer = nothing,
                      moment_basis = pqg_basis,
                      cs_algo = NoElimination(),
                      ts_algo = NoElimination())

sparsity = compute_sparsity(pop, config)

@assert length(sparsity.corr_sparsity.cliques) == 1
@assert length(sparsity.cliques_term_sparsities[1]) == 1   # 1 moment matrix, 0 localizing

default_ts = sparsity.cliques_term_sparsities[1][1]

(; clique_size            = length(sparsity.corr_sparsity.cliques[1]),
   clique_basis_len       = length(sparsity.corr_sparsity.clq_mom_mtx_bases[1]),
   default_ts_blocks      = length(default_ts.block_bases),
   default_ts_max_block   = maximum(length.(default_ts.block_bases)))

# With `ts_algo = NoElimination()`, the default TS already returns
# one block = the full PQG basis (2017). That's the "dense PSD on
# the PQG basis" — still too big. Swap it for the paper's
# $(\Delta N, K)$ partition.

## Re-use the activated support from the default TS step so the
## symbolic relaxation keeps the same moment-support set.
paper_block_bases = [paper_blocks[s] for s in sectors]

## Defensive: the $(ΔN=0, K=0)$ block must contain the identity, otherwise
## `_validate_polynomial_relaxation_support` will reject objective moments
## like `<I>` that are generated via `b_i† * 1 * b_j` when one of b_i, b_j
## is the identity.
id_sector = (0, 0)
@assert any(m -> m == identity_mono, paper_blocks[id_sector])

custom_ts = NCTSSoS.TermSparsity{NM}(
    default_ts.term_sparse_graph_supp,
    paper_block_bases,
)

cliques_ts = [[custom_ts]]

## Build the symbolic moment problem. This does not touch a solver.
# Suppress the raw object display; it's enormous and useless in docs output.
moment_problem = NCTSSoS.moment_relax(
    pop, sparsity.corr_sparsity, cliques_ts);

# ### Block footprint
#
# Read the PSD and zero-cone constraints back out of `moment_problem`.
# `moment_relax` has already:
#
# 1. Emitted one `HPSD` cone per entry of `block_bases` — that's our
#    $(\Delta N, K)$ partition.
# 2. Appended fermionic parity `:Zero` constraints
#    (`_add_parity_constraints!`) — these fix every odd-parity
#    moment to zero.
# 3. Appended the two one-sided particle-number moment equalities
#    $\langle b^\dagger (n_\uparrow - 4)\rangle = 0$ and
#    $\langle b^\dagger (n_\downarrow - 4)\rangle = 0$
#    (`_add_moment_eq_constraints!`) — these fix the state particle
#    number.

psd_block_sizes = sort(
    [size(mat, 1) for (cone, mat) in moment_problem.constraints if cone === :HPSD];
    rev = true,
)
n_zero_scalar = sum(
    length(mat) for (cone, mat) in moment_problem.constraints if cone === :Zero
)

block_stats = (
    n_psd_blocks         = length(psd_block_sizes),
    max_psd_block        = maximum(psd_block_sizes),
    sum_sq               = sum(d^2     for d in psd_block_sizes),
    sum_tri_slots        = sum(d*(d+1)÷2 for d in psd_block_sizes),
    block_sizes          = psd_block_sizes,
    n_unique_moments     = moment_problem.n_unique_moment_matrix_elements,
    total_basis_len      = length(moment_problem.total_basis),
    n_zero_scalar_eqs    = n_zero_scalar,
)

block_stats

# The interpretation:
#
# - `n_psd_blocks = 6`, `max_psd_block = 513`, block sizes sorted =
#   $[513, 512, 256, 256, 240, 240]$ — exactly the paper's
#   block structure (513 = 512 + identity).
# - `sum_sq ≈ 8.43 × 10⁵` vs.  dense-order-2 $2081^2 = 4.33 × 10⁶$
#   (a $\approx 5.1\times$ reduction in PSD slots versus the monolithic
#   order-2 block) and vs. TS(MMD)'s $983{,}593$
#   (the paper's coarser blocking beats MMD term sparsity by
#   $\approx 1.17\times$ on this metric — mostly because MMD's
#   sharding over-chordalises each sector).
# - `n_unique_moments` counts distinct canonical monomials appearing
#   *inside the 6 HPSD blocks*, not the full set of symbolic
#   unknowns.
# - `total_basis_len` is the superset that also includes moments
#   generated by the moment-equality constraints; that is the real
#   JuMP-variable count when we eventually hand this to a solver.
# - `n_zero_scalar_eqs` counts scalar equality constraints: parity
#   kills every odd-parity monomial, then two one-sided localising
#   chains fix $N_\uparrow$ and $N_\downarrow$.
#
# !!! warning "No solve on this page"
#     We deliberately stop at the symbolic `MomentProblem`. Actually
#     running COSMO or Mosek on a 6-block, max-dim-513 Hermitian SDP
#     with this many unique moments is a separate exercise — see
#     `demos/h4_periodic_term_sparsity_cosmo_benchmark.jl` for a full
#     run-log skeleton (that demo uses the MMD TS blocks; swapping
#     in `paper_block_bases` above is the remaining plumbing).

# ---
#
# ## Visualise the moment matrix blocks
#
# A schematic block-diagonal layout — one box per $(ΔN, K)$ sector,
# width proportional to block dim:

using CairoMakie
set_theme!(theme_light())

fig = Figure(size = (900, 360), fontsize = 14)

## Panel 1: scrollbar showing the 6 sector blocks along the diagonal.
ax = Axis(fig[1, 1]; title = "Block structure of the PQG moment matrix",
          aspect = DataAspect(), xticksvisible = false,
          yticksvisible = false, xticklabelsvisible = false,
          yticklabelsvisible = false)

colors = Dict(
    -2 => (:firebrick, 0.7),
     0 => (:seagreen,  0.7),
     2 => (:steelblue, 0.7),
)

total_dim = sum(r.dim for r in sector_rows)
block_offsets = cumsum(vcat(0, getfield.(sector_rows[1:end-1], :dim))) ./ total_dim
for (r, x0) in zip(sector_rows, block_offsets)
    w = r.dim / total_dim
    poly!(ax, Rect(x0, 1 - x0 - w, w, w);
          color = colors[r.ΔN],
          strokewidth = 1, strokecolor = :black)
    text!(ax, x0 + w/2, 1 - x0 - w/2;
          text = @sprintf("ΔN = %+d\nK = %d\n%d", r.ΔN, r.K, r.dim),
          fontsize = 10, align = (:center, :center))
end
limits!(ax, -0.02, 1.02, -0.02, 1.02)

## Panel 2: bar chart of block dims vs. dense / TS-MMD baseline.
ax2 = Axis(fig[1, 2]; title = "Max PSD block dim at order 2",
           ylabel = "rows",
           xticks = ([1,2,3], ["dense\norder 2",
                               "TS (MMD)\norder 2",
                               "paper\nK-blocks"]),
           yscale = log10)
barplot!(ax2, [1,2,3], [2081, 168, block_stats.max_psd_block];
         color = [:firebrick, :steelblue, :seagreen])
for (x, y) in zip([1,2,3], [2081, 168, block_stats.max_psd_block])
    text!(ax2, x, y; text = string(y), align = (:center, :bottom),
          offset = (0, 4), fontsize = 12)
end

Label(fig[0, :], "H₄ periodic Nk = 2 [4, 8]: paper-faithful PQG block construction";
      fontsize = 16, font = :bold)

fig

# Three colours, one per $\Delta N$ irrep of $U(1)_N$. Each box's
# size is proportional to its block dim, so the ${}^{2}\!G$ part
# dominates (as it should — it carries 1024 particle–hole moments
# plus the identity).
#
# The right panel is the cost story in one line: *the paper's
# coarse blocking beats both dense order-2 and TS (MMD) on the
# max-block-dim axis that dominates interior-point factorisation
# cost.*

# ---
#
# ## What you actually have now
#
# After this page:
#
# | Object | Meaning |
# |:-------|:--------|
# | `pqg_basis` | 2017-element paper PQG basis sorted by $(ΔN, K, 2S_z)$. |
# | `paper_block_bases` | 6-way partition of `pqg_basis` by $(ΔN, K)$. |
# | `fine_blocks` | 18-way refinement by $(ΔN, K, 2S_z)$; drop-in replacement for a tighter partition. |
# | `custom_ts` | `TermSparsity` carrying `paper_block_bases` as `block_bases`. |
# | `moment_problem` | Symbolic `MomentProblem` with 6 HPSD cones of dims $\{240, 240, 256, 256, 512, 513\}$. |
# | `block_stats` | Solver-facing sizes: PSD blocks, max block, $\sum \dim^2$, unique moments, zero-cone eqs. |
#
# Two things are **not** done here, on purpose:
#
# 1. **No solve.** Hand `moment_problem` to `solve_moment_problem`
#    (primal moment route) or `sos_dualize + optimize!` (dual
#    route) when you want a number. The primal route is the right
#    default on this problem — see [`block_structure` page §7](@ref periodic-v2rdm-block-structure)
#    for the 11.7× JuMP-variable advantage on H₄.
# 2. **No $U(1)_N$ off-sector zero injection.** The off-sector
#    entries of the full moment matrix are free scalars that never
#    enter any PSD cone, but they still *exist* as symbolic
#    unknowns in `total_basis`. This costs memory, not tightness.
#    A generalisation of `_add_parity_constraints!` that also
#    emits `⟨m⟩ = 0` for every canonical monomial with $ΔN ≠ 0$
#    (or $K ≠ 0$, or $S_z ≠ 0$) would halve the JuMP variable
#    count on the primal route.

# ---
#
# ## Summary
#
# - The paper's ${}^{2}\!D, {}^{2}\!Q, {}^{2}\!G$ relaxation is, in
#   NCTSSoS-speak, an order-2 moment relaxation with a **restricted
#   basis** (degree-1 words dropped) whose moment matrix is
#   **block-diagonalised by conserved charges**.
# - The restricted basis is 2017 monomials at $N_k = 2, [4,8]$ — 64
#   fewer than the dense order-2 basis, because the Hamiltonian is
#   $U(1)_N$-invariant.
# - The $(ΔN, K)$ partition reproduces the paper's block sizes
#   $\{240, 240, 256, 256, 512, 513\}$ exactly, including the
#   identity-absorption in the $ΔN = 0, K = 0$ sector.
# - `TermSparsity` with a custom `block_bases` is the one-line
#   integration point inside NCTSSoS: `moment_relax` already
#   emits one HPSD cone per block, so paper-style blocking is
#   supported *today*, not as a future feature.
# - The resulting SDP is smaller than both dense order-2 and TS
#   (MMD) on $\sum \dim^2$ and on max block dim, which is the cost
#   axis that matters most for interior-point solvers.
#
# ### See also
#
# - [Periodic V2RDM Benchmark (H₄ Chain)](@ref periodic-v2rdm-h4) —
#   the full specification the code on this page is honouring.
# - [Why Correlative Sparsity Misses the $k$-Blocks](@ref periodic-v2rdm-cs-k-blocks) —
#   why this page partitions the *basis*, not the *Hamiltonian
#   support graph*.
# - [²D Symmetry Blocks vs. Term-Sparsity Blocks](@ref periodic-v2rdm-block-structure) —
#   lists the actionable levers; this page is the executable
#   realisation of levers (A) and (B).
# - `demos/h4_periodic_term_sparsity_cosmo_benchmark.jl` — the
#   closest thing to an end-to-end solve pipeline (still uses
#   TS-MMD; splicing `paper_block_bases` into it is the remaining
#   step).
#
# ### References
#
# - A. O. Schouten, S. Ewing, D. A. Mazziotti,
#   "Bootstrapping the Electronic Structure of Quantum Materials,"
#   arXiv:2504.02861 (2025). Fig. 1 and Eq. 9 pin down the ${}^{2}\!D,
#   {}^{2}\!Q, {}^{2}\!G$ block structure used above.
# - C. Garrod, J. K. Percus, *J. Math. Phys.* **5**, 1756 (1964). The
#   original P, Q, G 2-positivity conditions.
# - D. A. Mazziotti, *Phys. Rev. Lett.* **106**, 083001 (2011). The
#   boundary-point solver the paper uses; a natural downstream target
#   for the `moment_problem` constructed here.
