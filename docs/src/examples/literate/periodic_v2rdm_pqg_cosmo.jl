# # [Periodic V2RDM: PQG Basis, Trace Constraints, and COSMO](@id periodic-v2rdm-pqg-cosmo)
#
# This page is the short, honest follow-through to the earlier periodic-V2RDM pages:
# build the paper-style PQG basis, add the same trace targets the paper fixes,
# compare the **single fat PQG block** against the **manual paper-style
# $(\Delta N, K)$ blocks**, and record what happens when we hand both routes to
# `COSMO`.
#
# There is one important correction up front.
#
# The CAR-derived PQG linear maps are **not** extra constraints in this
# formulation. They are already implicit because `${}^{2}\!D`, `${}^{2}\!Q`, and
# `${}^{2}\!G` are not modeled as independent matrices here; they are different
# slices of one shared canonical moment map. What we inject explicitly through
# `moment_eq_constraints` is only:
#
# - `N_up = 4`,
# - `N_dn = 4`,
# - `Tr(²D)`, `Tr(²Q)`, and `Tr(²G)`.
#
# The heavy lifting lives in the companion demo:
# `demos/h4_periodic_pqg_cosmo_benchmark.jl`.
# This docs page stays lightweight on purpose; the full H₄ benchmark allocates
# hundreds of GiB transiently during symbolic construction, which would be a dumb
# thing to run inside docs generation.
#
# See also:
#
# - [Periodic V2RDM Benchmark (H₄ Chain)](@ref periodic-v2rdm-h4)
# - [Periodic V2RDM: ²D Symmetry Blocks vs. Term-Sparsity Blocks](@ref periodic-v2rdm-block-structure)
# - [Periodic V2RDM: Building the Paper's PQG Relaxation](@ref periodic-v2rdm-pqg-construction)

# ## Trace targets for H₄ Nk = 2 [4,8]
#
# In the spin-orbital picture used by the demo, the periodic H₄ asset has:
#
# - `M = 32` spin-orbital modes,
# - `N = 8` electrons,
# - `M - N = 24` holes.
#
# The paper-style trace targets are therefore
#
# ```math
# \mathrm{Tr}(^{2}D) = \binom{N}{2}, \qquad
# \mathrm{Tr}(^{2}Q) = \binom{M-N}{2}, \qquad
# \mathrm{Tr}(^{2}G) = N(M-N+1).
# ```
#
# For this asset those numbers are tiny and exact.

M = 32
N = 8
n_holes = M - N

trace_targets = (
    D = N * (N - 1) ÷ 2,
    Q = n_holes * (n_holes - 1) ÷ 2,
    G = N * (n_holes + 1),
)
# `trace_targets` are the explicit scalar targets injected by the demo via
# `moment_eq_constraints`.

trace_targets

# So the H₄ benchmark fixes:
#
# - `Tr(²D) = 28`,
# - `Tr(²Q) = 276`,
# - `Tr(²G) = 200`.
#
# Those are the paper-facing equality targets. They are redundant with fixed
# particle number in the exact algebra, but still worth injecting explicitly
# because they make the benchmark intent unambiguous.

# ## Static block comparison
#
# The two routes compared by the companion demo are:
#
# 1. **single PQG block** — bypass TS completely and keep one 2017 × 2017 block;
# 2. **paper-style `(ΔN, K)` blocks** — manually partition the same PQG basis into
#    the 6 blocks the paper suggests.
#
# The paper-style block sizes for H₄ Nk = 2 are
# `[513, 512, 256, 256, 240, 240]`, where the `513` block is the `K = 0`
# particle-hole block plus the identity.

fat_block_sizes = [2017]
paper_block_sizes = [513, 512, 256, 256, 240, 240]

block_summary(block_sizes) = (
    n_psd_blocks = length(block_sizes),
    max_psd_block = maximum(block_sizes),
    sum_sq = sum(d^2 for d in block_sizes),
    sum_tri_slots = sum(d * (d + 1) ÷ 2 for d in block_sizes),
)

fat_summary = block_summary(fat_block_sizes)
paper_summary = block_summary(paper_block_sizes)

fat_summary, paper_summary

# Reading the two summaries:
#
# - the fat route keeps one PSD block of size `2017`, with `Σ dim² = 4,068,289`;
# - the paper route drops to 6 PSD blocks with max size `513`, with
#   `Σ dim² = 771,585`.
#
# That is the whole structural point. Manual sectoring cuts the PSD footprint by
# about **5.3×** on `Σ dim²` and by about **3.9×** on max block dimension.
#
# It does **not** magically make the total moment problem small — but it does stop
# the solver from eating the full 2017 × 2017 cone.

reduction_factors = (
    max_block = fat_summary.max_psd_block / paper_summary.max_psd_block,
    sum_sq = fat_summary.sum_sq / paper_summary.sum_sq,
    tri_slots = fat_summary.sum_tri_slots / paper_summary.sum_tri_slots,
)

reduction_factors

# ## Measured COSMO benchmark results
#
# The benchmark script was run on both routes with:
#
# - direct `COSMO` assemble/optimize,
# - `max_iter = 10_000`,
# - `eps_abs = eps_rel = 1e-4`,
# - `decompose = false`.
#
# The measured outcomes are summarised here as reviewed constants from the demo
# run logs in `demos/results/h4_periodic_pqg_cosmo_benchmark.md`.

benchmark_results = (
    fat = (
        shared_build_s = 208.63,
        moment_relax_s = 28.19,
        direct_vars = 1_274_786,
        direct_rows = 8_154_607,
        direct_nnz = 8_401_311,
        assemble_s = 40.79,
        reached_iter1 = false,
        note = "`COSMO.optimize!` started but printed no iteration banner before timeout.",
    ),
    paper = (
        shared_build_s = 207.98,
        moment_relax_s = 6.03,
        direct_vars = 340_546,
        direct_rows = 1_561_199,
        direct_nnz = 2_061_919,
        assemble_s = 8.05,
        reached_iter1 = false,
        note = "`COSMO.optimize!` started but printed no iteration banner before timeout.",
    ),
)

benchmark_results

# The key comparison is brutal and useful:
#
# - paper-style blocking shrinks the direct conic system from about
#   `8.15e6 × 1.27e6` to about `1.56e6 × 3.41e5`,
# - `moment_relax` drops from `28.19 s` to `6.03 s`,
# - `COSMO.assemble!` drops from `40.79 s` to `8.05 s`,
# - **but neither route reaches iteration 1 in the timed runs**.
#
# So the manual block construction is a real structural improvement, not theater.
# It just still is not enough for `COSMO` on this benchmark.

system_reduction = (
    vars = benchmark_results.fat.direct_vars / benchmark_results.paper.direct_vars,
    rows = benchmark_results.fat.direct_rows / benchmark_results.paper.direct_rows,
    nnz = benchmark_results.fat.direct_nnz / benchmark_results.paper.direct_nnz,
    moment_relax = benchmark_results.fat.moment_relax_s / benchmark_results.paper.moment_relax_s,
    assemble = benchmark_results.fat.assemble_s / benchmark_results.paper.assemble_s,
)

system_reduction

# ## How to run the companion demo
#
# The exact commands are:
#
# ```julia
# # paper-style 6-block route
# julia --project demos/h4_periodic_pqg_cosmo_benchmark.jl --route=paper --direct=true
#
# # single-fat-block comparison route
# julia --project demos/h4_periodic_pqg_cosmo_benchmark.jl --route=fat --direct=true
# ```
#
# The demo builds the Bloch Hamiltonian, injects:
#
# - `N_up = 4`,
# - `N_dn = 4`,
# - `Tr(²D) = 28`,
# - `Tr(²Q) = 276`,
# - `Tr(²G) = 200`,
#
# and then either:
#
# - keeps one fat PQG block, or
# - swaps in the manual paper-style `(ΔN, K)` blocks.
#
# It reports the symbolic moment counts, direct `COSMO` conic dimensions, and
# whether `COSMO` ever starts iterating.

run_commands = [
    "julia --project demos/h4_periodic_pqg_cosmo_benchmark.jl --route=paper --direct=true",
    "julia --project demos/h4_periodic_pqg_cosmo_benchmark.jl --route=fat --direct=true",
]

run_commands

# ## Bottom line
#
# Three conclusions survived contact with the benchmark.
#
# 1. **Explicit PQG traces are easy to add** with `moment_eq_constraints` once the
#    basis is fixed.
# 2. **Manual `(ΔN, K)` blocks beat the single fat block by a wide margin** on
#    every structural metric that matters.
# 3. **`COSMO` still stalls before iteration 1 on both routes** in the timed runs.
#
# So this work does earn something important: the package now has a concrete,
# reproducible PQG-basis benchmark path and a fair apples-to-apples comparison
# between the fat-block and paper-block formulations. What it does **not** yet earn
# is a solved H₄ periodic COSMO run.
#
# That is disappointing, but it is the useful kind of disappointing.
# At least now the disappointment has numbers.
