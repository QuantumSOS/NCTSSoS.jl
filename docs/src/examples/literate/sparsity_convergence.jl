# # [Stabilization vs. Exactness in the Sparse Hierarchy](@id sparsity-convergence)
#
# !!! note "Prerequisites"
#     This page assumes familiarity with the
#     [CS-TSSOS hierarchy](@ref cs-tssos-hierarchy), where the bound
#     ``\lambda^{\text{cs-ts}}_{d,k}`` is indexed by relaxation order ``d``
#     and term-sparsity refinement order ``k``.
#
# You have been climbing the [CS-TSSOS hierarchy](@ref cs-tssos-hierarchy),
# tightening the bound step by step — and then the graph stops changing.
# No new edges appear. The structure has *converged*. Surely the bound
# is now tight?
#
# **Not necessarily.** This is one of the most tempting misreadings
# in sparse polynomial optimization, and the distinction matters every
# time you interpret a numerical result from the hierarchy.
#
# Two entirely different things can happen when the graph freezes:
#
# | Property | What it means | What it tells you |
# |:---|:---|:---|
# | **Graph stabilization** | The term-sparsity graph no longer gains edges when ``k`` increases. | The *structure* of the relaxation has settled. |
# | **Bound exactness** | ``\lambda^{\text{cs-ts}}_{d,k}(f, S) = \lambda_{\min}(f, S)``. | The *numerical value* matches the true optimum. |
#
# The first is a structural checkpoint. The second is what you actually
# want. One does not imply the other.
#
# Wang and Magron make this precise in
# [wangExploitingTermSparsity2020](@cite) with two clean
# examples. In Example 3.4, the chordal extension chosen in the paper's
# Figure 3 stabilizes with a loose bound. That result depends on the
# extension: the package's `MMD()` completion reaches zero at the second
# step. Example 3.8 reaches the exact bound already at the first step.
# We compare these bounds below without imposing the paper's completion.

# ## Setup
#
# We use a silent Mosek backend for the SDP solves and `CairoMakie` for
# the graph figures. Example 3.4 uses the paper basis for both the
# solve and the graph, while Example 3.8 keeps the package basis for
# the numerical solve and switches to the paper basis only for the
# figure.

using NCTSSoS, Mosek, CairoMakie
import NCTSSoS

const SOLVER = MosekBackend()

# ### Reference values
#
# We pin the paper's known answers up front so every check below
# compares a live solve against a reviewed, stable target.

const PAPER_EX34_LAMBDA_MIN = 0.0
const PAPER_EX34_SPARSE = -0.0035
include(joinpath(pkgdir(NCTSSoS), "docs", "src", "examples", "literate", "data", "sparsity_convergence_helpers.jl")) #hide
nothing #hide

# ## Chordal extension matters (Example 3.4)
#
# Here is the polynomial that exposes the gap. It is unconstrained, so
# the true minimum ``\lambda_{\min}`` is determined entirely by the
# algebra. The polynomial comes from
# [wangExploitingTermSparsity2020](@cite), Example 3.4:
#
# ```math
# \begin{aligned}
# f ={}& X^2 - XY - YX + 3Y^2 - 2XYX + 2XY^2X - YZ - ZY + 6Z^2 \\
#      & + 9Y^2Z + 9ZY^2 - 54ZYZ + 142ZY^2Z.
# \end{aligned}
# ```
#
# The `rows` keyword of `Relaxation` lets us fix the basis to
# ``\{1, X, Y, Z, YX, YZ\}`` — the same six monomials used in
# Figure 3 of the paper — so our graph and theirs sit on identical
# node sets.
#
# ### Building the problem
#
# First, we set up the variables, the polynomial, and the explicit
# paper basis.

registry_34, (vars,) = site_variables(Free(), 1, 3)
X, Y, Z = vars

f_34 = X^2 - X * Y - Y * X + 3.0 * Y^2 -
    2.0 * X * Y * X + 2.0 * X * Y^2 * X -
    Y * Z - Z * Y + 6.0 * Z^2 + 9.0 * Y^2 * Z +
    9.0 * Z * Y^2 - 54.0 * Z * Y * Z +
    142.0 * Z * Y^2 * Z

problem_34 = Problem(f_34)
paper_basis_34_labels = ["1", "X", "Y", "Z", "YX", "YZ"]
paper_basis_34 = [one(X), X, Y, Z, Y*X, Y*Z]
nothing #hide

config_34_dense = Relaxation(rows=paper_basis_34)
config_34_sparse = Relaxation(rows=paper_basis_34,
    split=TermSplit(algorithm=NCTSSoS.MMD(), level=1))
config_34_higher = Relaxation(rows=paper_basis_34,
    split=TermSplit(algorithm=NCTSSoS.MMD(), level=2))
nothing #hide

const PAPER_GRAPH_34 = (
    labels = ["1", "X", "Y", "Z", "YX", "YZ"],
    coords = [(1.0, 1.0), (0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (0.0, 1.0), (2.0, 1.0)],
    initial_edges = [(1, 5), (1, 6), (2, 3), (2, 5), (3, 4), (3, 6), (4, 6)],
)
nothing #hide

# ### The dense solve: ground truth
#
# Without any sparsity exploitation, the dense relaxation on this basis
# recovers the true minimum exactly. This is the number the sparse
# hierarchy should be chasing.

built_34_dense = NCTSSoS.build(problem_34, config_34_dense);
result_34_dense = solve(built_34_dense; solver=SOLVER)
check_isapprox(
    "Example 3.4 dense bound matches λ_min",
    result_34_dense.bound,
    PAPER_EX34_LAMBDA_MIN;
    atol=1e-6,
)

# The dense bound:
result_34_dense.bound
#
# ### The first sparse step: where the gap appears
#
# Now we turn on term sparsity via `MMD()`. The solver exploits the
# block structure revealed by the term-sparsity graph — and the
# bound *drops below zero*. This is the gap: the sparse relaxation
# is strictly weaker than the dense one on this problem.

built_34_sparse = NCTSSoS.build(problem_34, config_34_sparse);
result_34_sparse = solve(built_34_sparse; solver=SOLVER)
check_isapprox(
    "Example 3.4 sparse bound matches the expected rounded sparse value",
    result_34_sparse.bound,
    PAPER_EX34_SPARSE;
    atol=1e-4,
)

# The first sparse-step bound:
result_34_sparse.bound
#
# ### Trying harder: the next sparse step
#
# `TermSplit(level=2)` propagates the support of the completed graph,
# including fill edges, as in Equation (23) of the paper. With `MMD()`,
# this enlarges the blocks and reaches zero at the second step.

built_34_higher = NCTSSoS.build(problem_34, config_34_higher);
result_34_higher = solve(built_34_higher; solver=SOLVER)
check_isapprox(
    "Example 3.4 MMD higher step reaches the exact bound",
    result_34_higher.bound,
    PAPER_EX34_LAMBDA_MIN;
    atol=1e-6,
)

# The higher sparse-step bound:
result_34_higher.bound
#
# The dense bound and the second sparse bound are both approximately
# zero; only the first sparse step stays near ``-0.0035``.

ex34_bounds = (
    dense = result_34_dense.bound,
    sparse = result_34_sparse.bound,
    sparse_higher = result_34_higher.bound,
)
nothing #hide

# The three objectives:
ex34_bounds

# ### Visualizing support propagation
#
# Completing a chordal graph again adds no edges, but that is not a
# hierarchy step. A step first extends support using the previous
# completed graph, rebuilds the graph, and then completes it.

# The seven initial edges in `PAPER_GRAPH_34` are data from Example 3.4,
# Figure 3 of [wangExploitingTermSparsity2020](@cite). Indices refer to the
# listed row words. They are not computed from package internals.
# `NCTSSoS.block_rows` gives the completed blocks: two rows are adjacent when
# one block holds both. Fill edges are completed edges minus the paper data.

labeled_edges(PAPER_GRAPH_34.initial_edges, PAPER_GRAPH_34.labels)

# Check that every paper edge belongs to a completed block.

package_edges_34_completed = completed_edges(built_34_sparse, paper_basis_34, registry_34)
missing_34 = setdiff(PAPER_GRAPH_34.initial_edges, package_edges_34_completed)
@assert isempty(missing_34) "Missing paper edges: $missing_34; blocks: $(completed_blocks(built_34_sparse, paper_basis_34, registry_34))"
package_fill_34 = sort(setdiff(package_edges_34_completed, PAPER_GRAPH_34.initial_edges))
check_equal(
    "Example 3.4 first-step block sides",
    sort(NCTSSoS.block_sides(built_34_sparse)),
    [3, 3, 3, 3],
)
check_equal(
    "Example 3.4 higher-step block sides",
    sort(NCTSSoS.block_sides(built_34_higher)),
    [4, 4, 4],
)

figure_34_graph = draw_term_graph(
    paper_basis_34_labels,
    PAPER_GRAPH_34.coords,
    PAPER_GRAPH_34.initial_edges;
    dashed_edges=package_fill_34,
    title="Example 3.4: term-sparsity graph and chordal extension",
)
nothing #hide

# Solid lines are the paper's initial edges; dashed lines are the first fill.
figure_34_graph

# The next graph uses the level-2 completed blocks. New edges appear,
# giving three blocks of side four. No solver moments are needed to build it.

package_edges_34_k2 = completed_edges(built_34_higher, paper_basis_34, registry_34)
figure_34_higher = draw_term_graph(
    paper_basis_34_labels,
    PAPER_GRAPH_34.coords,
    package_edges_34_completed;
    dashed_edges=sort(setdiff(package_edges_34_k2, package_edges_34_completed)),
    title="Example 3.4: second step with MMD",
)
figure_34_higher

# The dashed edges in this second figure are additions to the previous
# completed graph. The paper's Figure 3 chooses a different completion;
# its stabilized loose bound is not the result of this `MMD()` sequence.

# Check two further levels against the dense bound and print every completed
# block using the paper's labels. Levels 2–4 have the same completed graph.

later_34 = [NCTSSoS.build(problem_34, Relaxation(rows=paper_basis_34,
    split=TermSplit(algorithm=NCTSSoS.MMD(), level=k))) for k in 3:4];
builds_34 = [built_34_sparse, built_34_higher, later_34...];
results_34 = [result_34_sparse, result_34_higher, [solve(b; solver=SOLVER) for b in later_34]...];
[(level=k, bound=r.bound, blocks=[PAPER_GRAPH_34.labels[ids] for ids in
    completed_blocks(builds_34[k], paper_basis_34, registry_34)]) for (k,r) in enumerate(results_34)]
@assert all(isapprox(r.bound, result_34_dense.bound; atol=1e-6) for r in results_34[2:4]) #src
@assert all(completed_edges(b, paper_basis_34, registry_34) == package_edges_34_k2 for b in later_34) #src

# ## The happy path: stabilization with exactness (Example 3.8)
#
# Not every problem falls into the trap. Example 3.8 from
# [wangExploitingTermSparsity2020](@cite) is a constrained problem
# where the sparse hierarchy reaches the true optimum on the very
# first step:
#
# ```math
# f = 2 - X^2 + XY^2X - Y^2, \qquad S = \{4 - X^2 - Y^2,\; XY + YX - 2\}.
# ```
#
# The ball constraint ``4 - X^2 - Y^2 \geq 0`` keeps the variables
# bounded; the link constraint ``XY + YX = 2`` couples them. Together
# they give a well-structured feasible set where sparsity does no harm.
#
# We separate two roles below:
#
# - The **numerical solve** uses the package's own bases (dense and
#   sparse), so the bounds are package-native.
# - The **graph figure** uses the paper basis
#   ``\{1, X, Y, X^2, XY, YX, Y^2\}`` to match Figure 4.
#
# ### Building the problem

registry_38, (vars_38,) = site_variables(Free(), 1, 2)
X_38, Y_38 = vars_38

f_38 = 2.0 - X_38^2 + X_38 * Y_38^2 * X_38 - Y_38^2
g_ball_38 = 4.0 - X_38^2 - Y_38^2
g_link_38 = X_38 * Y_38 + Y_38 * X_38 - 2.0

problem_38 = Problem(f_38; constraints=(Positive(g_ball_38), Zero(g_link_38)))

config_38_dense = Relaxation(order=2)
config_38_sparse = Relaxation(order=2, split=TermSplit(algorithm=NCTSSoS.MMD(), level=1))
config_38_higher = Relaxation(order=2, split=TermSplit(algorithm=NCTSSoS.MMD(), level=2))
nothing #hide

# ### Dense bound: the target
#
# The dense relaxation at order 2 gives us the baseline. On this
# problem it already reaches the true minimum ``\lambda_{\min} = -1``.

built_38_dense = NCTSSoS.build(problem_38, config_38_dense);
result_38_dense = solve(built_38_dense; solver=SOLVER)
check_isapprox(
    "Example 3.8 dense bound matches λ_min",
    result_38_dense.bound,
    -1.0;
    atol=1e-6,
)

# The dense bound:
result_38_dense.bound

# ### First sparse step: matching immediately
#
# Unlike Example 3.4, turning on term sparsity here does *not*
# degrade the bound. The sparse relaxation reaches ``-1`` as well —
# the block structure exposed by the term-sparsity graph is rich
# enough to preserve the full strength of the SDP.

built_38_sparse = NCTSSoS.build(problem_38, config_38_sparse);
result_38_sparse = solve(built_38_sparse; solver=SOLVER)
check_isapprox(
    "Example 3.8 first sparse-step bound matches λ_min",
    result_38_sparse.bound,
    -1.0;
    atol=1e-6,
)

# The first sparse-step bound:
result_38_sparse.bound

# ### Higher sparse step: confirming stability
#
# One more refinement step confirms the picture: the bound stays
# at ``-1``, and the internal block structure may regroup but the
# numerical answer is unchanged.

built_38_higher = NCTSSoS.build(problem_38, config_38_higher);
result_38_higher = solve(built_38_higher; solver=SOLVER)
check_isapprox(
    "Example 3.8 higher sparse-step bound still matches λ_min",
    result_38_higher.bound,
    -1.0;
    atol=1e-6,
)
check_isapprox(
    "Example 3.8 higher sparse-step bound matches the first sparse step",
    result_38_higher.bound,
    result_38_sparse.bound;
    atol=1e-6,
)

# The second sparse-step bound:
result_38_higher.bound

# All three numbers land on the same value:

ex38_bounds = (
    lambda_min = -1.0,
    dense_bound = result_38_dense.bound,
    sparse_bound = result_38_sparse.bound,
    higher_bound = result_38_higher.bound,
    higher_flat = isapprox(result_38_higher.bound, result_38_sparse.bound; atol=1e-6),
)
nothing #hide

# The summary:
ex38_bounds

# ### Reconstructing Figure 4
#
# On the paper basis ``\{1, X, Y, X^2, XY, YX, Y^2\}``, the five initial
# edges below are data from Figure 4 of [wangExploitingTermSparsity2020](@cite).
# The first completed graph needs no fill. However, a full constrained
# refinement also propagates support from localizers. Level 2 adds six edges
# between the quadratic rows; levels 2–4 then agree. The bound is already
# exact at level 1, before this graph stabilizes.
#
# The earlier isolated moment-graph cycle did not include that localizer
# propagation. Here both figures come from the full public build.

const PAPER_GRAPH_38 = (
    labels = ["1", "X", "Y", "X^2", "XY", "YX", "Y^2"],
    coords = [(1.0, 2.0), (2.0, 1.5), (2.0, 0.5), (0.0, 1.0), (0.5, 0.0), (1.5, 0.0), (1.0, 1.0)],
    solid_edges = [(1, 4), (1, 5), (1, 6), (1, 7), (2, 3)],
)
nothing #hide

paper_basis_38 = [one(X_38), X_38, Y_38, X_38^2, X_38*Y_38, Y_38*X_38, Y_38^2]
nothing #hide

# The paper basis:
paper_basis_38

# Starting from these seven rows, build each level with both constraints.
# The complete localizer basis is ``\{1,X,Y\}``: largest row degree 2,
# minus half the constraint degree 2, leaves degree 1.

paper_builds_38 = [NCTSSoS.build(problem_38, Relaxation(rows=paper_basis_38,
    split=TermSplit(algorithm=NCTSSoS.MMD(), level=k))) for k in 1:4];
paper_results_38 = [solve(b; solver=SOLVER) for b in paper_builds_38];
@assert all(isapprox(r.bound, -1; atol=1e-6) for r in paper_results_38) #src
@assert all(r.sizes.localizing_rows == [3,3] for r in paper_results_38) #src
paper_edges_38_k1 = completed_edges(paper_builds_38[1], paper_basis_38, registry_38)
missing_38 = setdiff(PAPER_GRAPH_38.solid_edges, paper_edges_38_k1)
@assert isempty(missing_38) "Missing paper edges: $missing_38; blocks: $(completed_blocks(paper_builds_38[1], paper_basis_38, registry_38))"
labeled_edges(PAPER_GRAPH_38.solid_edges, PAPER_GRAPH_38.labels)

# Compare the next completed graph from the same public row inspection.

paper_edges_38_k2 = completed_edges(paper_builds_38[2], paper_basis_38, registry_38)
nothing #hide

# The first completed edge set matches the paper. The second has eleven
# edges, and the next two levels do not add any more.

check_equal(
    "Example 3.8 paper graph matches the reviewed Figure 4 edge set",
    paper_edges_38_k1,
    PAPER_GRAPH_38.solid_edges,
)
check_equal(
    "Example 3.8 level 2 has eleven completed edges",
    length(paper_edges_38_k2),
    11,
)
@assert all(completed_edges(b, paper_basis_38, registry_38) == paper_edges_38_k2 for b in paper_builds_38[3:4]) #src

[(level=k, bound=r.bound, blocks=[PAPER_GRAPH_38.labels[ids] for ids in
    completed_blocks(paper_builds_38[k], paper_basis_38, registry_38)]) for (k,r) in enumerate(paper_results_38)]

# The graph, rendered in the paper's layout:

draw_term_graph(
    PAPER_GRAPH_38.labels,
    PAPER_GRAPH_38.coords,
    PAPER_GRAPH_38.solid_edges;
    dashed_edges=setdiff(paper_edges_38_k1, PAPER_GRAPH_38.solid_edges),
    title="Example 3.8: Figure 4 reconstruction",
    markersize=46,
    fontsize=16,
)

# The second completed graph shows the six new edges as dashed lines:

draw_term_graph(PAPER_GRAPH_38.labels, PAPER_GRAPH_38.coords, paper_edges_38_k1;
    dashed_edges=setdiff(paper_edges_38_k2, paper_edges_38_k1),
    title="Example 3.8: second step with MMD")

# Five edges at level 1, no fill needed, and the bound is tight. Example 3.8 is
# the happy path: the algebraic structure of the problem allows the
# sparse relaxation to be just as strong as the dense one.

# ## What to take away
#
# Chordal completion and support extension are separate operations.
# The numerical behavior also depends on which completion is chosen:
#
# | Case | First bound | Higher bound |
# |:---|:---|:---|
# | Example 3.4, paper Figure 3 | ≈ −0.0035 | Stays loose with that completion |
# | Example 3.4, package `MMD()` | ≈ −0.0035 | ≈ 0, exact |
# | Example 3.8, package solve | ≈ −1 | ≈ −1, exact |
#
# **Stabilization means a full support-extension and completion cycle
# adds nothing**, not merely that the current graph is chordal. If a
# stabilized bound is not satisfactory, increase the relaxation order
# or choose a less sparse completion. The algebraic
# conditions under which stabilization *does* guarantee exactness are
# discussed in [wangExploitingTermSparsity2020](@cite); they must be
# verified for each problem, not assumed.
