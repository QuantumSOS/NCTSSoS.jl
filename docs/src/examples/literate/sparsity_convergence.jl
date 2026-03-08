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
# counterexamples. Example 3.4 is the **trap**: the graph stabilizes
# immediately yet the sparse bound stays strictly loose. Example 3.8
# is the **happy path**: stabilization and exactness arrive together.
# We reproduce both below, live, so you can see the gap for yourself.

# ## Setup
#
# We use a silent Mosek backend for the SDP solves and `CairoMakie` for
# the graph figures. Example 3.4 uses the paper basis for both the
# solve and the graph, while Example 3.8 keeps the package basis for
# the numerical solve and switches to the paper basis only for the
# figure.

using NCTSSoS, MosekTools, CairoMakie

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)

# ### Reference values
#
# We pin the paper's known answers up front so every check below
# compares a live solve against a reviewed, stable target.

const PAPER_EX34_LAMBDA_MIN = 0.0
const PAPER_EX34_SPARSE = -0.0035
include(joinpath(pkgdir(NCTSSoS), "docs", "src", "examples", "literate", "data", "sparsity_convergence_helpers.jl")) #hide
nothing #hide

# ## The trap: stabilization without exactness (Example 3.4)
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
# The `moment_basis` keyword lets us fix the basis to
# ``\{1, X, Y, Z, YX, YZ\}`` — the same six monomials used in
# Figure 3 of the paper — so our graph and theirs sit on identical
# node sets.
#
# ### Building the problem
#
# First, we set up the variables, the polynomial, and the explicit
# paper basis.

registry_34, (vars,) = create_noncommutative_variables([("X", 1:3)])
X, Y, Z = vars

f_34 = X^2 - X * Y - Y * X + 3.0 * Y^2 -
    2.0 * X * Y * X + 2.0 * X * Y^2 * X -
    Y * Z - Z * Y + 6.0 * Z^2 + 9.0 * Y^2 * Z +
    9.0 * Z * Y^2 - 54.0 * Z * Y * Z +
    142.0 * Z * Y^2 * Z

pop_34 = polyopt(f_34, registry_34)
paper_basis_34_labels, paper_basis_34 = let
    NM = typeof(one(X))
    yx = only(monomials(Y * X))
    yz = only(monomials(Y * Z))
    basis_entries = [
        ("1", one(X)),
        ("X", X),
        ("Y", Y),
        ("Z", Z),
        ("YX", yx),
        ("YZ", yz),
    ]
    (first.(basis_entries), NM[last(entry) for entry in basis_entries])
end
nothing #hide

config_34_dense = SolverConfig(
    optimizer=SILENT_MOSEK,
    moment_basis=paper_basis_34,
    cs_algo=NoElimination(),
    ts_algo=NoElimination(),
)
config_34_sparse = SolverConfig(
    optimizer=SILENT_MOSEK,
    moment_basis=paper_basis_34,
    cs_algo=NoElimination(),
    ts_algo=MMD(),
)
config_34_higher = SolverConfig(
    optimizer=SILENT_MOSEK,
    cs_algo=NoElimination(),
    ts_algo=MMD(),
)
nothing #hide

const PAPER_GRAPH_34 = (
    labels = ["1", "X", "Y", "Z", "YX", "YZ"],
    coords = [(1.0, 1.0), (0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (0.0, 1.0), (2.0, 1.0)],
)
nothing #hide

# ### The dense solve: ground truth
#
# Without any sparsity exploitation, the dense relaxation on this basis
# recovers the true minimum exactly. This is the number the sparse
# hierarchy should be chasing.

result_34_dense = cs_nctssos(pop_34, config_34_dense)
check_isapprox(
    "Example 3.4 dense bound matches λ_min",
    result_34_dense.objective,
    PAPER_EX34_LAMBDA_MIN;
    atol=1e-6,
)

# The dense bound:
result_34_dense.objective
#
# ### The first sparse step: where the gap appears
#
# Now we turn on term sparsity via `MMD()`. The solver exploits the
# block structure revealed by the term-sparsity graph — and the
# bound *drops below zero*. This is the gap: the sparse relaxation
# is strictly weaker than the dense one on this problem.

result_34_sparse = cs_nctssos(pop_34, config_34_sparse)
check_isapprox(
    "Example 3.4 sparse bound matches the expected rounded sparse value",
    result_34_sparse.objective,
    PAPER_EX34_SPARSE;
    atol=1e-4,
)

# The first sparse-step bound:
result_34_sparse.objective
#
# ### Trying harder: the next sparse step
#
# Can climbing one more rung help? We call `cs_nctssos_higher` to
# refine the term-sparsity graph. The answer: no — the bound does
# not budge.

result_34_higher = cs_nctssos_higher(pop_34, result_34_sparse, config_34_higher)
check_isapprox(
    "Example 3.4 higher sparse step matches the same rounded sparse value",
    result_34_higher.objective,
    PAPER_EX34_SPARSE;
    atol=1e-4,
)

# The higher sparse-step bound:
result_34_higher.objective
#
# All three values side by side tell the whole story: dense hits ``0``,
# while both sparse steps are stuck near ``-0.0035``.

ex34_bounds = (
    dense = result_34_dense.objective,
    sparse = result_34_sparse.objective,
    sparse_higher = result_34_higher.objective,
)
nothing #hide

# The three objectives:
ex34_bounds

# The graph has stabilized — we will confirm that visually next — but the
# bound has not reached the true optimum. That is the trap.

# ### Visualizing stabilization
#
# To see *why* the hierarchy stops, look at the term-sparsity graph.
# Starting from the initial graph at ``k = 1``, a chordal completion
# adds fill edges (dashed below). A second completion attempt adds
# nothing: the graph is already chordal, so the structure is frozen.

sparsity_34 = compute_sparsity(pop_34, config_34_sparse)
package_basis_34 = only(sparsity_34.corr_sparsity.clq_mom_mtx_bases)
package_init_34 = only(sparsity_34.initial_activated_supps)

package_graph_34_k1 = NCTSSoS.get_term_sparsity_graph([one(X)], package_init_34, package_basis_34)
package_edges_34_k1 = edge_pairs(package_graph_34_k1)
completed_graph_34 = deepcopy(package_graph_34_k1)
for block in NCTSSoS.clique_decomp(package_graph_34_k1, MMD())
    NCTSSoS.add_clique!(completed_graph_34, block)
end
package_edges_34_completed = edge_pairs(completed_graph_34)
package_fill_34 = sort(setdiff(package_edges_34_completed, package_edges_34_k1))
check_equal(
    "Example 3.4 chordal extension is already stable after one completion",
    chordal_completion_edges(completed_graph_34, MMD()),
    package_edges_34_completed,
)

figure_34_graph = draw_term_graph(
    paper_basis_34_labels,
    PAPER_GRAPH_34.coords,
    package_edges_34_k1;
    dashed_edges=package_fill_34,
    title="Example 3.4: term-sparsity graph and chordal extension",
)
nothing #hide

# Solid lines are the original edges; dashed lines are the fill from
# chordal completion. A second completion adds nothing new.
figure_34_graph

# The graph is frozen, yet the objective is still ``-0.0035`` — not
# ``0``. Structural convergence happened; numerical convergence did
# not. The missing edges that *would* tighten the bound simply cannot
# be discovered by the term-sparsity refinement on this basis.

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

registry_38, (vars_38,) = create_noncommutative_variables([("X", 1:2)])
X_38, Y_38 = vars_38

f_38 = 2.0 - X_38^2 + X_38 * Y_38^2 * X_38 - Y_38^2
g_ball_38 = 4.0 - X_38^2 - Y_38^2
g_link_38 = X_38 * Y_38 + Y_38 * X_38 - 2.0

pop_38 = polyopt(
    f_38,
    registry_38;
    ineq_constraints=[g_ball_38],
    eq_constraints=[g_link_38],
)

config_38_dense = SolverConfig(
    optimizer=SILENT_MOSEK,
    order=2,
    cs_algo=NoElimination(),
    ts_algo=NoElimination(),
)
config_38_sparse = SolverConfig(
    optimizer=SILENT_MOSEK,
    order=2,
    cs_algo=NoElimination(),
    ts_algo=MMD(),
)
nothing #hide

# ### Dense bound: the target
#
# The dense relaxation at order 2 gives us the baseline. On this
# problem it already reaches the true minimum ``\lambda_{\min} = -1``.

result_38_dense = cs_nctssos(pop_38, config_38_dense)
check_isapprox(
    "Example 3.8 dense bound matches λ_min",
    result_38_dense.objective,
    -1.0;
    atol=1e-6,
)

# The dense bound:
result_38_dense.objective

# ### First sparse step: matching immediately
#
# Unlike Example 3.4, turning on term sparsity here does *not*
# degrade the bound. The sparse relaxation reaches ``-1`` as well —
# the block structure exposed by the term-sparsity graph is rich
# enough to preserve the full strength of the SDP.

result_38_sparse = cs_nctssos(pop_38, config_38_sparse)
check_isapprox(
    "Example 3.8 first sparse-step bound matches λ_min",
    result_38_sparse.objective,
    -1.0;
    atol=1e-6,
)

# The first sparse-step bound:
result_38_sparse.objective

# ### Higher sparse step: confirming stability
#
# One more refinement step confirms the picture: the bound stays
# at ``-1``, and the internal block structure may regroup but the
# numerical answer is unchanged.

result_38_higher = cs_nctssos_higher(pop_38, result_38_sparse, config_38_sparse)
check_isapprox(
    "Example 3.8 higher sparse-step bound still matches λ_min",
    result_38_higher.objective,
    -1.0;
    atol=1e-6,
)
check_isapprox(
    "Example 3.8 higher sparse-step bound matches the first sparse step",
    result_38_higher.objective,
    result_38_sparse.objective;
    atol=1e-6,
)

# The second sparse-step bound:
result_38_higher.objective

# All three numbers land on the same value:

ex38_bounds = (
    lambda_min = -1.0,
    dense_bound = result_38_dense.objective,
    sparse_bound = result_38_sparse.objective,
    higher_bound = result_38_higher.objective,
    higher_flat = isapprox(result_38_higher.objective, result_38_sparse.objective; atol=1e-6),
)
nothing #hide

# The summary:
ex38_bounds

# ### Reconstructing Figure 4
#
# On the paper basis ``\{1, X, Y, X^2, XY, YX, Y^2\}``, the
# term-sparsity graph turns out to be *already chordal* — no fill
# edges needed. That means the very first graph is its own chordal
# completion, and a refinement step produces the same graph. The
# structure is a fixed point from the start.
#
# We verify this claim in four steps:
#
# 1. Build the initial activated support from the objective and both
#    constraints.
# 2. Construct the first term-sparsity graph (``k = 1``).
# 3. Push through one refinement cycle (graph → support → graph).
# 4. Check that the edge set is unchanged.

const PAPER_GRAPH_38 = (
    labels = ["1", "X", "Y", "X^2", "XY", "YX", "Y^2"],
    coords = [(1.0, 2.0), (2.0, 1.5), (2.0, 0.5), (0.0, 1.0), (0.5, 0.0), (1.5, 0.0), (1.0, 1.0)],
    solid_edges = [(1, 4), (1, 5), (1, 6), (1, 7), (2, 3)],
)
nothing #hide

paper_basis_38 = let
    NM = typeof(one(X_38))
    x2 = only(monomials(X_38^2))
    xy = only(monomials(X_38 * Y_38))
    yx = only(monomials(Y_38 * X_38))
    y2 = only(monomials(Y_38^2))
    NM[one(X_38), X_38, Y_38, x2, xy, yx, y2]
end
nothing #hide

# The paper basis:
paper_basis_38

# Starting from these seven monomials, we build the activated support.

paper_init_38 = NCTSSoS.init_activated_supp(f_38, [g_ball_38, g_link_38], paper_basis_38)
nothing #hide

# The first term-sparsity graph:

paper_graph_38_k1 = NCTSSoS.get_term_sparsity_graph([one(X_38)], paper_init_38, paper_basis_38)
paper_edges_38_k1 = edge_pairs(paper_graph_38_k1)
paper_edges_38_k1

# One refinement cycle: graph → refined support → rebuilt graph.

paper_supp_38_k1 = NCTSSoS.term_sparsity_graph_supp(paper_graph_38_k1, paper_basis_38, one(f_38))
nothing #hide

paper_graph_38_k2 = NCTSSoS.get_term_sparsity_graph([one(X_38)], paper_supp_38_k1, paper_basis_38)
paper_edges_38_k2 = edge_pairs(paper_graph_38_k2)
nothing #hide

# Now the payoff: both edge sets match the paper, and the refinement
# changes nothing.

check_equal(
    "Example 3.8 paper graph matches the reviewed Figure 4 edge set",
    paper_edges_38_k1,
    PAPER_GRAPH_38.solid_edges,
)
check_equal(
    "Example 3.8 paper graph is unchanged after one refinement step",
    paper_edges_38_k2,
    paper_edges_38_k1,
)

# The graph, rendered in the paper's layout:

draw_term_graph(
    PAPER_GRAPH_38.labels,
    PAPER_GRAPH_38.coords,
    PAPER_GRAPH_38.solid_edges;
    title="Example 3.8: Figure 4 reconstruction",
    markersize=46,
    fontsize=16,
)

# Five edges, no fill needed, and the bound is tight. Example 3.8 is
# the happy path: the algebraic structure of the problem allows the
# sparse relaxation to be just as strong as the dense one.

# ## What to take away
#
# These two examples bracket the range of outcomes when the
# term-sparsity graph stabilizes:
#
# | | Example 3.4 | Example 3.8 |
# |:---|:---|:---|
# | **Graph stabilizes?** | Yes (one chordal completion) | Yes (already a fixed point) |
# | **Bound exact?** | No — gap of ≈ 0.0035 | Yes — dense and sparse both hit ``-1`` |
# | **Lesson** | Stabilization alone proves nothing about the value. | When the algebra cooperates, sparsity is free. |
#
# The practical rule: **graph stabilization tells you the hierarchy
# has done all it can at the current relaxation order.** If the bound
# is not yet satisfactory, the next move is to increase the relaxation
# order ``d``, not the sparse refinement order ``k``. The algebraic
# conditions under which stabilization *does* guarantee exactness are
# discussed in [wangExploitingTermSparsity2020](@cite); they must be
# verified for each problem, not assumed.
