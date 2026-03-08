# [Sparsities](@id sparsities)
The main goal is to reduce the number of monomials used in indexing the moment
matrix. This is a crucial step in making the semidefinite programs more
efficient, as the size of the moment matrix directly affects the computational
cost of solving the problem. By exploiting the structure of the problem, we can
often significantly reduce the number of monomials needed, leading to
substantial performance improvements.

![Overview of how sparsity is exploited](../assets/sparsity.typ.svg)

## [Correlative Sparsity](@id correlative-sparsity)

Correlative sparsity, also known as chordal sparsity, arises from the underlying
structure of the problem's variables. In many physical systems, not all
variables are directly coupled. This lack of coupling can be represented by a
graph, where the vertices are the variables and the edges represent direct
interactions. The correlative sparsity pattern is then determined by the maximal
cliques of this graph. A clique is a subset of vertices where every two distinct
vertices are adjacent. By only considering monomials within these maximal
cliques, we can significantly reduce the size of the moment matrix. This is
because the moment matrix will be block-diagonal, with each block corresponding
to a maximal clique.

![Hand waving effect of Correlative Sparsity](../assets/sub_problem.typ.svg)

## [Term Sparsity](@id term-sparsity)

Term sparsity, also known as ideal sparsity, is a more direct way of exploiting
the structure of the polynomials involved in the problem. If a certain monomial
does not appear in any of the polynomial constraints, then it can be safely
removed from the basis of the moment matrix. This is because the corresponding
entry in the moment matrix will not be constrained by the problem, and thus can
be set to zero without affecting the solution. This type of sparsity is
particularly effective when the polynomials are sparse, i.e., they have only a
few non-zero terms.

![Hand waving effect of Term Sparsity](../assets/sub_basis.typ.svg)

### [Chordal Graphs](@id chordal-graph-def)

A graph is **chordal** if every cycle of length at least four has a
**chord**, meaning an edge joining two non-consecutive vertices of that cycle.
Chordal graphs matter here because their maximal cliques can be read off
directly and used as SDP blocks. If a term-sparsity graph is already chordal,
no extra fill-in edges are needed before clique decomposition. For a concrete
example, see Example 3.8 in
[Stabilization vs. Exactness](@ref sparsity-convergence).

### [CS-TSSOS Hierarchy](@id cs-tssos-hierarchy)

Combining both layers produces the **CS-TSSOS hierarchy**: a two-level grid of
lower bounds indexed by the relaxation order ``d`` (rows) and the term-sparsity
refinement order ``k`` (columns)
[wangExploitingTermSparsity2021](@cite):

```math
\begin{array}{ccccccc}
\lambda^{\text{cs-ts}}_{d,\,1}(f, S) & \leq & \lambda^{\text{cs-ts}}_{d,\,2}(f, S) & \leq & \cdots & \leq & \lambda^{\text{cs}}_{d}(f, S) \\[4pt]
\leq & & \leq & & & & \leq \\[4pt]
\lambda^{\text{cs-ts}}_{d{+}1,\,1}(f, S) & \leq & \lambda^{\text{cs-ts}}_{d{+}1,\,2}(f, S) & \leq & \cdots & \leq & \lambda^{\text{cs}}_{d{+}1}(f, S) \\[4pt]
\leq & & \leq & & & & \leq \\[4pt]
\vdots & & \vdots & & \ddots & & \vdots
\end{array}
```

Reading this grid:

- **Along a row** (fixed ``d``, increasing ``k``): each term-sparsity refinement
  step can only add edges to the adjacency graph, so the SDP grows weakly larger
  and the bound tightens toward the correlative-sparsity-only bound
  ``\lambda^{\text{cs}}_d``.
- **Along a column** (fixed ``k``, increasing ``d``): raising the relaxation
  order enlarges the moment matrix and tightens the bound.
- The rightmost column ``\lambda^{\text{cs}}_d`` is the bound obtained with
  correlative sparsity alone (no term-sparsity reduction).

Under certain conditions the hierarchy converges: as ``d \to \infty``, the
bounds approach ``\lambda_{\min}(f, S)``
[wangExploitingTermSparsity2020](@cite). However, at any *finite* ``(d, k)``,
graph stabilization (the term-sparsity graph becoming a fixed point in ``k``)
does not guarantee that the bound reaches ``\lambda_{\min}(f, S)``; it may even
stay below the dense SDP bound at the same ``d``.
For concrete examples illustrating this distinction, see
[Stabilization vs. Exactness](@ref sparsity-convergence).
