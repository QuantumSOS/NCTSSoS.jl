# Sparsities

## Correlative Sparsity

```@docs
NCTSSoS.CorrelativeSparsity
NCTSSoS.get_correlative_graph
NCTSSoS.assign_constraint
NCTSSoS.correlative_sparsity
```

## [Term Sparsity](@id term-sparsity-graph-def)

Given a set of basis monomials ``\{b_1, \ldots, b_n\}`` indexing a
moment or localizing matrix, the **term-sparsity graph** ``G`` is the
simple graph on ``n`` vertices defined by
[wangExploitingTermSparsity2021](@cite):

```math
(i, j) \in E(G) \iff \exists\; s \in \mathcal{S} \text{ such that }
    b_i^{\dagger}\, s\, b_j \;\in\; \mathrm{supp}_{\text{act}}
    \;\text{ or }\;
    b_j^{\dagger}\, s\, b_i \;\in\; \mathrm{supp}_{\text{act}}
```

where ``\mathcal{S}`` is the support of the relevant polynomial
(objective or constraint) and ``\mathrm{supp}_{\text{act}}`` is the
**activated support** — the set of monomials that have appeared in
previous iterations (initialized from the objective, constraints, and
diagonal moment entries).

A [chordal](@ref chordal-graph-def) completion of ``G`` decomposes it into
cliques, each defining an independent block in the SDP matrix. Iterating
— computing the graph,
extracting blocks, updating the activated support — produces the
term-sparsity hierarchy at sparse orders ``k = 1, 2, \ldots``. See
[`NCTSSoS.get_term_sparsity_graph`](@ref) for the implementation and
[Stabilization vs. Exactness](@ref sparsity-convergence) for the
convergence semantics.

```@docs
NCTSSoS.TermSparsity
NCTSSoS.term_sparsities
NCTSSoS.get_term_sparsity_graph
NCTSSoS.iterate_term_sparse_supp
NCTSSoS.term_sparsity_graph_supp
```

# Eliminations

```@docs
NCTSSoS.clique_decomp
```
