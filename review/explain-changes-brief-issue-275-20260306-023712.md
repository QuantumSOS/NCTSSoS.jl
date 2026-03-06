Hi.  
Small diffs; strong guarantees.

# Change Brief

## Repo Purpose

`NCTSSoS.jl` provides tools for sparse noncommutative polynomial optimization via
structured moment-SOHS/SOS relaxations, with correlative + term sparsity support
(`README.md`).

## Diff Baseline

- Repo: `/Users/exaclior/QuantumSOS/NCTSSoS.jl`
- Base: `origin/main` at `6556c5680ca59460460b9ace951426a881574651`
- Head: `29dbbdc9553caac6086347e748d42daf03ed2ed9`
- Commits:
  - `29dbbdc test(correlated_sparsity): cover Wang-Magron example 3`
- Touched files (from `review/explain-changes-diff-issue-275-20260306-023708.md`):
  - `test/correlated_sparsity/graph_and_cliques.jl` (+58)
  - `test/data/expectations/correlated_structure.json` (+17)

## Site Index

- Site: Correlated sparsity tests
  - `test/correlated_sparsity/graph_and_cliques.jl`
- Site: Expectations fixtures
  - `test/data/expectations/correlated_structure.json`

## Site: Correlated sparsity tests

### Files

- `test/correlated_sparsity/graph_and_cliques.jl`

### Change Summary

- Add local helpers to make graph edge comparisons deterministic:
  - `_normalized_edge_pairs(graph::SimpleGraph)` sorts undirected edge pairs
    (`test/correlated_sparsity/graph_and_cliques.jl:3`).
  - `_graph_support(edges_uv, basis)` expands seed edges into activated support
    monomials via `monomials(basis[u] * basis[v])`
    (`test/correlated_sparsity/graph_and_cliques.jl:7`).
- Add a clique-decomposition regression test on a 6-cycle (“Wang-Magron Example
  3.2 concept”):
  - Compare `NoElimination()` vs `MF()` by clique-size stats and clique count
    (`test/correlated_sparsity/graph_and_cliques.jl:180`).
- Add a term-sparsity “support extension” regression test (“Wang-Magron Example
  3.1 concept”):
  - Build a small noncommutative basis and show `get_term_sparsity_graph` adds a
    specific new edge implied by activated support closure
    (`test/correlated_sparsity/graph_and_cliques.jl:196`).

### Why It Matters

- Guards two brittle-but-important behaviors:
  - Chordal extension / clique decomposition sensitivity on non-chordal graphs
    (cycle graphs).
  - Term-sparsity graph edge creation induced by activated support (ensures the
    graph grows beyond the initial “seed edges” when warranted).

### Behavior/Risk

- `MF()` decomposition outcomes can vary with upstream graph/chordal toolchain
  changes; the test intentionally checks only:
  - max clique size for `NoElimination()` vs `MF()`
  - number of cliques for `MF()`
  This reduces brittleness but can miss structural differences in the actual
  clique sets.
- The “support extension” test assumes:
  - Basis ordering is fixed (`basis = [1, x1, x2, x3, x2x3, x3x1, x1x2]`).
  - Edge ordering is normalized via `_normalized_edge_pairs`, and JSON fixtures
    match that sorted order.

### Evidence

- New tests + helpers:
  - `test/correlated_sparsity/graph_and_cliques.jl:3`
  - `test/correlated_sparsity/graph_and_cliques.jl:180`
  - `test/correlated_sparsity/graph_and_cliques.jl:196`
- Expectations lookup path for these tests:
  - JSON loaded and cases accessed by id via `correlated_structure_case(...)` in
    `test/correlated_sparsity/runtests.jl:31` and `test/correlated_sparsity/runtests.jl:34`.
- Underlying behaviors being exercised (unchanged, but directly implicated):
  - Clique decomposition entrypoint:
    `src/optimization/elimination.jl:57`
  - Term sparsity graph construction:
    `src/optimization/sparsity.jl:381`
- Verification run (local):
  - `julia --project -e 'using Pkg; Pkg.test()'`
  - Summary: `Pass 2309  Broken 2  Total 2311` (completed successfully).

### Open Questions

- The new JSON fixture note says “Provenance: issue #292 …”, while the diff file
  name references issue 275; is the provenance pointer intentional, or should
  it reference the same issue thread as the branch context?

## Site: Expectations fixtures

### Files

- `test/data/expectations/correlated_structure.json`

### Change Summary

- Add new expectations case:
  - id: `wang_magron_example_3_support_and_chordal`
  - expected:
    - `support_extension`: seed/extended/added edges
    - `minimum_chordal_extension`: clique-size stats for `NoElimination()` and
      `MF()` on a 6-cycle
  (`test/data/expectations/correlated_structure.json:182`)

### Why It Matters

- Centralizes “golden” structural expectations so tests remain readable and
  future adjustments are data-only (update fixture, not test logic).

### Behavior/Risk

- Fixture encodes *ordered* edge pairs; tests compare full vectors for equality
  after sorting observed edges, so fixture ordering must remain sorted to match.
- Fixture encodes algorithm-dependent numbers (`mf_max_clique_size`,
  `mf_n_cliques`); changes in Graphs/CliqueTrees chordal routines could require
  updating these numbers even if behavior is still “reasonable”.

### Evidence

- New case block:
  - `test/data/expectations/correlated_structure.json:182`
- Consumption of this fixture:
  - `test/correlated_sparsity/runtests.jl:31`
  - `test/correlated_sparsity/runtests.jl:34`

### Open Questions

- Should the fixture split “Example 3.1” and “Example 3.2” into separate ids, or
  is bundling them under one id preferred to keep provenance and maintenance
  together?

## Cross-Site Notes

- The new tests are wired to the new JSON fixture via the shared helper
  `correlated_structure_case(id)` in `test/correlated_sparsity/runtests.jl:34`.
- “Support extension” test is a focused probe of `get_term_sparsity_graph(...)`
  logic: edges appear when products of basis terms (with `cons_support=[1]`)
  land in `activated_supp` (`src/optimization/sparsity.jl:381`).
- “Chordal” test is a focused probe of `clique_decomp(G, algo)` behavior across
  `NoElimination()` (forced dense) and `MF()` (fill/triangulation path)
  (`src/optimization/elimination.jl:57`).

## Suggested Walkthrough Order

1. `review/explain-changes-diff-issue-275-20260306-023708.md` (what changed)
2. `test/data/expectations/correlated_structure.json:182` (new fixture contract)
3. `test/correlated_sparsity/graph_and_cliques.jl:180` (cycle/clique expectations)
4. `test/correlated_sparsity/graph_and_cliques.jl:196` (support-extension logic)
5. `src/optimization/elimination.jl:57` (why `NoElimination()` is dense)
6. `src/optimization/sparsity.jl:381` (why seed edges can induce new edges)

## Unknowns

- Whether the specific `MF()` clique counts/sizes for a 6-cycle remain stable
  across future dependency updates (currently verified by local `Pkg.test()`,
  but upstream changes could legitimately shift these stats).
- Whether “Wang-Magron Example 3.* concept” is intended as a minimal sanity
  check (as implemented) or should match a published example more literally
  (would require clarifying what exact structure is canonical).