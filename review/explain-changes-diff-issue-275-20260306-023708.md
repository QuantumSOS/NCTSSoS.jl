# Diff Context

## Repo

- Path: /Users/exaclior/QuantumSOS/NCTSSoS.jl
- Base: origin/main (6556c5680ca59460460b9ace951426a881574651)
- Head: 29dbbdc9553caac6086347e748d42daf03ed2ed9

## Name Status

```text
M	test/correlated_sparsity/graph_and_cliques.jl
M	test/data/expectations/correlated_structure.json
```

## Diff Stat

```text
test/correlated_sparsity/graph_and_cliques.jl    | 58 ++++++++++++++++++++++++
 test/data/expectations/correlated_structure.json | 17 +++++++
 2 files changed, 75 insertions(+)
```

## Commit Log

```text
29dbbdc test(correlated_sparsity): cover Wang-Magron example 3
```

## Unified Diff

```text
diff --git a/test/correlated_sparsity/graph_and_cliques.jl b/test/correlated_sparsity/graph_and_cliques.jl
index fff47fc..b97d287 100644
--- a/test/correlated_sparsity/graph_and_cliques.jl
+++ b/test/correlated_sparsity/graph_and_cliques.jl
@@ -1,5 +1,17 @@
 # test/correlated_sparsity/graph_and_cliques.jl
 
+function _normalized_edge_pairs(graph::SimpleGraph)
+    sort([(min(src(e), dst(e)), max(src(e), dst(e))) for e in edges(graph)])
+end
+
+function _graph_support(edges_uv::Vector{Tuple{Int,Int}}, basis::Vector{NormalMonomial{A,T}}) where {A<:AlgebraType,T<:Integer}
+    support = NormalMonomial{A,T}[]
+    for (u, v) in edges_uv
+        append!(support, monomials(basis[u] * basis[v]))
+    end
+    return sort(unique(support))
+end
+
 @testset "Graph and Cliques" begin
     @testset "get_correlative_graph" begin
         @testset "Ring graph (n=4)" begin
@@ -164,6 +176,52 @@
             @test normalize_cliques(clique_decomp(graph, MaximalElimination())) ==
                 json_int_vec_vec(expected["maximal_elimination"])
         end
+
+        @testset "Wang-Magron Example 3.2 concept (n=6 cycle)" begin
+            expected = correlated_structure_case("wang_magron_example_3_support_and_chordal")
+            expected_chordal = expected["minimum_chordal_extension"]
+            graph = cycle_graph(json_int(expected_chordal["n_vertices"]))
+
+            cliques_dense = clique_decomp(graph, NoElimination())
+            cliques_mf = clique_decomp(graph, MF())
+
+            @test maximum(length.(cliques_dense)) ==
+                json_int(expected_chordal["dense_max_clique_size"])
+            @test maximum(length.(cliques_mf)) ==
+                json_int(expected_chordal["mf_max_clique_size"])
+            @test length(cliques_mf) == json_int(expected_chordal["mf_n_cliques"])
+        end
+    end
+
+    @testset "support extension" begin
+        @testset "Wang-Magron Example 3.1 concept" begin
+            expected = correlated_structure_case("wang_magron_example_3_support_and_chordal")
+            expected_support = expected["support_extension"]
+            expected_seed_edges = Tuple{Int,Int}[
+                (row[1], row[2]) for row in json_int_vec_vec(expected_support["seed_edges"])
+            ]
+            expected_extended_edges = Tuple{Int,Int}[
+                (row[1], row[2]) for row in json_int_vec_vec(expected_support["extended_edges"])
+            ]
+            expected_added_edges = Tuple{Int,Int}[
+                (row[1], row[2]) for row in json_int_vec_vec(expected_support["added_edges"])
+            ]
+
+            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
+            one_mono = one(typeof(x[1]))
+            x23 = monomials(x[2] * x[3])[1]
+            x31 = monomials(x[3] * x[1])[1]
+            x12 = monomials(x[1] * x[2])[1]
+            basis = [one_mono, x[1], x[2], x[3], x23, x31, x12]
+            activated_supp = _graph_support(expected_seed_edges, basis)
+
+            extended_graph = NCTSSoS.get_term_sparsity_graph([one_mono], activated_supp, basis)
+            observed_extended_edges = _normalized_edge_pairs(extended_graph)
+            observed_added_edges = setdiff(observed_extended_edges, expected_seed_edges)
+
+            @test observed_extended_edges == expected_extended_edges
+            @test observed_added_edges == expected_added_edges
+        end
     end
 
     @testset "assign_constraint" begin
diff --git a/test/data/expectations/correlated_structure.json b/test/data/expectations/correlated_structure.json
index 27cd104..bcaac53 100644
--- a/test/data/expectations/correlated_structure.json
+++ b/test/data/expectations/correlated_structure.json
@@ -178,6 +178,23 @@
       },
       "notes": "Provenance: issue #257 (constrained n=2 structure expectations)."
     },
+    {
+      "id": "wang_magron_example_3_support_and_chordal",
+      "expected": {
+        "support_extension": {
+          "seed_edges": [[1, 5], [3, 6]],
+          "extended_edges": [[1, 5], [3, 4], [3, 6]],
+          "added_edges": [[3, 4]]
+        },
+        "minimum_chordal_extension": {
+          "n_vertices": 6,
+          "dense_max_clique_size": 6,
+          "mf_max_clique_size": 3,
+          "mf_n_cliques": 4
+        }
+      },
+      "notes": "Provenance: issue #292 (Wang-Magron 2021 Example 3.1/3.2 concepts)."
+    },
     {
       "id": "term_sparsity_empty_basis_nonstate",
       "expected": {
```
