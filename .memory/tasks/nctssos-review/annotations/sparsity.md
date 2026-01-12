# Review: src/optimization/sparsity.jl

**Lines**: 984
**Purpose**: Correlative + Term sparsity exploitation for SDP decomposition

## Structure Overview

```
# Correlative Sparsity (lines 1-360)
CorrelativeSparsity{A,T,P,M}
├── get_correlative_graph()      # Build variable correlation graph
├── clique_decomp()              # External: CliqueTrees.jl
├── assign_constraint()          # Assign constraints to cliques
└── correlative_sparsity()       # Main entry point

# Term Sparsity (lines 360-600)
TermSparsity{M}
├── init_activated_supp()        # Initialize support
├── get_term_sparsity_graph()    # Build term graph
├── iterate_term_sparse_supp()   # Iterative refinement
└── term_sparsity_graph_supp()   # Eq 10.4 from book

# State Polynomial versions (lines 600-984)
StateCorrelativeSparsity{A,ST,T,P,M}
└── Parallel implementations for NCStatePolynomial
```

## Correlative Sparsity Algorithm (lines 263-359)

1. **Build graph** (`get_correlative_graph`): Variables connected if appear in same monomial
2. **Decompose** (`clique_decomp`): Chordal completion → maximal cliques
3. **Assign constraints**: Each constraint → clique containing its variables
4. **Generate bases**: Per-clique moment matrix bases via `get_ncbasis`
5. **Localizing bases**: Truncated by constraint degree: `order - ceil(deg/2)`

Key insight: Uses CliqueTrees.jl for chordal decomposition - mature external package.

## Term Sparsity Algorithm (lines 460-597)

Core idea: Further decompose moment matrices within each clique.

1. **Initialize support** (`init_activated_supp`, line 460):
   - Symmetric-canonicalized objective monomials
   - Constraint monomials
   - Diagonal moment entries (b†b)

2. **Build term graph** (`get_term_sparsity_graph`, line 524):
   - Nodes = basis elements
   - Edge (i,j) if basis[i]† * supp * basis[j] ∈ activated_support

3. **Iterate** (`iterate_term_sparse_supp`, line 563):
   - Clique decomposition on term graph
   - Update support via eq 10.4

## Quality Assessment

| Aspect | Rating | Notes |
|--------|--------|-------|
| Documentation | ✓✓✓ | Clear docstrings, references book equation |
| Correctness | ✓✓✓ | Matches NCTSSOS algorithm description |
| Performance | ✓✓ | O(n²) graph construction per clique |
| Code structure | ✓✓ | Duplication between Polynomial/NCStatePolynomial |

## Key Observations

### Strengths

1. **Algebra-aware simplification**: `_neat_dot3` + `simplify()` handles algebra-specific products
2. **Clean separation**: Correlative vs term sparsity clearly separated
3. **Reference implementation**: eq 10.4 citation aids verification

### Concerns

1. **Graph lookup inefficiency** (lines 539-543):
   ```julia
   if any(m in sorted_activated_supp for m in monos_lr)
   ```
   Linear scan of sorted array - could use `searchsorted` for O(log n).

2. **Index type handling** (line 115-124):
   ```julia
   abs_idx = T(abs(idx))  # Fermionic: +n = creation, -n = annihilation
   ```
   Works correctly but non-obvious. Comment explains signed convention.

3. **Duplication** between Polynomial and NCStatePolynomial versions (~200 lines).
   Could potentially abstract with a trait for "polynomial-like" types.

4. **`extract_monomials_from_basis`** (lines 183-261):
   Multiple dispatch implementations for different algebra types.
   Necessary but adds complexity. Well-documented.

5. **State version monosquare comment** (line 851-854):
   Important note about NOT including diagonal entries by default.
   Matches NCTSSOS's `monosquare=false` - crucial for correctness.

## Data Flow

```
PolyOpt
   │
   ▼
get_correlative_graph() → SimpleGraph
   │
   ▼
clique_decomp(G, algo) → Vector{Vector{Int}}  [clique indices]
   │
   ▼
assign_constraint() → (clq_cons, global_cons)
   │
   ▼
get_ncbasis() → bases per clique
   │
   ▼
CorrelativeSparsity{A,T,P,M}
   │
   ▼
init_activated_supp() → Vector{M}
   │
   ▼
get_term_sparsity_graph() → SimpleGraph (per basis)
   │
   ▼
clique_decomp(G, algo) → blocks
   │
   ▼
TermSparsity{M} per poly (moment + localizing)
```

## Code Quality: Good

Algorithm implementation correct. Some performance opportunities and code duplication, but well-structured and documented.
