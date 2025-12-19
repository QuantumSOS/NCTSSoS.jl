# Research: StatePolyOpt Solver Pipeline Support

## Executive Summary

- **Goal**: Enable `cs_nctssos()` to solve `StatePolyOpt` problems (state polynomial optimization with `ς(...)` and `tr(...)` expressions)
- **Current Status**: `StatePolyOpt` type exists but solver pipeline rejects it
- **Root Cause**: `correlative_sparsity()` and `moment_relax()` only accept `Polynomial`, not `NCStatePolynomial`
- **Key Finding**: Need to add ~5 new functions and extend 2-3 existing ones
- **Estimated Complexity**: MEDIUM-HIGH (well-defined path but touches core solver)

## Research Questions Answered

### 1. What accessor functions exist for NCStatePolynomial?

**Source**: `src/FastPolynomials/src/state_polynomial.jl` (lines 565-611)

```julia
coefficients(ncsp::NCStatePolynomial) = ncsp.coeffs
monomials(ncsp::NCStatePolynomial) = ncsp.nc_state_words  # Returns Vector{NCStateWord}
terms(ncsp::NCStatePolynomial) = zip(ncsp.coeffs, ncsp.nc_state_words)
degree(ncsp::NCStatePolynomial) = maximum(degree, ncsp.nc_state_words)
variables(ncsp::NCStatePolynomial) = union of variables from all nc_state_words
```

**Key Insight**: `monomials()` returns `NCStateWord` objects, not `Monomial`. This affects how `correlative_sparsity` extracts variable indices.

### 2. How do NCStateWord work?

**Source**: `src/FastPolynomials/src/state_word.jl` (lines 380-638)

```julia
struct NCStateWord{ST<:StateType,A<:AlgebraType,T<:Integer}
    sw::StateWord{ST,A,T}     # Commutative part: <M1><M2>...<Mk>
    nc_word::Monomial{A,T}    # Non-commutative operator part: Onc
    hash::UInt64
end
```

**Key Functions**:
- `degree(ncsw)` = `degree(sw)` + `degree(nc_word)` (line 458)
- `variables(ncsw)` = union of StateWord variables and nc_word variables (lines 465-471)
- `expval(ncsw)` = converts to StateWord by adding nc_word as expectation (lines 593-595)

**Structure Details**:
- `StateWord.state_monos` = `Vector{Monomial{A,T}}` - sorted, canonicalized monomials representing expectations
- Each `Monomial` in StateWord applies involution canonicalization: `min(m, adjoint(m))`

### 3. How does correlative_sparsity work for PolyOpt?

**Source**: `src/sparse.jl` (lines 239-310)

**Key Steps**:
1. Extract all constraints: `all_cons = vcat(pop.eq_constraints, pop.ineq_constraints)`
2. Build correlative graph via `get_correlative_graph(registry, objective, all_cons)`
   - Uses `variable_indices(obj)` and `variable_indices(con)` to extract variable indices
   - For each polynomial, extracts variable indices from monomials
3. Decompose graph into cliques via `clique_decomp(G, elim_algo)`
4. Assign constraints to cliques via `assign_constraint(cliques, all_cons, registry)`
5. Generate moment matrix bases using `subregistry()` + `get_ncbasis()`

**What needs to change for StatePolyOpt**:
- Need `variable_indices(::NCStatePolynomial)` to extract indices
- Need specialized basis generation for state words

### 4. How does moment_relax build constraint matrices?

**Source**: `src/moment_solver.jl` (lines 133-194)

**Key Function**: `moment_relax(pop::PolyOpt{A,P}, corr_sparsity, term_sparsities)`

**Process**:
1. Compute total basis from all clique term sparsities (lines 139-150)
   - Uses `_neat_dot3(row_idx, mono, col_idx)` which expects `Monomial`
   - Calls `monomials(poly)` to extract constraint monomials
2. Build constraint matrices via `_build_constraint_matrix(poly, basis, cone)` (lines 72-92)
   - For each (i,j) pair, computes `sum(coef * _neat_dot3(row_idx, mono, col_idx))`
3. Select cone type based on algebra (PSD for real, HPSD for complex)

**What needs to change for StatePolyOpt**:
- `_build_constraint_matrix` needs to handle `NCStatePolynomial`
- Total basis computation needs to iterate over `NCStateWord` not `Monomial`
- `_neat_dot3` already has overloads for `NCStateWord` (see utils.jl lines 75-90)

### 5. Are there existing variable_indices implementations?

**Source**: `src/FastPolynomials/src/polynomial.jl` (lines 512-549)

```julia
function variable_indices(p::Polynomial{A,T,C}) where {A,T,C}
    result = Set{T}()
    for t in p.terms
        for idx in t.monomial.word
            push!(result, abs(idx))  # abs for fermionic/bosonic
        end
    end
    return result
end

function variable_indices(m::Monomial{A,T}) where {A<:AlgebraType, T<:Integer}
    result = Set{T}()
    for idx in m.word
        push!(result, T(abs(idx)))
    end
    return result
end
```

**Pattern to follow**: Iterate over structure, extract indices, return `Set{T}`

### 6. What functions need NCStatePolynomial equivalents?

**Functions called on polynomials in solver pipeline**:

| Function | Current Implementation | Needed for NCStatePolynomial |
|----------|----------------------|------------------------------|
| `variable_indices(p)` | Polynomial only | **NEED TO ADD** |
| `maxdegree(p)` / `degree(p)` | Polynomial only | Exists as `degree(ncsp)` |
| `monomials(p)` | Returns `Vector{Monomial}` | Returns `Vector{NCStateWord}` ✓ |
| `coefficients(p)` | Returns `Vector{C}` | Exists ✓ |
| `terms(p)` | Returns iterator | Exists ✓ |
| `symmetric_canon(m)` | Monomial only | **NEED FOR NCStateWord** |
| `expval(m)` | Monomial only | Exists for NCStateWord ✓ |
| `_neat_dot3(a, m, b)` | Monomial only | Exists for NCStateWord ✓ |

## Detailed Implementation Plan

### Phase 1: Add Missing Accessor Functions

#### 1.1 `variable_indices(::NCStatePolynomial)`

**File**: `src/FastPolynomials/src/state_polynomial.jl` (after line 611)

```julia
"""
    variable_indices(ncsp::NCStatePolynomial) -> Set{T}

Get the set of all variable indices used in the NCStatePolynomial.
"""
function variable_indices(ncsp::NCStatePolynomial{C,ST,A,T}) where {C,ST,A,T}
    result = Set{T}()
    for ncsw in ncsp.nc_state_words
        union!(result, variable_indices(ncsw))
    end
    return result
end
```

#### 1.2 `variable_indices(::NCStateWord)`

**File**: `src/FastPolynomials/src/state_word.jl` (after line 471)

Already exists! Lines 465-471:
```julia
function variables(ncsw::NCStateWord{ST,A,T}) where {ST,A,T}
    sw_vars = variables(ncsw.sw)
    for idx in ncsw.nc_word.word
        push!(sw_vars, abs(idx))
    end
    sw_vars
end
```

But need to add `variable_indices` wrapper that calls this:
```julia
variable_indices(ncsw::NCStateWord) = variables(ncsw)
```

#### 1.3 `maxdegree(::NCStatePolynomial)`

**File**: `src/FastPolynomials/src/state_polynomial.jl` (after degree function)

```julia
maxdegree(ncsp::NCStatePolynomial) = degree(ncsp)
```

### Phase 2: Add symmetric_canon for NCStateWord

**File**: `src/FastPolynomials/src/canonicalization.jl`

For state polynomial optimization, `symmetric_canon` is used to canonicalize moment matrix entries. For NCStateWord, this means:
1. Canonicalize the StateWord part (expectations are already canonicalized via involution)
2. Canonicalize the nc_word Monomial part

```julia
"""
    symmetric_canon(ncsw::NCStateWord{ST,A,T}) -> NCStateWord

Return a canonicalized NCStateWord with symmetric canonical nc_word.
"""
function symmetric_canon(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # StateWord is already canonicalized (involution invariant)
    # Just canonicalize the nc_word part
    canon_nc = symmetric_canon(ncsw.nc_word)
    NCStateWord(ncsw.sw, canon_nc)
end
```

### Phase 3: Extend correlative_sparsity for StatePolyOpt

**File**: `src/sparse.jl`

Add new method after the existing `correlative_sparsity` (around line 310):

```julia
"""
    correlative_sparsity(pop::StatePolyOpt{A,ST,P}, order::Int, elim_algo::EliminationAlgorithm)

Build correlative sparsity structure for state polynomial optimization problems.
"""
function correlative_sparsity(
    pop::StatePolyOpt{A,ST,P},
    order::Int,
    elim_algo::EliminationAlgorithm
) where {A<:AlgebraType, ST<:StateType, C<:Number, T<:Integer, P<:NCStatePolynomial{C,ST,A,T}}
    registry = pop.registry
    all_cons = vcat(pop.eq_constraints, pop.ineq_constraints)

    # Build correlative graph using variable_indices for NCStatePolynomial
    G, sorted_indices, idx_to_node = get_state_correlative_graph(registry, pop.objective, all_cons)

    # Decompose graph into cliques
    clique_node_sets = clique_decomp(G, elim_algo)

    # Convert node positions back to variable indices
    node_to_idx = Dict{Int,T}(pos => idx for (idx, pos) in idx_to_node)
    cliques = map(clique_node_sets) do node_set
        sort([node_to_idx[node] for node in node_set])
    end

    # Assign constraints to cliques
    cliques_cons, global_cons = assign_state_constraint(cliques, all_cons, registry)

    # Generate moment matrix bases using state basis
    M = NCStateWord{ST,A,T}  # Basis elements are NCStateWord, not Monomial!
    cliques_moment_matrix_bases = Vector{Vector{M}}()

    for clique_indices in cliques
        sub_reg = subregistry(registry, clique_indices)
        basis = get_state_basis(sub_reg, order; state_type=ST)
        push!(cliques_moment_matrix_bases, sorted_unique(basis))
    end

    # Compute degrees for localizing matrix basis truncation
    cliques_moment_matrix_bases_dg = [degree.(bs) for bs in cliques_moment_matrix_bases]

    # Generate localizing matrix bases
    cliques_localizing_bases = map(zip(eachindex(cliques), cliques_cons)) do (clique_idx, clique_cons_indices)
        if isempty(clique_cons_indices)
            return Vector{M}[]
        end

        cur_orders = order .- cld.(maxdegree.(all_cons[clique_cons_indices]), 2)
        cur_lengths = map(cur_orders) do o
            searchsortedfirst(cliques_moment_matrix_bases_dg[clique_idx], o + 1) - 1
        end

        map(cur_lengths) do len
            iszero(len) ? M[] : cliques_moment_matrix_bases[clique_idx][1:len]
        end
    end

    # NOTE: Return type needs StateCorrelativeSparsity or we need to generalize CorrelativeSparsity
    return StateCorrelativeSparsity{A,ST,T,P,M}(
        cliques, registry, all_cons, cliques_cons, global_cons,
        cliques_moment_matrix_bases, cliques_localizing_bases
    )
end
```

#### Helper functions needed:

```julia
"""
    get_state_correlative_graph(registry, obj::NCStatePolynomial, cons::Vector{NCStatePolynomial})

Build correlative graph for state polynomial optimization.
"""
function get_state_correlative_graph(
    registry::VariableRegistry{A,T},
    obj::NCStatePolynomial{C,ST,A,T},
    cons::Vector{<:NCStatePolynomial{C,ST,A,T}}
) where {A<:AlgebraType, T<:Integer, C<:Number, ST<:StateType}
    # Collect all unique variable indices
    all_indices = Set{T}()
    union!(all_indices, variable_indices(obj))
    for c in cons
        union!(all_indices, variable_indices(c))
    end

    sorted_indices = sort(collect(all_indices))
    nvars = length(sorted_indices)
    idx_to_node = Dict{T,Int}(idx => pos for (pos, idx) in enumerate(sorted_indices))

    G = SimpleGraph(nvars)

    # Add cliques from objective NCStateWords
    for ncsw in monomials(obj)
        positions = [idx_to_node[idx] for idx in variable_indices(ncsw)]
        add_clique!(G, positions)
    end

    # Add cliques from constraint NCStateWords
    for con in cons
        positions = [idx_to_node[idx] for idx in variable_indices(con)]
        add_clique!(G, positions)
    end

    return G, sorted_indices, idx_to_node
end

"""
    assign_state_constraint(cliques, cons::Vector{NCStatePolynomial}, registry)

Assign state polynomial constraints to cliques.
"""
function assign_state_constraint(
    cliques::Vector{Vector{T}},
    cons::Vector{<:NCStatePolynomial},
    registry::VariableRegistry{A,T}
) where {A<:AlgebraType, T<:Integer}
    clique_sets = [Set{T}(clq) for clq in cliques]

    clique_cons = map(clique_sets) do clique_set
        findall(cons) do con
            con_indices = variable_indices(con)
            issubset(con_indices, clique_set)
        end
    end

    assigned = isempty(clique_cons) ? Int[] : union(clique_cons...)
    global_cons = setdiff(1:length(cons), assigned)

    return clique_cons, collect(global_cons)
end
```

### Phase 4: New CorrelativeSparsity Type for State Polynomials

**File**: `src/sparse.jl`

The existing `CorrelativeSparsity` has:
```julia
clq_mom_mtx_bases::Vector{Vector{M}}  # where M<:Monomial
```

For state polynomials, we need bases of `NCStateWord`. Options:
1. **Parameterize** `CorrelativeSparsity` to accept different basis types
2. **Create new type** `StateCorrelativeSparsity`

**Recommendation**: Create new type (cleaner, doesn't break existing code):

```julia
"""
    StateCorrelativeSparsity{A, ST, T, P, M}

Correlative sparsity for state polynomial optimization.
Similar to CorrelativeSparsity but bases are NCStateWord instead of Monomial.
"""
struct StateCorrelativeSparsity{A<:AlgebraType, ST<:StateType, T<:Integer, P<:NCStatePolynomial, M<:NCStateWord}
    cliques::Vector{Vector{T}}
    registry::VariableRegistry{A,T}
    cons::Vector{P}
    clq_cons::Vector{Vector{Int}}
    global_cons::Vector{Int}
    clq_mom_mtx_bases::Vector{Vector{M}}
    clq_localizing_mtx_bases::Vector{Vector{Vector{M}}}
end
```

### Phase 5: Extend moment_relax for StatePolyOpt

**File**: `src/moment_solver.jl`

Need new method that works with `NCStatePolynomial` and `StateCorrelativeSparsity`:

```julia
"""
    moment_relax(pop::StatePolyOpt, corr_sparsity::StateCorrelativeSparsity, term_sparsities)

Construct moment relaxation for state polynomial optimization.
"""
function moment_relax(
    pop::StatePolyOpt{A,ST,P},
    corr_sparsity::StateCorrelativeSparsity{A,ST,TI,P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {A<:AlgebraType, ST<:StateType, TI<:Integer, C<:Number, P<:NCStatePolynomial{C,ST,A,TI}, M<:NCStateWord{ST,A,TI}}

    # Compute total basis using _neat_dot3 for NCStateWord
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        reduce(vcat, [
            # _neat_dot3 works with NCStateWord (returns NCStateWord)
            [ncsw for ncsw_m in monomials(poly) for ncsw in [_neat_dot3(rol_idx, ncsw_m, col_idx)]]
            for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities)
            for basis in term_sparsity.block_bases
            for rol_idx in basis
            for col_idx in basis
        ])
    end...)

    is_complex = _is_complex_problem(A)
    psd_cone = is_complex ? :HPSD : :PSD

    constraints = Vector{Tuple{Symbol, Matrix{P}}}()

    for (term_sparsities, cons_idx) in zip(cliques_term_sparsities, corr_sparsity.clq_cons)
        polys = [one(pop.objective); corr_sparsity.cons[cons_idx]...]

        for (term_sparsity, poly) in zip(term_sparsities, polys)
            for ts_sub_basis in term_sparsity.block_bases
                cone = poly in pop.eq_constraints ? :Zero : psd_cone
                constraint = _build_state_constraint_matrix(poly, ts_sub_basis, cone)
                push!(constraints, constraint)
            end
        end
    end

    # Global constraints
    for global_con in corr_sparsity.global_cons
        poly = corr_sparsity.cons[global_con]
        cone = poly in pop.eq_constraints ? :Zero : psd_cone
        constraint = _build_state_constraint_matrix(poly, [one(M)], cone)
        push!(constraints, constraint)
    end

    return StateMomentProblem{A, ST, TI, M, P}(pop.objective, constraints, total_basis)
end
```

#### Helper function for building constraint matrices:

```julia
"""
    _build_state_constraint_matrix(poly::NCStatePolynomial, basis::Vector{NCStateWord}, cone)

Build constraint matrix for state polynomial optimization.
"""
function _build_state_constraint_matrix(
    poly::P,
    local_basis::Vector{M},
    cone::Symbol
) where {C, ST, A, T, P<:NCStatePolynomial{C,ST,A,T}, M<:NCStateWord{ST,A,T}}
    moment_mtx = Matrix{P}(undef, length(local_basis), length(local_basis))

    for (i, row_idx) in enumerate(local_basis)
        for (j, col_idx) in enumerate(local_basis)
            # _neat_dot3 handles NCStateWord multiplication
            element_poly = sum(
                coef * NCStatePolynomial([1.0], [_neat_dot3(row_idx, ncsw, col_idx)])
                for (coef, ncsw) in zip(coefficients(poly), monomials(poly))
            )
            moment_mtx[i, j] = element_poly
        end
    end

    return (cone, moment_mtx)
end
```

### Phase 6: Term Sparsity Adaptation

The existing `term_sparsities()` and related functions in `sparse.jl` work with `Monomial` bases. They need to be extended or overloaded for `NCStateWord` bases.

**Key functions to extend**:
- `init_activated_supp` - needs NCStatePolynomial/NCStateWord versions
- `term_sparsities` - needs NCStateWord basis handling
- `get_term_sparsity_graph` - needs NCStateWord support
- `iterate_term_sparse_supp` - needs NCStateWord support

The `_neat_dot3` function already has NCStateWord overloads, which is helpful.

### Phase 7: Interface Updates

**File**: `src/interface.jl`

Update `cs_nctssos` to handle StatePolyOpt properly. The current dispatch should work because `StatePolyOpt <: OptimizationProblem{A,P}`, but internal calls need the new methods.

**Key change**: Partial objective computation (lines 99-105) uses `variable_indices(mono)` on `Monomial`, needs version for `NCStateWord`:

```julia
# Current code (Polynomial):
cliques_objective = map(corr_sparsity.cliques) do clique_indices
    clique_set = Set(clique_indices)
    reduce(+, [
        issubset(variable_indices(mono), clique_set) ? coef * mono : zero(coef) * one(mono)
        for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))
    ])
end

# Needed for NCStatePolynomial:
cliques_objective = map(corr_sparsity.cliques) do clique_indices
    clique_set = Set(clique_indices)
    reduce(+, [
        issubset(variable_indices(ncsw), clique_set) ? coef * ncsw : zero(coef) * one(ncsw)
        for (coef, ncsw) in zip(coefficients(pop.objective), monomials(pop.objective))
    ])
end
```

## Implementation Order (Recommended)

1. **Phase 1**: Add accessor functions (30 min)
   - `variable_indices(::NCStatePolynomial)`
   - `variable_indices(::NCStateWord)` alias
   - `maxdegree(::NCStatePolynomial)` alias

2. **Phase 2**: Add canonicalization (15 min)
   - `symmetric_canon(::NCStateWord)`

3. **Phase 3-4**: CorrelativeSparsity for StatePolyOpt (2 hours)
   - New `StateCorrelativeSparsity` type
   - `correlative_sparsity(::StatePolyOpt, ...)` method
   - Helper functions

4. **Phase 5**: Moment relaxation for StatePolyOpt (1.5 hours)
   - `_build_state_constraint_matrix`
   - `moment_relax(::StatePolyOpt, ...)` method

5. **Phase 6**: Term sparsity (1.5 hours)
   - Extend `init_activated_supp` for NCStateWord
   - Extend `term_sparsities` for NCStateWord
   - Extend helper functions

6. **Phase 7**: Interface updates (1 hour)
   - Update partial objective computation in `cs_nctssos`
   - Test end-to-end

## Edge Cases and Complications

### 1. Mixed State Word Types
NCStatePolynomial combines StateWord (commutative expectations) with Monomial (non-commutative operator). Need to ensure both parts contribute to sparsity graphs.

### 2. One Function for NCStateWord
Need `one(::Type{NCStateWord{ST,A,T}})` for identity elements. Already exists (state_word.jl lines 619-630).

### 3. NCStateWord Multiplication Returns NCStateWord
Unlike Polynomial multiplication which can produce multiple terms, NCStateWord multiplication returns a single NCStateWord (line 513-519). This simplifies constraint matrix construction.

### 4. Degree Computation
For StatePolyOpt order calculation, `maxdegree` on NCStatePolynomial sums StateWord and nc_word degrees. This may require adjusting how `order` is computed in `cs_nctssos`.

### 5. StatePolynomial vs NCStatePolynomial
Users create `StatePolynomial` via `ς(m)` then multiply by `one(Monomial)` to get `NCStatePolynomial`. This conversion is already implemented (state_polynomial.jl lines 373-393).

## Testing Strategy

1. **Unit tests for new accessors**:
   - `variable_indices(::NCStatePolynomial)` returns correct indices
   - `maxdegree(::NCStatePolynomial)` matches expected values

2. **Unit tests for correlative sparsity**:
   - Graph construction captures all variable correlations
   - Clique decomposition works correctly

3. **Integration tests** (from existing disabled tests):
   - State Poly Opt 7.2.0: Expected `-2.8284271321623202`
   - State Poly Opt 7.2.1: Expected `-4.0`
   - State Poly Opt 7.2.2: Expected `-5.0`

## Sources

- `src/FastPolynomials/src/state_polynomial.jl` - NCStatePolynomial type and accessors
- `src/FastPolynomials/src/state_word.jl` - NCStateWord type and operations
- `src/sparse.jl` - correlative_sparsity implementation
- `src/moment_solver.jl` - moment_relax implementation
- `src/pop.jl` - StatePolyOpt type definition
- `src/interface.jl` - cs_nctssos entry point
- `src/FastPolynomials/src/utils.jl` - _neat_dot3 for NCStateWord
- `test/state_poly_opt.jl` - Disabled tests to enable
