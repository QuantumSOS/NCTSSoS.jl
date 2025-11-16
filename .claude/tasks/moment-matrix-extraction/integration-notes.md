# Moment Matrix Extraction Integration - Implementation Notes

## Summary

Successfully integrated the MomentSupport structure building into the NCTSSoS solver pipeline. The moment support structure is now built during problem construction and stored in `PolyOptResult`, enabling efficient moment matrix extraction after solving.

## Changes Made

### 1. Updated `build_moment_support()` (src/moment_extraction.jl)

**Key Changes:**
- Now accepts `global_support::Vector{M}` as a parameter instead of building it internally
- This ensures exact alignment between the moment support structure and the constraint ordering
- Only processes MOMENT MATRIX blocks (first TermSparsity in each clique), not localizing matrices
- Uses `canonicalize(simplify!(mono, sa), sa)` for consistent monomial canonicalization

**Signature:**
```julia
function build_moment_support(
    corr_sparsity::CorrelativeSparsity{P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
    global_support::Vector{M},  # NEW: pre-built global support
    sa::SimplifyAlgorithm
) where {P,M}
```

**Algorithm:**
1. Create lookup dict: `support_map = Dict(mono => idx for (idx, mono) in enumerate(global_support))`
2. For each clique, iterate only over the first TermSparsity (moment matrix)
3. For each block basis and each pair (i, j):
   - Compute: `mono = neat_dot(basis[i], basis[j])`
   - Canonicalize: `mono = canonicalize(simplify!(mono, sa), sa)`
   - Look up: `dual_indices[i, j] = support_map[mono]`
4. Return `MomentSupport(cliques, global_support)`

### 2. Modified `moment_relax()` (src/moment_solver.jl)

**Key Changes:**
- Now returns a tuple: `(MomentProblem, MomentSupport{M})`
- Builds canonicalized global support after creating `monomap`
- Calls `build_moment_support()` before returning

**Implementation:**
```julia
# Build canonicalized global support for moment extraction
unsymmetrized_basis = collect(keys(monomap))
global_support = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))

# Build moment support structure
moment_support = build_moment_support(corr_sparsity, cliques_term_sparsities, global_support, sa)

return (MomentProblem(model, constraint_matrices, monomap, sa), moment_support)
```

### 3. Modified `sos_dualize()` - Real Version (src/sos_solver.jl)

**Key Changes:**
- Now accepts additional parameters: `corr_sparsity` and `cliques_term_sparsities`
- Returns a tuple: `(SOSProblem, MomentSupport{M})`
- Uses `symmetric_basis` as the global support (already canonicalized)

**Signature:**
```julia
function sos_dualize(
    moment_problem::MomentProblem{T,M},
    corr_sparsity::CorrelativeSparsity{P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {T,M,P}
```

**Implementation:**
```julia
# symmetric_basis is already built and canonicalized
symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(moment_problem.sa)))

# ... constraint creation ...

# Build moment support structure
moment_support = build_moment_support(corr_sparsity, cliques_term_sparsities, symmetric_basis, moment_problem.sa)

return (SOSProblem(dual_model), moment_support)
```

### 4. Modified `sos_dualize()` - Complex Version (src/sos_solver.jl)

**Key Changes:**
- Similar to real version
- Accepts `corr_sparsity` and `cliques_term_sparsities`
- Returns tuple with moment support
- Uses `symmetric_basis = sort(cmp.total_basis)` as global support

**Signature:**
```julia
function sos_dualize(
    cmp::ComplexMomentProblem{T,M,P},
    corr_sparsity::CorrelativeSparsity{P2,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {T,M,P,P2}
```

### 5. Modified `moment_relax()` for Complex (src/complex_moment_solver.jl)

**Key Changes:**
- Returns tuple: `(ComplexMomentProblem, MomentSupport{M})`
- Builds canonicalized global support from `total_basis`

**Implementation:**
```julia
# Build canonicalized global support for moment extraction
global_support = sorted_unique(canonicalize.(total_basis, Ref(sa)))

# Build moment support structure
moment_support = build_moment_support(corr_sparsity, cliques_term_sparsities, global_support, sa)

return (ComplexMomentProblem(simplified_objective, constraints, total_basis, sa), moment_support)
```

### 6. Updated `cs_nctssos()` (src/interface.jl)

**Key Changes:**
- Unpacks tuple from `moment_relax()`: `(moment_problem, moment_support_primal)`
- Unpacks tuple from `sos_dualize()`: `(problem_to_solve, moment_support)`
- Passes `corr_sparsity` and `cliques_term_sparsities` to `sos_dualize()`
- Passes `moment_support` to `PolyOptResult` constructor

**Implementation:**
```julia
(moment_problem, moment_support_primal) = moment_relax(pop, corr_sparsity, cliques_term_sparsities)

(problem_to_solve, moment_support) = !dualize ?
    (moment_problem, moment_support_primal) :
    sos_dualize(moment_problem, corr_sparsity, cliques_term_sparsities)

# ... solve ...

return PolyOptResult(objective_value(problem_to_solve.model), corr_sparsity, cliques_term_sparsities, problem_to_solve.model, moment_support)
```

### 7. Updated `cs_nctssos_higher()` (src/interface.jl)

**Key Changes:**
- Same pattern as `cs_nctssos()`
- Unpacks tuples from `moment_relax()` and `sos_dualize()`
- Passes moment_support to result constructor

### 8. Module Organization (src/NCTSSoS.jl)

**Key Changes:**
- Added `canonicalize` to imports from FastPolynomials
- Added `get_moment_matrices` to exports
- Split moment extraction into two files:
  - `moment_extraction.jl`: Data structures and `build_moment_support()` (included before `interface.jl`)
  - `moment_extraction_api.jl`: Public API functions that need `PolyOptResult` (included after `interface.jl`)

**Include order:**
```julia
include("moment_extraction.jl")     # Data structures + build function
include("interface.jl")              # PolyOptResult definition
include("moment_extraction_api.jl") # API functions using PolyOptResult
```

## Critical Implementation Details

### Monomial Canonicalization

The key to correctness is ensuring the same canonicalization process is used everywhere:
1. In constraint construction: monomials are canonicalized
2. In `build_moment_support()`: same canonicalization must be applied

**Correct pattern:**
```julia
mono = neat_dot(basis[i], basis[j])        # Compute rol_idx^† * col_idx
mono = canonicalize(simplify!(mono, sa), sa)  # Canonicalize
```

### Global Support Alignment

The `global_support` parameter to `build_moment_support()` is critical:
- **In `moment_relax()`**: Built from `monomap` keys after canonicalization
- **In `sos_dualize()` (real)**: Uses `symmetric_basis` (already canonicalized)
- **In `sos_dualize()` (complex)**: Uses `symmetric_basis = sort(cmp.total_basis)`

This ensures the dual variable indices in `MomentSupport` match the constraint ordering exactly.

### Moment Matrix vs Localizing Matrices

**Important:** `build_moment_support()` only processes MOMENT MATRIX blocks:
```julia
# Only process the FIRST TermSparsity (moment matrix)
moment_matrix_ts = clique_ts_vec[1]
```

The remaining TermSparsities in each clique are localizing matrices, which are not needed for moment extraction.

## Testing Verification

The module compiles successfully:
```bash
julia --project=. -e 'using NCTSSoS'
# ✓ Success
```

## Usage Example

After these changes, users can extract moment matrices as follows:

```julia
using NCTSSoS, Clarabel

# Define problem
@ncpolyvar x y z
pop = polyopt(x^2 + y^2, [x*y >= 1])

# Solve
config = SolverConfig(optimizer=Clarabel.Optimizer, order=2)
result = cs_nctssos(pop, config)

# Extract moment matrices
moments = get_moment_matrices(result)

# Access specific moment matrix
M11 = moments[1][1]  # Clique 1, Block 1
```

## Files Modified

1. `/Users/yushengzhao/projects/NCTSSoS.jl/src/moment_extraction.jl`
2. `/Users/yushengzhao/projects/NCTSSoS.jl/src/moment_solver.jl`
3. `/Users/yushengzhao/projects/NCTSSoS.jl/src/complex_moment_solver.jl`
4. `/Users/yushengzhao/projects/NCTSSoS.jl/src/sos_solver.jl`
5. `/Users/yushengzhao/projects/NCTSSoS.jl/src/interface.jl`
6. `/Users/yushengzhao/projects/NCTSSoS.jl/src/NCTSSoS.jl`

## Files Created

1. `/Users/yushengzhao/projects/NCTSSoS.jl/src/moment_extraction_api.jl`

## Next Steps

To complete the feature:
1. Write tests for the integration (test moment matrix extraction on real problems)
2. Verify correctness with known test cases
3. Add integration tests for both primal and dual formulations
4. Test with complex polynomial optimization problems
5. Performance benchmarking
