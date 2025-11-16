# Primal/Dual Formulation Fix for Moment Matrix Extraction

## Problem Summary

The moment matrix extraction API only worked for the dual SOS formulation but failed for the primal moment formulation with a BoundsError. This was because:

- **Dual formulation**: Moment values are stored as DUAL VARIABLES of equality constraints
- **Primal formulation**: Moment values are stored as PRIMAL VARIABLES in the JuMP model

The original `extract_dual_variables()` function only extracted dual variables, causing the extraction to fail for primal formulations.

## Solution

We implemented a unified approach that detects which formulation was used and extracts values accordingly.

### Changes Made

#### 1. Modified `PolyOptResult` struct (src/interface.jl)

Added a `is_dual::Bool` field to track which formulation was used, and an additional type parameter `M2` to handle cases where the moment support uses a different monomial type than the correlative sparsity (e.g., StateWord vs NCStateWord):

```julia
struct PolyOptResult{T,P,M,M2}
    objective::T
    corr_sparsity::CorrelativeSparsity{P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{T}
    moment_support::MomentSupport{M2}
    is_dual::Bool  # NEW: true for SOS dual, false for moment primal
end
```

#### 2. Updated `cs_nctssos()` and `cs_nctssos_higher()` (src/interface.jl)

Both functions now pass the `dualize` parameter to `PolyOptResult`:

```julia
return PolyOptResult(objective_value(problem_to_solve.model),
                     corr_sparsity,
                     cliques_term_sparsities,
                     problem_to_solve.model,
                     moment_support,
                     dualize)  # Pass dualization flag
```

#### 3. Implemented `extract_primal_variables()` (src/moment_extraction_api.jl)

New function to extract primal variable values from the JuMP model:

```julia
function extract_primal_variables(model::GenericModel{T}) where {T}
    # Get all variables from the model
    all_vars = all_variables(model)

    # Extract primal values
    primal_vals = value.(all_vars)

    return primal_vals
end
```

#### 4. Updated `get_moment_matrices()` (src/moment_extraction_api.jl)

Modified to branch based on the formulation type:

```julia
function get_moment_matrices(result::PolyOptResult{T,P,M}) where {T,P,M}
    # 1. Extract values from model based on formulation type
    if result.is_dual
        # Dual SOS formulation: extract dual variables from equality constraints
        moment_values = extract_dual_variables(result.model)
    else
        # Primal moment formulation: extract primal variable values
        moment_values = extract_primal_variables(result.model)
    end

    # 2. Reconstruct moment matrices using hierarchical structure
    moment_mats = reconstruct_moment_matrices(
        moment_values,
        result.moment_support
    )

    return moment_mats
end
```

#### 5. Fixed Test Suite (test/test_moment_extraction.jl)

Updated type checking in tests - Julia's `isa` doesn't work with `Vector{Vector{<:Matrix}}` syntax. Changed to:

```julia
# Before (incorrect):
@test moments isa Vector{Vector{<:Matrix}}

# After (correct):
@test moments isa Vector
@test moments[1] isa Vector
@test moments[1][1] isa Matrix
```

## Testing Results

All tests now pass for both formulations:

```
Test Summary:            | Pass  Total
Moment Matrix Extraction |   23     23
```

Tests verified:
1. ✓ Dual SOS formulation works correctly
2. ✓ Primal moment formulation works correctly
3. ✓ Pauli algebra problems work
4. ✓ Moment matrix properties (symmetry, etc.)
5. ✓ Higher order relaxations work

## Key Design Decisions

1. **Flag-based approach**: We chose to store `is_dual::Bool` in `PolyOptResult` rather than trying to detect formulation from model structure. This is explicit, fast, and reliable.

2. **Unified extraction**: Both primal and dual values flow through the same `reconstruct_moment_matrices()` function, which only cares about having a vector of values indexed correctly.

3. **No breaking changes**: The API remains the same - `get_moment_matrices(result)` works for both formulations transparently.

## Additional Technical Details

### Type Parameter Flexibility

During implementation, we discovered that some polynomial types (like StateWord) have different representations before and after canonicalization (NCStateWord vs StateWord). To handle this:

1. **`PolyOptResult`** now has two monomial type parameters: `M` for correlative sparsity and `M2` for moment support
2. **`build_moment_support`** accepts `global_support::Vector{M2}` instead of requiring exact type match
3. **`sos_dualize`** accepts different monomial types for moment problem vs correlative sparsity

This flexibility ensures the code works correctly with all polynomial types (standard, complex, state, trace).

## Files Modified

1. `/Users/yushengzhao/projects/NCTSSoS.jl/src/interface.jl`
   - Added `is_dual` field and `M2` type parameter to `PolyOptResult`
   - Updated `cs_nctssos()` to pass `dualize` flag
   - Updated `cs_nctssos_higher()` to pass `dualize` flag

2. `/Users/yushengzhao/projects/NCTSSoS.jl/src/moment_extraction_api.jl`
   - Implemented `extract_primal_variables()` function
   - Updated `get_moment_matrices()` to handle both formulations and new type parameters
   - Updated docstrings to reflect generalized behavior

3. `/Users/yushengzhao/projects/NCTSSoS.jl/src/moment_extraction.jl`
   - Updated `build_moment_support` to accept flexible monomial type for global_support
   - Fixed basis product computation to use `expval(_neat_dot3(...))` matching moment_solver.jl

4. `/Users/yushengzhao/projects/NCTSSoS.jl/src/moment_solver.jl`
   - Modified to return tuple `(MomentProblem, MomentSupport)`
   - Builds moment support structure for primal formulation

5. `/Users/yushengzhao/projects/NCTSSoS.jl/src/sos_solver.jl`
   - Modified `sos_dualize` to accept correlative sparsity and term sparsities
   - Updated to return tuple `(SOSProblem, MomentSupport)`
   - Builds moment support structure for dual formulation
   - Added flexible type parameters to handle StateWord/NCStateWord differences

6. `/Users/yushengzhao/projects/NCTSSoS.jl/src/complex_moment_solver.jl`
   - Modified to return tuple `(ComplexMomentProblem, MomentSupport)`
   - Builds moment support structure for complex problems

7. `/Users/yushengzhao/projects/NCTSSoS.jl/test/test_moment_extraction.jl`
   - Fixed type checking assertions

## Performance Impact

Negligible - the only overhead is a single boolean check in `get_moment_matrices()`. The additional type parameters have no runtime cost.
