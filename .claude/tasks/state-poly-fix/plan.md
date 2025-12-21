# State Polynomial Optimization Fix Plan (Option B - Revised)

## Status: ✅ COMPLETE

All three state polynomial optimization tests now pass:
- Test 7.2.0 (CHSH): PASS (-2.828...)
- Test 7.2.1 (squared expectations): PASS (-4.0) ← was broken
- Test 7.2.2 (covariance): PASS (-5.0)

## Problem Summary

Test 7.2.1 (squared expectations like `<xy>²`) returns ~0 instead of -4.0. The basis generation is missing key NCStateWord forms.

## Key Discovery from Main Branch

The main branch generates **all possible (StateWord, Monomial) combinations** up to degree d:

```
Degree 2 basis on main branch:
  <I>*I          (identity both parts)
  <M>*I          (single expectation, identity operator)
  <M1><M2>*I     (compound expectation, identity operator)
  <I>*M          (identity expectation, operator monomial)
  <M1>*M2        (single expectation, operator monomial) ← MIXED FORM!
```

The **mixed forms** like `<x₁>*x₂` are critical. They allow proper constraint structure:
- Test 7.2.0 needs `<I>*M` for operator constraints
- Test 7.2.1 needs `<M>*I` for squared expectations
- Mixed forms `<M1>*M2` bridge these requirements

## Root Cause (Revised)

Current `get_state_basis` only generates `<I>*M` forms. The attempted fix of adding `<M>*I` forms failed because we need the **full combinatorial structure** of all (StateWord, Monomial) pairs.

The main branch's approach:
```julia
mapreduce(vcat, 0:d) do nc_deg
    nc_basis = monomials(variables, nc_deg)      # degree nc_deg monomials
    cw_deg = d - nc_deg                           # remaining degree for StateWord
    cw_basis = [compound StateWords up to cw_deg] # includes <M>, <M1><M2>, etc.
    Iterators.product(cw_basis, nc_basis)         # ALL combinations
end
```

This generates NCStateWords where `degree(sw) + degree(nc_word) <= d`.

## Fix Strategy (Revised)

Generate **all (StateWord, Monomial) combinations** where `degree(sw) + degree(nc_word) <= d`.

This matches main branch's approach but fits into our current architecture.

## Implementation Plan

### Step 1: Rewrite `get_state_basis` to Match Main Branch Logic

Location: `src/FastPolynomials/src/state_word.jl`, lines 759-799

New approach - generate all combinations:
```julia
function get_state_basis(
    registry::VariableRegistry{A,T},
    d::Int;
    state_type::Type{ST}=Arbitrary
) where {A<:AlgebraType, T<:Integer, ST<:StateType}
    
    # Get all monomials for each degree 0..d
    all_monos_by_deg = [Monomial{A,T}[] for _ in 0:d]
    for deg in 0:d
        poly_basis = get_ncbasis(registry, deg)
        for poly in poly_basis
            for term in poly.terms
                push!(all_monos_by_deg[deg+1], term.monomial)
            end
        end
        unique!(sort!(all_monos_by_deg[deg+1]))
    end
    
    all_monos = reduce(vcat, all_monos_by_deg)
    unique!(sort!(all_monos))
    
    result = NCStateWord{ST,A,T}[]
    identity_sw = one(StateWord{ST,A,T})
    identity_mono = one(Monomial{A,T})
    
    # Generate all (StateWord, Monomial) combinations with total degree <= d
    for nc_deg in 0:d
        nc_monos = all_monos_by_deg[nc_deg+1]
        sw_max_deg = d - nc_deg
        
        # Get all StateWords up to sw_max_deg
        sw_basis = _generate_statewords_up_to_degree(all_monos, sw_max_deg, ST)
        
        # Create NCStateWords for all combinations
        for sw in sw_basis
            for nc_m in nc_monos
                if degree(sw) + degree(nc_m) <= d
                    push!(result, NCStateWord(sw, nc_m))
                end
            end
        end
    end
    
    unique!(sort!(result))
    return result
end
```

### Step 2: Implement `_generate_statewords_up_to_degree`

```julia
"""
    _generate_statewords_up_to_degree(monos, max_deg, ST) -> Vector{StateWord}

Generate all StateWords with total degree <= max_deg.
Includes:
- Identity StateWord <I> (degree 0)
- Single expectations <M> 
- Compound expectations <M1><M2>, <M1><M2><M3>, etc.

# Arguments
- `monos`: Vector of all available monomials
- `max_deg`: Maximum total degree for the StateWord
- `ST`: StateType (Arbitrary or MaxEntangled)
"""
function _generate_statewords_up_to_degree(
    monos::Vector{Monomial{A,T}},
    max_deg::Int,
    ::Type{ST}
) where {A<:AlgebraType, T<:Integer, ST<:StateType}
    
    result = StateWord{ST,A,T}[]
    
    # Always include identity StateWord
    push!(result, one(StateWord{ST,A,T}))
    
    if max_deg == 0
        return result
    end
    
    # Filter to non-identity monomials
    valid_monos = filter(m -> !isone(m) && degree(m) <= max_deg, monos)
    
    isempty(valid_monos) && return result
    
    # Use BFS to build StateWords level by level
    # Level k contains StateWords that are products of k expectations
    current_level = [(StateWord{ST}([m]), degree(m)) for m in valid_monos if degree(m) <= max_deg]
    
    # Add single expectations
    for (sw, _) in current_level
        push!(result, sw)
    end
    
    # Build compound expectations
    while !isempty(current_level)
        next_level = Tuple{StateWord{ST,A,T}, Int}[]
        seen = Set{StateWord{ST,A,T}}()
        
        for (sw, sw_deg) in current_level
            for m in valid_monos
                new_deg = sw_deg + degree(m)
                if new_deg <= max_deg
                    new_sw = StateWord{ST}(vcat(sw.state_monos, [m]))
                    if !(new_sw in seen)
                        push!(seen, new_sw)
                        push!(next_level, (new_sw, new_deg))
                        push!(result, new_sw)
                    end
                end
            end
        end
        
        current_level = next_level
    end
    
    unique!(result)
    return result
end
```

### Step 3: Remove `required_state_words` Parameter

Remove the now-unnecessary parameter from the function signature.

### Step 4: Verify Constraint Matrix Construction

The moment matrix should now include:
- Entries from `<I>*M` basis elements (for operator constraints)
- Entries from `<M>*I` basis elements (for expectation constraints)  
- Entries from `<M1>*M2` mixed basis elements (bridging both)

`_neat_dot3` already handles all these correctly - no changes needed.

### Step 5: Update Tests

Remove `@test_skip` from test 7.2.1.

## Expected Basis Structure

For degree d=2 with variables x, y:

| Form | Examples | Count |
|------|----------|-------|
| `<I>*I` | identity | 1 |
| `<I>*M` | `<I>*x`, `<I>*y`, `<I>*xy` | ~n + n² |
| `<M>*I` | `<x>*I`, `<y>*I`, `<xy>*I` | ~n + n² |
| `<M1><M2>*I` | `<x><x>*I`, `<x><y>*I` | ~n² |
| `<M>*N` | `<x>*y`, `<y>*x` | ~n² |

This matches main branch's output for degree 2.

## File Changes

| File | Change |
|------|--------|
| `src/FastPolynomials/src/state_word.jl` | Rewrite `get_state_basis`, add `_generate_statewords_up_to_degree` |
| `test/state_poly_opt.jl` | Remove `@test_skip` from test 7.2.1 |

## Testing

```bash
julia --project -e 'include("test/state_poly_opt.jl")'
```

Expected:
- Test 7.2.0: PASS (-2.828...)
- Test 7.2.1: PASS (-4.0)
- Test 7.2.2: PASS (-5.0)

## Verification Checklist

- [x] Basis includes `<I>*M` forms (existing)
- [x] Basis includes `<M>*I` forms (new)
- [x] Basis includes `<M1><M2>*I` compound forms (new)
- [x] Basis includes `<M>*N` mixed forms (new)
- [x] Total basis size matches main branch for same inputs
- [x] Test 7.2.0 passes
- [x] Test 7.2.1 passes
- [x] Test 7.2.2 passes

## Implementation Summary

### Files Changed

1. `src/FastPolynomials/src/state_word.jl`
   - Rewrote `get_state_basis` to generate all (StateWord, Monomial) combinations
   - Added `_generate_statewords_up_to_degree` helper function
   - Removed `required_state_words` parameter (no longer needed)

2. `src/sparse.jl`
   - Updated `correlative_sparsity` for StatePolyOpt to not pass `required_state_words`

3. `test/state_poly_opt.jl`
   - Changed `@test_skip` to `@test` for test 7.2.1
