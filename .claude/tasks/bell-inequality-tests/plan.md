# Implementation Plan: Bell Inequality Tests

## Current State
- `test/bell_ineq.jl` exists but uses **non-existent** API (`@ncpolyvar`, `comm_gps`, `is_projective`)
- Tests are broken - `@ncpolyvar` is not defined in current codebase

## Approach
Rewrite `test/bell_ineq.jl` to use new FastPolynomials API while preserving test functionality.

## Steps

### 1. [x] Research complete - understand API mapping

### 2. [ ] Replace variable creation
**Old:**
```julia
@ncpolyvar X[1:5] Y[1:5]
```

**New:**
```julia
registry, (X, Y) = create_noncommutative_variables([("X", 1:5), ("Y", 1:5)])
```

### 3. [ ] Update polynomial construction
The string parsing `eval(Meta.parse(equations[idx]))` should still work because:
- `X` and `Y` will be vectors of `Monomial` objects
- `X[1]`, `X[2]`, etc. will access individual monomials
- Multiplication `*` and addition `+` work on monomials/polynomials

### 4. [ ] Replace polyopt call
**Old:**
```julia
pop = polyopt(p, comm_gps=[X[1:ms[idx]], Y[1:ns[idx]]], is_projective=true)
```

**New:**
```julia
pop = polyopt(p, registry)
```

### 5. [ ] Handle dynamic variable ranges
Original code uses `ms[idx]` and `ns[idx]` for variable ranges.
Need to create registry with max range (1:5 for both) since that's what equations use.

### 6. [ ] Test and verify
- Run tests to check objective values match expected `λd[i]`
- Skip A10 (index 58) as in original

## Testing Strategy
- Run single test instance first
- Verify objective matches expected value
- Run full test suite

## Files to Modify
- `test/bell_ineq.jl` - main changes

## Citations
- Research findings in `.claude/tasks/bell-inequality-tests/research.md`
- Pattern from `test/interface.jl` for `create_noncommutative_variables()` usage
