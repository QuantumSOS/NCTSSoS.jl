# Research: Bell Inequality Tests Migration

## 1. Bell Inequality Mathematical Background

Bell inequalities are tests for nonlocal correlations in quantum mechanics. The tests optimize over measurement correlations between two parties (Alice with X measurements, Bob with Y measurements).

Key characteristics:
- **Two measurement parties**: X[1:m] (Alice), Y[1:n] (Bob)
- **Commuting groups**: X's commute with each other, Y's commute with each other
- **Projective measurements**: Operators satisfy P² = P (projectors)
- **Non-commuting across parties**: X and Y operators don't commute

## 2. Old API Analysis

```julia
# Variable creation - implicit NonCommutativeAlgebra
@ncpolyvar X[1:5] Y[1:5]

# Polynomial from string equation
p = 1.0 * eval(Meta.parse(equations[idx]))

# Optimization with special flags
pop = polyopt(p, comm_gps=[X[1:ms[idx]], Y[1:ns[idx]]], is_projective=true)

# Solving
solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=d[i])
res = cs_nctssos(pop, solver_config)
```

### Key Old Parameters
- **`comm_gps`**: Specifies commuting groups for correlative sparsity
- **`is_projective`**: Flag for projective measurement semantics

## 3. New API Mapping

### Variable Creation
```julia
# NEW: Explicit registry + typed variables
registry, (X, Y) = create_noncommutative_variables([
    ("X", 1:ms[idx]),
    ("Y", 1:ns[idx])
])
```

### polyopt() Signature
```julia
# NEW: Just objective + registry, no comm_gps/is_projective
pop = polyopt(p, registry)
```

### How Old Flags Map to New API
| Old API | New API |
|---------|---------|
| `@ncpolyvar X[1:5]` | `create_noncommutative_variables([("X", 1:5), ...])` |
| `comm_gps=[X, Y]` | Site-based index encoding in registry |
| `is_projective=true` | Not needed - algebra type handles semantics |

## 4. Algebra Type Decision

**Use `NonCommutativeAlgebra`** for Bell inequality tests:
- Preserves word order (no simplification)
- Site-based encoding handles group structure
- Correlative sparsity exploits index encoding

Alternative considered: `ProjectorAlgebra` for P² = P semantics, but the old tests didn't use algebraic simplification.

## 5. Key Code Patterns from Existing Tests

### From heisenberg.jl
```julia
registry, (x, y, z) = create_pauli_variables(1:N)
ham = sum(T(J1/4) * op[i] * op[mod1(i+1,N)] for op in [x,y,z] for i in 1:N)
pop = polyopt(ham, registry)
solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())
res = cs_nctssos(pop, solver_config)
```

### From interface.jl
```julia
registry, (x,) = create_noncommutative_variables([("x", 1:n)])
f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
pop = polyopt(f, registry; ineq_constraints=[g], eq_constraints=[h1])
```

## 6. Implementation Approach

### Challenge: Building Polynomial from String
The old code uses:
```julia
p = 1.0 * eval(Meta.parse(equations[idx]))
```

This parses strings like `"1*X[2]+1*Y[1]-1*X[1]*Y[1]..."` directly into polynomial expressions.

**Problem**: The new API creates `Monomial{NonCommutativeAlgebra, T}` objects, not Julia symbols.

**Solution**: Build polynomials programmatically instead of string parsing:
1. Parse the equation string to extract coefficients and variable indices
2. Build polynomial using the new API's variable monomials

### Alternative: Keep @ncpolyvar for backward compatibility
Check if `@ncpolyvar` still works and creates compatible types.

## 7. Files to Reference

- `src/interface.jl:62-112` - polyopt() signature
- `src/FastPolynomials/src/variable_registry.jl:266-727` - create_*_variables()
- `test/heisenberg.jl` - Pauli algebra test pattern
- `test/interface.jl` - NonCommutativeAlgebra test patterns
