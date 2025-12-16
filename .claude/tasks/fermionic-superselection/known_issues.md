# Known Issues: Fermionic Moment/SOS Hierarchy

## Issue Summary

The fermionic moment/SOS relaxation hierarchy produces incorrect bounds for ground state energy calculations. This is a fundamental issue with how fermionic algebras are integrated with the moment relaxation framework.

## Evidence

### Simple 2-Site Hopping Problem

```julia
registry, (a, a_dag) = create_fermionic_variables(1:2)
ham = a_dag[1] * a[2] + a_dag[2] * a[1]  # Expected ground state: E₀ = -1.0
```

**Results (order=2, NoElimination):**
- Primal (moment): -1.17 (valid lower bound, but should be tight at -1.0)
- Dual (SOS): -0.81 (INVALID - greater than true minimum)

### N=4 Fermionic Chain (XY Model)

```julia
# Anti-periodic boundary conditions
# Expected ground state: E₀ = -4.0
```

**Results:**
- Primal: ~-5.0 (loose bound)
- Dual: ~-1.5 with severe numerical instability

## Symptoms

1. **Primal-Dual Gap**: Strong duality should give equal primal/dual values, but they differ significantly
2. **Numerical Instability**: Solver shows `PRSTATUS = -1.00` at every iteration
3. **Oscillating Bounds**: Objective values oscillate wildly (10^9 to 10^14) during optimization
4. **Invalid Dual Bounds**: SOS dual gives bounds greater than the true minimum (violates weak duality)

## What We Ruled Out

1. **Parity constraints**: Correctly implemented (constraint injection, not basis filtering)
2. **Basis canonicalization**: Fixed in Hermitian SOS dualization (`symmetric_canon`)
3. **Sparsity algorithms**: Issue persists with `NoElimination` (no sparsity)
4. **Complex vs Real cone**: Issue persists with both HPSD and PSD cones

## Possible Root Causes

1. **Moment matrix construction**: The fermionic `_neat_dot3` using Wick's theorem may produce constraint matrices that are structurally different from what the SOS dualization expects

2. **Missing constraints**: Fermionic systems may require additional constraints beyond parity (e.g., anticommutation relations as explicit constraints)

3. **Coefficient type mismatch**: Fermionic uses `Float64` but is marked as `_is_complex_problem = true`, causing potential issues in the Hermitian embedding

4. **Localizing matrix structure**: The localizing matrices for fermionic constraints may need special handling

## Recommended Investigation

1. **Manual verification**: Construct a 2-site fermionic moment matrix by hand and compare with code output

2. **Constraint enumeration**: Print all constraints generated and verify they match physical expectations

3. **Compare with working algebra**: Trace through Pauli algebra (which works correctly) and identify differences in fermionic path

4. **Literature review**: Check how fermionic moment relaxations are formulated in academic papers (e.g., the paper cited in the original verdict)

## Files Involved

- `src/moment_solver.jl`: Moment matrix construction
- `src/sos_solver.jl`: SOS dualization (both real and Hermitian)
- `src/FastPolynomials/src/simplification/fermionic.jl`: Wick's theorem implementation
- `src/FastPolynomials/src/utils.jl`: `_neat_dot3` function

## References

- Original verdict on parity superselection (user-provided)
- Fermionic polynomial optimization literature
- NCTSSoS moment-SOS hierarchy implementation
