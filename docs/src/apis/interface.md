# User Interface

Load the modeling API with `using NCTSSoS`. Names prefixed with `NCTSSoS`
are public, qualified calls; use `import NCTSSoS` for that prefix.
A wrong declaration gives an `ArgumentError` whose text names the fix.

## Problem and constraints

```@docs
NCTSSoS.Problem
NCTSSoS.Positive
NCTSSoS.Zero
NCTSSoS.ExpectationPositive
NCTSSoS.ExpectationZero
```

`ExpectationZero` and `ExpectationPositive` each add one scalar row, not an
operator localizing matrix. They combine with symmetry, grouping, term split
and strengthening. The selected rows must cover the constrained polynomial,
and any declared symmetry must preserve the constraints.

## Relaxation and solution

```@docs
NCTSSoS.Relaxation
NCTSSoS.solve
NCTSSoS.inspect
NCTSSoS.expectation
NCTSSoS.gns
NCTSSoS.GNSResult
NCTSSoS.evaluate
```

## [Reading a result](@id next-result)

See [`solve`](@ref NCTSSoS.solve) for the result fields, program sizes and backend diagnostics.

Use `expectation(result, observable)` to read represented moments. Use
`gns(result)` for checked reconstruction; a numerical optimum alone does not
prove flatness or produce an optimizer. Backend options, native task sizes and
diagnostic conventions are described in [Solver preparation](@ref next-solver-preparation).

## Symmetry

Declare symmetry on `Problem`, not on `Relaxation`. See
[Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for the pipeline and
[Clifford Symmetry Detection](@ref clifford-symmetry-detection) for detection.

```@docs
NCTSSoS.SitePermutation
NCTSSoS.LetterPermutation
NCTSSoS.CliffordMap
NCTSSoS.SpinRotation
NCTSSoS.spin_rotations
NCTSSoS.NumberRotation
NCTSSoS.Symmetry
NCTSSoS.detect_symmetry
```

## Strengthening

See [Strengthening](@ref strengthening) for validity domains and numerical
comparisons. Symmetry preserves an invariant problem; a sector restriction
changes which states are allowed.

```@docs
NCTSSoS.DensityMatrix
NCTSSoS.InvariantCone
NCTSSoS.Stationarity
NCTSSoS.StateOptimality
NCTSSoS.FermionicConditions
NCTSSoS.Sector
```

## Module scope

```@docs
NCTSSoS
```
