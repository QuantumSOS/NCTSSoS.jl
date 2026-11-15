# User Interface

Load the modeling API with `using NCTSSoS.Next`. Names prefixed with `Next`
are public, qualified calls; use `import NCTSSoS.Next as Next` for that prefix.
A wrong declaration gives an `ArgumentError` whose text names the fix.

## Problem and constraints

```@docs
NCTSSoS.Next.Problem
NCTSSoS.Next.Positive
NCTSSoS.Next.Zero
NCTSSoS.Next.ExpectationPositive
NCTSSoS.Next.ExpectationZero
```

## Relaxation and solution

```@docs
NCTSSoS.Next.Relaxation
NCTSSoS.Next.solve
NCTSSoS.Next.inspect
NCTSSoS.Next.expectation
NCTSSoS.Next.gns
NCTSSoS.Next.GNSResult
NCTSSoS.Next.evaluate
```

## [Reading a result](@id next-result)

See [`solve`](@ref NCTSSoS.Next.solve) for the result fields, program sizes and backend diagnostics.

Use `expectation(result, observable)` to read represented moments. Use
`gns(result)` for checked reconstruction; a numerical optimum alone does not
prove flatness or produce an optimizer. Backend options, native task sizes and
diagnostic conventions are described in [Solver preparation](@ref next-solver-preparation).

## Symmetry

Declare symmetry on `Problem`, not on `Relaxation`. See
[Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for the pipeline and
[Clifford Symmetry Detection](@ref clifford-symmetry-detection) for detection.

```@docs
NCTSSoS.Next.SitePermutation
NCTSSoS.Next.LetterPermutation
NCTSSoS.Next.CliffordMap
NCTSSoS.Next.SpinRotation
NCTSSoS.Next.spin_rotations
NCTSSoS.Next.NumberRotation
NCTSSoS.Next.Symmetry
NCTSSoS.Next.detect_symmetry
```

## Strengthening

See [Strengthening](@ref strengthening) for validity domains and numerical
comparisons. Symmetry preserves an invariant problem; a sector restriction
changes which states are allowed.

```@docs
NCTSSoS.Next.DensityMatrix
NCTSSoS.Next.InvariantCone
NCTSSoS.Next.Stationarity
NCTSSoS.Next.StateOptimality
NCTSSoS.Next.FermionicConditions
NCTSSoS.Next.Sector
```

## Module scope

```@docs
NCTSSoS.Next
```
