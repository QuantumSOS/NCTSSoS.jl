# Contracts for backend and moment-table authors

A user of `Problem`, `Relaxation` and `solve` does not need this page.

## Solver data and results

```@docs
NCTSSoS.DataCone
NCTSSoS.MomentStandardForm
NCTSSoS.ConicSolution
```

## Backend entry points

```@docs
NCTSSoS.prepare_standard
NCTSSoS.solve_standard
```

## Moment-table access

```@docs
NCTSSoS.moment_covered
NCTSSoS.moment_terms
NCTSSoS.moment_value
NCTSSoS.refresh_coverage!
```
