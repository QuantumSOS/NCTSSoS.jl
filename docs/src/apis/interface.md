# User Interface

## Problem Definition

```@docs
NCTSSoS.PolyOpt
NCTSSoS.polyopt
```

## Variable Creation

```@docs
NCTSSoS.create_pauli_variables
NCTSSoS.create_fermionic_variables
NCTSSoS.create_bosonic_variables
NCTSSoS.create_projector_variables
NCTSSoS.create_unipotent_variables
NCTSSoS.create_noncommutative_variables
```

## Solver Interface

```@docs
NCTSSoS.SolverConfig
NCTSSoS.cs_nctssos
NCTSSoS.cs_nctssos_higher
NCTSSoS.PolyOptResult
```

## Symmetry Reduction

See [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for an overview of
what is supported, the MVP scope, and the CHSH reference acceptance case.

```@docs
NCTSSoS.SignedPermutation
NCTSSoS.SymmetrySpec
NCTSSoS.SymmetryReport
```

### Fermionic Symmetry

Fermionic-specific helpers for sector splitting (parity / particle
number / `S_z` / Abelian orbital irrep) and Casimir-based SU(2) spin
adaptation. See the worked example
[H₂/STO-3G with Fermionic Symmetry Adaptation](@ref h2-fermionic-symmetry).

```@docs
NCTSSoS.FermionicModePermutation
NCTSSoS.FermionicModeLayout
NCTSSoS.AbelianIrrepTable
NCTSSoS.FermionicSectorSpec
NCTSSoS.FermionicSectorLabel
NCTSSoS.FermionicSpinAdaptationSpec
NCTSSoS.FermionicSpinBlockLabel
```

## Basis Selection

```@docs
NCTSSoS.newton_chip_basis
```
