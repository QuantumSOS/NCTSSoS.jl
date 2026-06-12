# User Interface

## Problem Definition

```@docs
NCTSSoS.PolyOpt
NCTSSoS.polyopt
NCTSSoS.particle_number_constraint
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

## Ground-State Physical Constraints

```@docs
NCTSSoS.PhysicalPSDConstraint
NCTSSoS.commutator_constraints
NCTSSoS.curvature_block
NCTSSoS.rdm_block
NCTSSoS.rdm_blocks
```

## Symmetry Reduction

See [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for an overview of
what is supported, the MVP scope, and the CHSH reference acceptance case.

```@docs
NCTSSoS.SignedPermutation
NCTSSoS.FermionicModePermutation
NCTSSoS.CliffordSymmetry
NCTSSoS.CliffordSymmetryGroup
NCTSSoS.SymmetrySpec
NCTSSoS.SymmetryReport
```

## SympleQ Automatic Symmetry Detection

The SympleQ pipeline recognizes Clifford symmetries of Pauli Hamiltonians
automatically from the binary symplectic tableau. See
[Pauli Symmetry Reduction](@ref pauli-clifford-symmetry) for a worked example.

```@docs
NCTSSoS.sympleq_symmetry_spec
NCTSSoS.sympleq_clifford_symmetry
NCTSSoS.SymplecticTableau
NCTSSoS.SymplecticMatrix
NCTSSoS.PhaseVector
NCTSSoS.SympleQGenerator
```

## Basis Selection

```@docs
NCTSSoS.newton_chip_basis
```
