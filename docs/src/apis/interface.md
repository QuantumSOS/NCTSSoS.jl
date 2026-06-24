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

## Symmetry Reduction

See [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for an overview of
what is supported, the MVP scope, and reference cases including CHSH and the
Pauli charge/spatial/singlet XXX path.

```@docs
NCTSSoS.SignedPermutation
NCTSSoS.FermionicModePermutation
NCTSSoS.CliffordSymmetry
NCTSSoS.CliffordSymmetryGroup
NCTSSoS.pauli_site_permutation
NCTSSoS.pauli_contiguous_chain_basis
NCTSSoS.pauli_sign_symmetry
NCTSSoS.PauliChargeSectorSpec
NCTSSoS.PauliSingletConstraintSpec
NCTSSoS.PauliChargeBlockLabel
NCTSSoS.heisenberg_chain_hamiltonian
NCTSSoS.pauli_chain_translation
NCTSSoS.pauli_chain_reflection
NCTSSoS.pauli_global_axis_rotation_generators
NCTSSoS.heisenberg_chain_symmetry_spec
NCTSSoS.TranslationInvariantReport
NCTSSoS.TranslationInvariantResult
NCTSSoS.pauli_translation_invariant_moment_relaxation
NCTSSoS.pauli_translation_invariant_nctssos
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
