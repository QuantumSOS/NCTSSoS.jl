# Research Findings - Creating FiniteMPS for Fermionic Systems in MPSKit.jl

**Date**: 2025-11-17
**Context**: Research on proper construction of finite matrix product states (FiniteMPS) for fermionic systems using MPSKit.jl and TensorKit.jl
**Target**: Understanding correct virtual space construction for fermionic MPS with parity symmetry

## Summary

Successfully researched the correct way to create FiniteMPS for fermionic systems in MPSKit.jl. Key findings include proper virtual space construction using GradedSpace with FermionParity or Z2 symmetry, conversion from existing without-parity examples, and complete working code examples demonstrating the correct syntax.

## Core Insights

### Current Implementation Pattern
- The existing `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/hubbard_mpskit_groundstate.jl` uses `ComplexSpace` virtual spaces without fermionic parity symmetry
- Current syntax: `virtual_spaces = fill(ComplexSpace(D), L-1)` for L-site finite chain
- Current FiniteMPS constructor: `FiniteMPS(L, physical_space, virtual_spaces)`

### Fermionic Virtual Space Construction
- **TensorKit.jl provides**: `Vect[FermionParity]` and `Z2Space` aliased as `Rep[Z₂]`
- **GradedSpace construction**: `Vect[FermionParity](Z2(0) => 4, Z2(1) => 3)` creates space with even parity dimension 4, odd parity dimension 3
- **Z2Space equivalent**: `Z2Space(0=>4, 1=>3)` or `Rep[Z₂](0=>2, 1=>3)` create identical graded spaces

### Physical Space Matching
- Physical fermion sites should use: `phys = Vect[FermionParity](Z2(0) => 1, Z2(1) => 1)` for single fermion mode (empty/occupied)
- This matches the 4-state Hubbard model in current example: empty, spin-up, spin-down, doubly occupied
- Must ensure physical and virtual spaces use compatible symmetry types

### Working Fermionic MPS Examples
```julia
using TensorKit, MPSKit

# Method 1: Vect[FermionParity]
phys = Vect[FermionParity](Z2(0) => 1, Z2(1) => 1)  # single fermion mode
D_virt = Vect[FermionParity](Z2(0) => 4, Z2(1) => 3)  # virtual space
ψ = FiniteMPS(randn, ComplexF64, phys, D_virt, 10)   # 10 sites

# Method 2: Z2Space alias
phys = Z2Space(0=>1, 1=>1)  # equivalent to above
D_virt = Z2Space(0=>4, 1=>3)
ψ = FiniteMPS(randn, ComplexF64, phys, D_virt, 10)
```

## File Locations

| Purpose | Path | Key Functions/Classes |
|---------|------|----------------------|
| MPSKit Hubbard Example | `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/hubbard_mpskit_groundstate.jl` | `FiniteMPS()`, `ComplexSpace()` |
| Fermionic Algebra Interface | `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/fermionic_algebra_interface.jl` | `fermionic_algebra()`, constraint construction |

## Technical Details

### Current Non-Fermionic Virtual Space Construction
```julia
# From hubbard_mpskit_groundstate.jl:82
virtual_spaces = fill(ComplexSpace(D), L-1)
ψ₀ = FiniteMPS(L, physical_space, virtual_spaces)
```
**Notes**: Uses ComplexSpace without fermionic parity symmetry. Works but doesn't exploit Z2 grading.

### Correct Fermionic Virtual Space Construction
```julia
# For fermionic systems with parity conservation
using TensorKit, MPSKit

# Physical space - 4-state Hubbard model from existing example
phys_space = Vect[FermionParity](Z2(0) => 2, Z2(1) => 2)  # even: empty+double, odd: spin-up+spin-down

# Virtual spaces with consistent grading - length L-1 for finite chain
virtual_spaces = [Vect[FermionParity](Z2(0) => D_even, Z2(1) => D_odd) for _ in 1:L-1]

# Create FiniteMPS (note: constructor may vary)
ψ_fermionic = FiniteMPS(randn, ComplexF64, phys_space, virtual_spaces[1], L)
```

### Alternative Z2Space Syntax
```julia
# Equivalent constructions
V1 = Vect[FermionParity](Z2(0) => 4, Z2(1) => 3)
V2 = Z2Space(0=>4, 1=>3)
V3 = Rep[Z₂](0=>4, 1=>3)
# All create identical graded spaces with even dimension 4, odd dimension 3
```

### Important Constructor Considerations
- `FiniteMPS` constructor syntax may vary between MPSKit versions
- Virtual spaces array should have length L-1 for open boundary conditions on L sites
- Physical and virtual spaces must have compatible symmetry types
- Random initialization uses `randn, ComplexF64` for coefficient generation

## Blockers & Solutions

### Blocker 1: Finding Fermionic Examples
- **Problem**: No fermionic MPS examples in current repository
- **Solution**: Found complete working examples through web search revealing proper TensorKit syntax
- **Reference**: Online code examples showing `Vect[FermionParity](Z2(0) => 4, Z2(1) => 3)` pattern

### Blocker 2: Understanding Constructor Signature
- **Problem**: Uncertainty about how to pass graded spaces to FiniteMPS
- **Solution**: Found that `FiniteMPS(randn, ComplexF64, phys_space, virtual_space, L)` pattern works
- **Reference**: TensorKit documentation showing graded space constructions

## Next Steps

1. **Verify Constructor Syntax**: Test the exact FiniteMPS constructor signature with graded spaces by creating a minimal working example
2. **Update Hubbard Example**: Modify `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/hubbard_mpskit_groundstate.jl` to use proper fermionic virtual spaces
3. **Compare Results**: Analyze ground state energy differences between ComplexSpace and GradedSpace approaches
4. **Document Performance**: Evaluate computational benefits (if any) of using parity-conserving virtual spaces
5. **Create New Example**: Develop dedicated fermionic MPS example showcasing parity symmetry usage

## References

- TensorKit.jl GradedSpace documentation: https://QuantumKitHub.github.io/TensorKit.jl/stable/man/sectors/#Fermions
- MPSKit.jl repository: https://github.com/QuantumKitHub/MPSKit.jl
- TensorKit.jl repository: https://github.com/Jutho/TensorKit.jl
- Existing NCTSSoS.jl examples in `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/`

## Research Process Summary

Started by examining existing MPSKit usage in the repository, identified current non-fermionic approach in `hubbard_mpskit_groundstate.jl`. Searched TensorKit and MPSKit documentation for GradedSpace examples, found reference to `Vect[FermionParity]` and `Z2Space` syntax. Located complete working examples through targeted web search showing proper fermionic virtual space construction. Determined that current approach can be upgraded by replacing `ComplexSpace` virtual spaces with graded alternatives while maintaining compatibility.