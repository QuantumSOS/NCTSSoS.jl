# MPSKit Hubbard Model Ground State Energy Example - Implementation Plan

## Overview
Based on research of MPSKit.jl and MPSKitModels.jl repositories, this plan outlines how to create a literate example that computes the ground state energy of a Fermi-Hubbard model at half-filling on a 1D chain lattice of length 4.

## Key Findings from Repository Analysis

### 1. Model Construction
The `hubbard_model` function is defined in `/Users/exaclior/projects/MPSKitModels.jl/src/models/hamiltonians.jl` (lines 321-351):

```julia
function hubbard_model(
    T::Type{<:Number} = ComplexF64,
    particle_symmetry::Type{<:Sector} = Trivial,
    spin_symmetry::Type{<:Sector} = Trivial,
    lattice::AbstractLattice = InfiniteChain(1);
    t = 1.0, U = 1.0, mu = 0.0, n::Integer = 0
)
```

**Parameters:**
- `t`: Hopping parameter (default: 1.0)
- `U`: On-site interaction parameter (default: 1.0)
- `mu`: Chemical potential (default: 0.0)
- `n`: Fixed particle number density (default: 0)

The Hamiltonian is implemented as:
```math
H = -t \sum_{\langle i,j \rangle} \sum_{\sigma} (e_{i,\sigma}^+ e_{j,\sigma}^- + c_{i,\sigma}^- c_{j,\sigma}^+) \ + U \sum_i n_{i,\uparrow}n_{i,\downarrow} - \sum_i \mu n_i
```

### 2. Ground State Computation Pattern
From the MPSKit example at `/Users/exaclior/projects/MPSKit.jl/examples/quantum1d/6.hubbard/main.jl`, the basic pattern is:

**For Infinite Systems:**
```julia
using TensorKit, MPSKit, MPSKitModels

# Set parameters
const t = 1.0
const mu = 0.0
const U = 3.0

# Create infinite chain model
H = hubbard_model(InfiniteChain(2); U, t, mu = U / 2)

# Set up tensor spaces
Vspaces = fill(Vect[FermionParity](0 => 10, 1 => 10), 2)
psi = InfiniteMPS(physicalspace(H), Vspaces)

# Find ground state
psi, = find_groundstate(psi, H; tol = 1e-3, verbosity)

# Compute expectation value
E = real(expectation_value(psi, H)) / 2
```

**For Finite Systems:**
From `/Users/exaclior/projects/MPSKit.jl/examples/quantum1d/2.haldane/main.jl`:
```julia
L = 11
chain = FiniteChain(L)
H = heisenberg_XXX(symmetry, chain; J, spin)

# Set up tensor spaces
physical_space = SU2Space(1 => 1)
virtual_space = SU2Space(0 => 12, 1 => 12, 2 => 5, 3 => 3)
ψ₀ = FiniteMPS(L, physical_space, virtual_space)

# Find ground state
ψ, envs, delta = find_groundstate(ψ₀, H, DMRG(; verbosity = 0))
E₀ = real(expectation_value(ψ, H))
```

### 3. Half-Filling Condition
From the MPSKit documentation, half-filling is achieved by:
- Setting `mu = 0`
- Using the transformed Hamiltonian: $H = - \sum_{\langle i, j \rangle, \sigma} c^{\dagger}_{i,\sigma} c_{j,\sigma} + u \sum_i (1 - 2 n_{i,\uparrow}) (1 - 2 n_{i,\downarrow})$ where $u = U/4$

### 4. Analytic Reference Solution
The example includes an analytic solution for comparison:
```julia
function hubbard_energy(u; rtol = 1.0e-12)
    integrandum(ω) = besselj0(ω) * besselj1(ω) / (1 + exp(2u * ω)) / ω
    int, err = quadgk(integrandum, 0, Inf; rtol = rtol)
    return -u - 4 * int
end
```

## Implementation Plan

### Phase 1: Basic Setup and Model Construction
1. **Import Required Packages**
   - `TensorKit`
   - `MPSKit`
   - `MPSKitModels`
   - `QuadGK` (for analytic reference)
   - `SpecialFunctions` (for Bessel functions)

2. **Define System Parameters**
   - Lattice length: L = 4
   - Hopping parameter: t = 1.0
   - On-site interaction: U = 3.0
   - Chemical potential: mu = 0.0 (for half-filling)

3. **Create Model and Lattice**
   - Use `FiniteChain(L)` for finite chain
   - Create Hubbard model with `hubbard_model(FiniteChain(4); t=t, U=U, mu=mu)`

### Phase 2: Tensor Space Configuration
1. **Determine Physical Space**
   - Use `physicalspace(H)` for the physical Hilbert space
   - This creates the fermionic Fock space $|0\rangle, |\uparrow\rangle, |\downarrow\rangle, |\uparrow\downarrow\rangle$

2. **Configure Virtual Spaces**
   - Choose appropriate bond dimension (e.g., 10-20 for small system)
   - Set up virtual spaces with `Vspace(n)` for each site

### Phase 3: Ground State Algorithm
1. **Initialize MPS**
   - Create `FiniteMPS(L, physical_space, virtual_spaces)`
   - Use random initialization for exploration

2. **Use DMRG Algorithm**
   - Apply `find_groundstate(psi, H, DMRG(; verbosity=2))`
   - Configure convergence parameters (tolerance, max iterations)

3. **Extract Energy**
   - Compute `expectation_value(psi, H)`
   - Extract real part and divide by number of sites for energy per site

### Phase 4: Verification and Documentation
1. **Analytic Comparison**
   - Implement the analytic formula for the infinite system
   - Note expected differences for finite vs infinite systems

2. **Parameter Scan (Optional)**
   - Vary U parameter to show energy dependence
   - Include a simple plot comparing numerical vs analytic results

3. **Literate Documentation**
   - Use markdown sections as shown in existing examples
   - Explain each step with mathematical background
   - Include convergence analysis

## Julia Idioms and Code Style

### Import Style for Literate Examples
```julia
# [Title](@id hubbard-ground-state)

# Description text in markdown

using TensorKit  # For tensor operations
using MPSKit     # For MPS algorithms
using MPSKitModels # For model constructions
```

### Configuration Pattern
```julia
# Configurable parameters with descriptive names
const HUBBARD_T = 1.0    # Hopping parameter
const HUBBARD_U = 3.0    # On-site interaction
const HUBBARD_MU = 0.0   # Chemical potential (half-filling)
const LATTICE_LENGTH = 4 # System size
const BOND_DIMENSION = 16 # Virtual bond dimension
```

### Math in Markdown
Use standard LaTeX syntax for equations:
```julia
md"""
The Hubbard model Hamiltonian is:
```math
H = -t \sum_{\langle i,j \rangle, \sigma} (c^\dagger_{i,\sigma} c_{j,\sigma} + h.c.) \ + U \sum_i n_{i,\uparrow}n_{i,\downarrow} - \mu \sum_{i,\sigma} n_{i,\sigma}
```
"""
```

## File Structure and Content Outline

Expected file: `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate/hubbard_mpskit_groundstate.jl`

### Sections:
1. **Introduction** - Purpose and background
2. **Mathematical Background** - Hubbard model definition
3. **Setup** - Package imports and parameters
4. **Model Construction** - Creating the Hamiltonian
5. **Tensor Space Configuration** - Setting up MPS spaces
6. **Ground State Algorithm** - DMRG computation
7. **Results** - Energy extraction and verification
8. **Conclusion** - Summary and extensions

## Testing Strategy

Verify the implementation by:
1. **Convergence Check** - Ensure energy converges with bond dimension
2. **Parameter Consistency** - Check energy scales correctly with t and U
3. **Symmetry Verification** - Confirm particle number is conserved
4. **Comparison with Exact** - For small L=4 system, compare with exact diagonalization if available

## Potential Issues and Solutions

### Issue 1: Convergence
- **Problem**: DMRG may not converge for certain parameters
- **Solution**: Adjust tolerance, increase bond dimension, use random initialization

### Issue 2: Wrong Sector
- **Problem**: May converge to wrong particle number sector
- **Solution**: Use appropriate tensor spaces, consider symmetry breaking for half-filling

### Issue 3: Performance
- **Problem**: Computation may be slow for larger bond dimensions
- **Solution**: Start with smaller bond dimensions, use term sparsity if available

## Completion Criteria

✅ **Requirements Met:**
- [ ] Uses MPSKit.jl and MPSKitModels.jl
- [ ] Computes Hubbard model ground state energy on 1D chain (L=4)
- [ ] Implements half-filling condition (mu=0)
- [ ] Follows literate example format matching existing examples
- [ ] Includes proper mathematical documentation
- [ ] Provides numerical results with verification

✅ **Code Quality:**
- [ ] Follows Julia idioms and best practices
- [ ] Clear variable names and documentation
- [ ] Proper error handling for convergence
- [ ] Appropriate comments for complex sections

✅ **Integration:**
- [ ] Fits into existing NCTSSoS.jl documentation structure
- [ ] Compatible with documentation build system
- [ ] References to external sources where appropriate