# Hubbard Model Exact Diagonalization Literate Example Plan

## Objective

Create a comprehensive literate example demonstrating exact diagonalization computation of the Fermi-Hubbard model ground state energy using ExactDiagonalization.jl, complementing the existing MPSKit example. This will provide a systematic comparison of exact vs. approximate methods and establish benchmark energies for validation.

## Target Audience

- Graduate students in condensed matter physics
- Researchers working with quantum many-body systems
- Users learning exact diagonalization methods in Julia
- Developers needing exact benchmarks for testing SDP relaxations

## Julia Language Idioms Summary

Based on the ExactDiagonalization.jl codebase, key Julia idioms to follow:

### Quantum Physics Domain-Specific Patterns
```julia
# Type-safe quantum number construction
quantumnumber = ℕ(length(lattice)) ⊠ 𝕊ᶻ(0)

# Constructor pattern with keyword arguments
unitcell = Lattice([0.0, 0.0]; name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])

# Generator comprehension for hilbert space
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

# Multiple dispatch on quantum statistics
Hopping(:t, -1.0, 1)  # Range parameter as final argument
```

### Performance-Oriented Patterns
```julia
# Sparse matrix representation throughout
using LinearAlgebra: eigen  # Only import needed functions

# Type stability in quantum operations
@inline ED(lattice, hilbert, terms, quantumnumbers) # Inline performance

# Memory management for large matrices
prepare!(ed)  # Build matrix
release!(ed)  # Clean memory
```

### Integration Patterns
```julia
# Package ecosystem compatibility
using QuantumLattices: Fock, Hilbert, Hopping, Hubbard
# Seamless data structure integration
```

## Example Structure and Content

### 1. Introduction and Mathematical Background

**Content**:
- Problem statement: Computing exact ground state energy of Fermi-Hubbard model
- Reference the existing MPSKit example and establish complementarity
- Mathematical formulation with proper LaTeX equations
- Physical interpretation of competing kinetic vs. interaction terms

**Code Elements**:
```julia
# Title block with proper metadata
# # [Fermi-Hubbard Model Ground State Energy with Exact Diagonalization]
# Cross-reference to MPSKit example: "Complement to tensor network approach"
```

### 2. Package Integration and Setup

**Content**:
- Demonstrate seamless integration between quantum lattice libraries
- Import strategy showing the complete tool ecosystem
- Environment preparation for reproducible results

**Implementation**:
```julia
import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))  # Match NCTSSoS pattern

# Core diagonalization
using QuantumLattices  # Fundamental quantum constructs
using ExactDiagonalization  # Our main computational framework
using LinearAlgebra: eigen  # Eigenvalue solvers

# Data analysis and visualization
using Plots  # Energy comparisons and parameter sweeps
using DataFrames  # Systematic result organization
using Statistics  # Error analysis
using Printf  # Formatted output
```

### 3. System Construction

**Content**:
- Step-by-step lattice construction with physical interpretation
- Hilbert space definition with proper quantum statistics
- Term-by-term Hamiltonian construction
- Symmetry exploitation via quantum numbers

**Implementation**: Include detailed educational comments explaining:
```julia
# Lattice construction strategy
unitcell = Lattice([0.0, 0.0]; name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])
# Explanation: Unit cell choice determines dispersion relations and symmetry

# Physical Hilbert space
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))
# Quantum numbers: Single-orbital spin-1/2 complex fermions at each site

# Conservation laws
quantumnumber = ℕ(length(lattice)) ⊠ 𝕊ᶻ(0)  # Half-filling, Sz=0
# Semantic: Particle number ⊗ Spin-z conservation for block diagonalization
```

### 4. Exact Diagonalization Workflow

**Content**:
- Proper ED algorithm construction
- Matrix preparation and memory optimization
- Eigenvalue computation with convergence criteria
- Results validation and sanity checks

**Implementation**:
```julia
# Algorithm construction with proper typing
ed = ED(lattice, hilbert, (t, U), quantumnumber)

# Explicit preparation for educational purposes
prepare!(ed)

# Specialized eigensolver with quantum sectors
eigensystem = eigen(ed; nev=1)  # Ground state only

# Results extraction with error handling
E0_exact = eigensystem.values[1]
ψ0_exact = eigensystem.vectors[1]
sector = eigensystem.sectors[1]
```

### 5. Systematic Parameter Study

**Content**:
- Comprehensive U/t parameter sweep showing physically relevant regimes
- Systematic comparison with analytical and MPS results
- Convergence analysis demonstrating finite-size scaling
- Energy decomposition (kinetic vs. interaction contributions)

**Implementation**:
```julia
# Define parameter grid spanning physical regimes
U_values = range(0.0, 10.0, length=21)  # 0 to 10t covering weak/intermediate/strong coupling

# Results dataframe for systematic analysis
results = DataFrame(
    U=@view U_values[:],
    exact_ground_state=Float64[],
    kinetic_energy=Float64[],
    interaction_energy=Float64[],
    total_particles=Int[],
    spin_z=Float64[]
)

# Efficient computation loop with conservation checks
for U_u in U_values
    # Build new Hamiltonian with updated U
    U_term = Hubbard(:U, U_u)
    ed_updated = update!(deepcopy(ed); U=U_u)

    # Solve and extract observables
    eigensystem = eigen(ed_updated; nev=1)

    # Extract ground state properties
    E0 = eigensystem.values[1]
    ψ0 = eigensystem.vectors[1]

    # Compute separate energy contributions
    E_kinetic = expectation_value(ψ0, t_term_matrix)
    E_interaction = expectation_value(ψ0, U_u * density_density_matrix)
end
```

### 6. Results Analysis and Visualization

**Content**:
- Ground state energy vs. U/t with proper physical interpretation
- Comparison plots showing exact vs. MPS vs. SDP results
- Energy decomposition showing competing kinetic/interaction effects
- Quality assessment comparing different computational approaches

**Implementation**:
```julia
# Publication-quality plotting consistent with NCTSSoS style
p1 = plot(U_values, E_results,
         label="Exact Diagonalization", linewidth=2,
         xlabel="U/t (interaction/hopping)", ylabel="E₀",
         title="Ground State Energy vs. Intensity")

# Add comparison with MPS results from existing example
p1 = plot!(p1, U_values, mps_results,
           label="Matrix Product States (D=16)", linewidth=2, linestyle=:dash,
           alpha=0.8)

# Include analytical reference when available
p1 = vline!(p1, [4.0], label="U=4t (benchmark)", linestyle=:dot)
```

### 7. Critical Analysis and Limitations

**Content**:
- Systematic comparison of computational complexity (ED vs. MPS vs. SDP)
- Memory requirements and scaling behavior analysis
- Accuracy validation against known analytical results
- Physical interpretation of obtained energies

**Implementation**:
```julia
# Computational complexity analysis
trues <- "ED complexity: O(dim³) where dim = C(4L, 2L) at half-filling"
trues <- "Physical memory scales as O(dim²) for sparse matrices"
trues <- "Quantum number conservation reduces dimension by ~4× for Sz=0"
```

### 8. Benchmark Comparison and Integration

**Content**:
- Direct comparison with the existing MPSKit example
- Cross-validation showing different computational approaches
- Integration with NCTSSoS testing framework
- Provision of exact reference values for future validation

**Implementation**:
```julia
# Reference values for NCTSSoS testing
exact_benchmarks = Dict(
    (2, U=0) => -2.8284,   # 2-site, free fermions: -2√2
    (2, U=4) => -1.0,     # 2-site, intermediate coupling
    (4, U=4) => -4.913,  # 4-site, typical benchmark
)

println("Exact reference energies for validation:")
for (system, energy) in exact_benchmarks
    println("  L=$(system[1]), U=$(system[2])t: E₀ = $energy")
end
```

## Advanced Features and Extensions

### 1. Correlation Functions

**Implementation**:
```julia
# Static spin-spin correlations
Szi_szj = zero(Pair{Int,Int})
for i in 1:L, j in 1:L
    Szi_szj = Szi_szj + Sz_operators[i] * Sz_operators[j]
end
correlation_matrix = expectation_value(ψ0_exact, Szi_szj)
```

### 2. Excitation Spectrum

**Implementation**:
```julia
# Multiple eigenstates for gap analysis
eigensystem_multi = eigen(ed; nev=4)
eigen_energies = eigensystem_multi.values
energy_gaps = diff(eigen_energies)
```

### 3. Systematic Convergence

**Implementation**:
```julia
# Finite-size scaling analysis
L_values = [2, 3, 4, 6, 8]
systematic_results = Dict()

for L in L_values
    # Build system
    lattice_L = Lattice(unitcell, (L, 1))  # L×1 chain
    # Compute ground state
    systematic_results[L] = E0_exact
end

# Thermodynamic limit extrapolation
thermodynamic_estimate = extrapolate_to_L_infinity(systematic_results)
```

## Integration with NCTSSoS Ecosystem

### 1. Testing Integration

Reference implementation for NCTSSoS testing:
```julia
# Exact diagonalization for systematic testing
function compute_exact_hubbard_energy(n_sites, U_value)
    # Complete implementation that can be called from tests
    return E0_exact
end
```

### 2. Benchmark Validation

Cross-validation with existing SDP bounds:
```julia
# Compare exact energies with SDP relaxation bounds
@assert E0_exact ≈ expected_sdp_lower_bound "Exact energy should match SDP"
```

### 3. Documentation Cross-Referencing

```markdown
# Links to related examples
- See [MPSKit ground state example](@ref hubbard-mpskit-groundstate) for tensor network approach
- Compare with [Fermionic algebra interface](@ref fermionic-algebra-interface) for SDP relaxation
- References exact benchmark values used in testing
```

## Quality Assurance Plan

### 1. Correctness Verification

- Cross-check energies against known analytical results
- Validate quantum number conservation
- Compare with independent computation (manual diagonalization)
- Verify symmetry properties of eigenstates

### 2. Performance Benchmarking

- Measure memory usage for different system sizes
- Profile computational time vs. Hilbert space dimension
- Compare efficiency with alternative approaches
- Document practical limitations

### 3. Integration Testing

- Test compatibility with existing NCTSSoS examples
- Verify references and cross-links work correctly
- Ensure reproducibility across different environments
- Validate documentation builds properly

## Technical Implementation Details

### 1. Code Organization

```julia
# Main script structure
include("helper_functions.jl")      # Shared utility functions
include("exact_diagonalization.jl")  # Core ED implementation
include("analysis_tools.jl")         # Result processing and visualization
```

### 2. Error Handling

```julia
# Robust error checking
try
    eigensystem = eigen(encoded_ed; nev=1, tolerance=1e-10)
    @assert eigensystem.values[1] ≲ 0 "Ground state energy should be negative"
catch e
    @error "ED failed: $(e)"
    println("System too large for exact diagonalization")
end
```

### 3. Reproducibility

```julia
# Set random seed for reproducible results
using Random
Random.seed!(42)  # Consistent across runs

# Benchmark exact values for testing
const REFERENCE_ENERGIES = Dict(
    (2, 4.0) => -1.0,      # 2 sites, U=4t
    (4, 8.0) => -4.913259,  # 4 sites, U=8t
)
```

This plan provides a comprehensive framework for implementing an educational and practically useful exact diagonalization example that integrates seamlessly with the existing NCTSSoS documentation ecosystem while providing exact benchmark values for validation and testing purposes. The example will demonstrate best practices for using the ExactDiagonalization.jl package while maintaining the high pedagogical standards established in the existing MPSKit example."}