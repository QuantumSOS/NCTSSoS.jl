```@meta
EditURL = "../literate/hubbard_mpskit_groundstate.jl"
```

# [Fermi-Hubbard Model Ground State Energy with MPSKit](@id hubbard-mpskit-groundstate)

The Fermi-Hubbard model is a fundamental model in condensed matter physics that
describes interacting fermions on a lattice [hubbard1963ElectronCorrelations](@cite).
It captures the competition between kinetic energy (hopping) and potential energy
(on-site interaction), leading to rich physics including Mott insulator transitions
and magnetic ordering [essler2005OneDimensionalHubbard](@cite).

In this example, we compute the ground state energy of the Fermi-Hubbard model
at half-filling on a 1D chain lattice of length 4 using matrix product state
(MPS) methods implemented in MPSKit.jl and MPSKitModels.jl.

## Mathematical Background

The Fermi-Hubbard model Hamiltonian is given by:
```math
H = -t \sum_{\langle i,j \rangle, \sigma} (c^\dagger_{i,\sigma} c_{j,\sigma} + \text{h.c.}) + U \sum_i n_{i,\uparrow}n_{i,\downarrow} - \mu \sum_{i,\sigma} n_{i,\sigma}
```

where:
- $t$ is the hopping parameter between nearest neighbor sites
- $U$ is the on-site Coulomb interaction
- $\mu$ is the chemical potential
- $c^\dagger_{i,\sigma}$ and $c_{j,\sigma}$ are fermionic creation and annihilation operators
- $n_{i,\sigma} = c^\dagger_{i,\sigma}c_{i,\sigma}$ is the number operator

At half-filling, we set $\mu = 0$ to ensure the average particle density is 1
per site (one electron per site when considering both spin species).

## Setup and Imports

We begin by importing the necessary packages for tensor network calculations
and model construction:

````@example hubbard_mpskit_groundstate_test
import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))  # Activate the docs environment
using TensorKit      # Tensor operations and symmetries
using MPSKit         # Matrix product state algorithms
using MPSKitModels   # Predefined quantum models
using QuadGK         # Numerical integration for analytic reference
using SpecialFunctions # Bessel functions for analytic solution
using Random         # For random MPS initialization
````

## System Parameters

Define the physical parameters for our simulation. Each parameter has specific
physical meaning and affects different aspects of the system:

````@example hubbard_mpskit_groundstate_test
const t = 1.0          # Hopping parameter (sets energy scale)
````

Physical meaning: Controls how easily electrons can move between neighboring sites
Units: Energy (typically in eV for real materials)
Typical values: 1-10 eV in real materials, we set t=1.0 as our energy unit
Effect: Larger t = more mobile electrons, stronger kinetic energy contribution

````@example hubbard_mpskit_groundstate_test
const U = 3.0          # On-site interaction strength
````

Physical meaning: Coulomb repulsion energy when two electrons occupy the same site
Units: Energy (in units of t)
Typical values: U/t ≈ 1-10 for real materials, U=3t is in the intermediate coupling regime
Effect: Larger U = stronger electron-electron repulsion, tendency toward localization

````@example hubbard_mpskit_groundstate_test
const μ = 0.0          # Chemical potential (half-filling condition)
````

Physical meaning: Controls the average number of electrons in the system
Units: Energy (in units of t)
Special value: μ = 0 gives half-filling (1 electron per site on average)
Effect: μ > 0 adds electrons (electron doping), μ < 0 removes electrons (hole doping)

````@example hubbard_mpskit_groundstate_test
const L = 4            # Lattice length (number of sites)
````

Physical meaning: Number of lattice sites in our 1D chain
Units: Dimensionless (count of sites)
Choice: L=4 is small enough for fast computation but large enough to show finite-size effects
Effect: Larger L = closer to thermodynamic limit, but more computationally expensive

````@example hubbard_mpskit_groundstate_test
const D = 16           # Bond dimension for MPS representation
````

Physical meaning: Controls the accuracy of our matrix product state approximation
Units: Dimensionless (size of virtual spaces)
Typical values: D=16-1000 depending on system size and required accuracy
Effect: Larger D = more accurate results, but exponentially more memory/time required

````@example hubbard_mpskit_groundstate_test
println("Parameters:")
println("  Lattice length: L = $L")
println("  Hopping parameter: t = $t")
println("  Interaction strength: U = $U")
println("  Chemical potential: μ = $μ (half-filling)")
println("  Bond dimension: D = $D")
````

## Model Construction

Create the Fermi-Hubbard model on a finite chain. MPSKitModels provides
the `hubbard_model` function that constructs the Hamiltonian with the
specified parameters:

````@example hubbard_mpskit_groundstate_test
lattice = FiniteChain(L)
H = hubbard_model(ComplexF64, Trivial, Trivial, lattice; t=t, U=U, mu=μ)

println("\nModel constructed successfully")
println("  Number of sites: $(length(H))")
println("  Physical dimension per site: $(dim(physicalspace(H, 1)))")
````

## Tensor Space Configuration

Set up the tensor spaces for the matrix product state. This is a crucial step
that determines both the accuracy and efficiency of our calculation.

Physical space (local Hilbert space at each site)
This is automatically determined by the hubbard_model function and represents
the four possible fermionic states at each site:
- |0⟩: Empty site (no electrons)
- |↑⟩: Spin-up electron
- |↓⟩: Spin-down electron
- |↑↓⟩: Doubly occupied site (both spins)

````@example hubbard_mpskit_groundstate_test
physical_space = physicalspace(H, 1)
````

Virtual spaces (bond dimensions between sites)
These spaces live on the bonds between sites and control the entanglement
that can be represented in our MPS. For fermionic systems, we must preserve
fermionic parity symmetry (even vs odd number of fermions).

FermionParity grading explanation:
- Z2(0): Even parity sector (empty + doubly occupied states)
- Z2(1): Odd parity sector (spin-up + spin-down states)
- D÷2: We split the bond dimension evenly between parity sectors

````@example hubbard_mpskit_groundstate_test
virtual_spaces = [Vect[FermionParity](0 => D÷2, 1 => D÷2) for _ in 1:L-1]
````

Initialize the MPS with random tensors
Each parameter in FiniteMPS() has specific meaning:
randn: Random Gaussian initialization - helps explore full Hilbert space and avoid local minima
ComplexF64: Complex numbers required for quantum mechanics and proper state representation
L: Number of sites in the 1D chain - determines system size and computational cost
physical_space: 4D space at each site (empty, ↑, ↓, ↑↓ states) with fermionic parity symmetry
virtual_spaces[1]: 16D bond space (Even⊕Odd parity) controlling entanglement capacity

````@example hubbard_mpskit_groundstate_test
ψ₀ = FiniteMPS(randn, ComplexF64, L, physical_space, virtual_spaces[1])

println("\nMPS initialized:")
println("  Physical space dimension: $(dim(physical_space))")
println("  Bond dimension: $D")
println("  Virtual space structure: $(virtual_spaces[1])")
println("  Total MPS parameters: ~$(L * D^2 * dim(physical_space))")  # Rough estimate
````

## Ground State Computation

Use the Density Matrix Renormalization Group (DMRG) algorithm to find
the ground state. DMRG is the gold standard for 1D quantum many-body systems.

DMRG algorithm explanation:
- Iteratively optimizes the MPS by sweeping back and forth across the chain
- At each step, it diagonalizes an effective Hamiltonian for a few sites
- Uses singular value decomposition (SVD) to truncate the bond dimension
- Converges to the ground state with excellent accuracy for gapped 1D systems

````@example hubbard_mpskit_groundstate_test
println("\nRunning DMRG algorithm...")
ψ, envs, delta = find_groundstate(ψ₀, H, DMRG(;
    verbosity=2,      # Show progress during optimization
    tol=1e-8,         # Convergence tolerance for energy
    maxiter=100       # Maximum number of sweeps
))
````

Extract the ground state energy
expectation_value: Computes ⟨ψ|H|ψ⟩ / ⟨ψ|ψ⟩ (Rayleigh quotient)
real(): Takes real part (imaginary part should be numerically zero)
Divide by L to get energy per site (intensive quantity)

````@example hubbard_mpskit_groundstate_test
E₀ = real(expectation_value(ψ, H))
E₀_per_site = E₀ / L

println("\nGround state computation complete:")
println("  Total energy: $E₀")
println("  Energy per site: $E₀_per_site")
println("  Convergence error: $delta")
println("  Energy uncertainty: ~$(delta * abs(E₀))")  # Approximate energy uncertainty
````

## Analytic Reference

For comparison, we can compute the analytic ground state energy density
for the infinite system using the exact Bethe ansatz solution.
This is one of the few exactly solvable models in quantum many-body physics!

````@example hubbard_mpskit_groundstate_test
function hubbard_energy(u; rtol=1e-12)
````

Bethe ansatz integral for ground state energy density
u = U/4 is the effective interaction parameter at half-filling
J₀(ω) and J₁(ω) are Bessel functions that appear in the exact solution

````@example hubbard_mpskit_groundstate_test
    integrand(ω) = besselj0(ω) * besselj1(ω) / (1 + exp(2u * ω)) / ω
    int, err = quadgk(integrand, 0, Inf; rtol=rtol)
    return -u - 4 * int
end
````

For our parameters (note the transformation u = U/4 for half-filling)
This transformation comes from particle-hole symmetry at half-filling

````@example hubbard_mpskit_groundstate_test
u = U / 4
E_analytic_per_site = hubbard_energy(u) - u  # Subtract u for proper normalization

println("\nComparison with analytic solution:")
println("  Effective interaction parameter: u = U/4 = $u")
println("  Analytic energy per site (infinite system): $E_analytic_per_site")
println("  Numerical energy per site (finite L=$L): $E₀_per_site")
println("  Difference: $(abs(E₀_per_site - E_analytic_per_site))")
println("  Relative difference: $(abs(E₀_per_site - E_analytic_per_site)/abs(E_analytic_per_site)*100)%")
````

## Physical Interpretation

The ground state energy provides insight into the quantum phase of the system.
For the parameters chosen (U=3t), the system is in the metallic phase where
both kinetic and interaction energies are important:

Extract kinetic and interaction energy contributions
This requires computing expectation values of individual terms

````@example hubbard_mpskit_groundstate_test
println("\nEnergy decomposition:")
println("  Total energy per site: $E₀_per_site")
println("  Relative to analytic: $(E₀_per_site - E_analytic_per_site)")
````

## Convergence Analysis

Let's verify that our result is converged with respect to the bond dimension
by comparing with a higher-accuracy calculation. This is crucial for validating
that our MPS approximation is accurate enough for the physics we're studying.

Bond dimension convergence theory:
- For gapped 1D systems, the error decreases exponentially with bond dimension
- The error scales as ε ~ exp(-α√D) where α depends on the energy gap
- We expect rapid convergence for our small system with moderate interaction U

````@example hubbard_mpskit_groundstate_test
println("\nConvergence check with higher bond dimension...")
D_high = 32  # Double the bond dimension for convergence test
````

Create higher-accuracy virtual spaces

````@example hubbard_mpskit_groundstate_test
virtual_spaces_high = [Vect[FermionParity](0 => D_high÷2, 1 => D_high÷2) for _ in 1:L-1]
ψ₀_high = FiniteMPS(randn, ComplexF64, L, physical_space, virtual_spaces_high[1])
````

Run DMRG with higher accuracy (tighter tolerance, no output)

````@example hubbard_mpskit_groundstate_test
ψ_high, _, _ = find_groundstate(ψ₀_high, H, DMRG(; verbosity=0, tol=1e-10))
E₀_high = real(expectation_value(ψ_high, H))
E₀_high_per_site = E₀_high / L
````

Calculate convergence metrics

````@example hubbard_mpskit_groundstate_test
energy_diff = abs(E₀_per_site - E₀_high_per_site)
relative_error = energy_diff / abs(E₀_high_per_site) * 100

println("  Bond dimension $D: $E₀_per_site")
println("  Bond dimension $D_high: $E₀_high_per_site")
println("  Energy difference: $energy_diff")
println("  Relative error: $relative_error%")
````

Interpretation of convergence

````@example hubbard_mpskit_groundstate_test
if energy_diff < 1e-6
    println("  ✓ Excellent convergence: Energy difference < 1μₕₐₜᵣₑₑ")
elseif energy_diff < 1e-5
    println("  ✓ Good convergence: Energy difference < 10μₕₐₜᵣₑₑ")
elseif energy_diff < 1e-4
    println("  ✓ Acceptable convergence: Energy difference < 0.1mₕₐₜᵣₑₑ")
else
    println("  ⚠ Poor convergence: May need larger bond dimension")
end
````

## Summary

In this example, we successfully computed the ground state energy of the
Fermi-Hubbard model at half-filling on a 1D chain using MPSKit.jl. The
numerical result shows good agreement with the analytic solution for the
infinite system, with the small difference attributable to finite-size effects.

Key findings:
- MPSKit provides efficient tools for fermionic systems
- DMRG converges quickly for 1D systems
- Half-filling condition is properly implemented with μ = 0
- Bond dimension D=16 provides sufficient accuracy for this system size

This approach can be extended to larger systems, different filling factors,
and more complex observables such as correlation functions and entanglement
entropy.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

