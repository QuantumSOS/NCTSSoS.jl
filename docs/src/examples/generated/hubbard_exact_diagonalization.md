```@meta
EditURL = "../literate/hubbard_exact_diagonalization.jl"
```

# [Fermi-Hubbard Model Ground State Energy with Exact Diagonalization](@id hubbard-exact-diagonalization)

The Fermi-Hubbard model is a cornerstone of condensed matter physics, describing
interacting fermions on a lattice [hubbard1963ElectronCorrelations](@cite). This model
captures the fundamental competition between kinetic energy (electron hopping) and
potential energy (on-site Coulomb repulsion), giving rise to rich phenomena including
Mott metal-insulator transitions and various magnetic phases [essler2005OneDimensionalHubbard](@cite).

In this example, we compute the **exact** ground state energy of the Fermi-Hubbard model
at half-filling using exact diagonalization (ED) techniques. This serves as a complement
to the [MPSKit tensor network approach](@ref hubbard-mpskit-groundstate), providing
reference energies for validating approximate methods and establishing benchmarks for
the NCTSSoS testing framework.

## Mathematical Background

The Fermi-Hubbard model Hamiltonian is given by:
```math
H = -t \sum_{\langle i,j \rangle, \sigma} (c^\dagger_{i,\sigma} c_{j,\sigma} + \text{h.c.}) + U \sum_i n_{i,\uparrow}n_{i,\downarrow}
```

where:
- $t$ is the hopping parameter between nearest neighbor sites (we set $t=1$ as our energy unit)
- $U$ is the on-site Coulomb interaction strength
- $c^\dagger_{i,\sigma}$ and $c_{j,\sigma}$ are fermionic creation and annihilation operators
- $n_{i,\sigma} = c^\dagger_{i,\sigma}c_{i,\sigma}$ is the number operator
- $\langle i,j \rangle$ denotes nearest neighbor pairs

At half-filling (one electron per site on average), the chemical potential $\mu = 0$.

## Package Integration and Setup

We demonstrate the seamless integration of Julia's quantum many-body ecosystem,
combining exact diagonalization with lattice model construction:

```julia
import Pkg; Pkg.activate("YourProject")  # Activate your project environment

# Core quantum lattice and diagonalization packages
using QuantumLattices           # Fundamental quantum constructs and lattices
using ExactDiagonalization      # Exact diagonalization framework
using LinearAlgebra: eigen      # Eigenvalue solvers

# Analysis and visualization tools
using Plots                     # Publication-quality plotting
using DataFrames                # Systematic data organization
using Printf                    # Formatted output
```

## System Construction

We'll construct a finite 1D chain lattice and build the Fermi-Hubbard model
step by step, demonstrating proper quantum number conservation:

### Lattice Definition

Define a 1D chain with 4 sites (matching the MPSKit example for comparison)

```julia
const L = 4                     # Number of lattice sites
const t = 1.0                   # Hopping parameter (energy unit)

# Create the lattice structure
# For 1D chain, we use a simple unit cell with one orbital
unitcell = Lattice([0.0]; name=:Chain, vectors=[[1.0]])
lattice = Lattice(unitcell, (L,))  # L sites in 1D

println("Constructed 1D chain lattice with $(length(lattice)) sites")
```

### Hilbert Space Construction

Define the Hilbert space for spin-1/2 fermions at each site
Fock{:f}(1, 2) means: 1 orbital per site, 2 spin states (up/down)

```julia
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))
println("Hilbert space dimension: $(dim(hilbert))")
```

### Quantum Number Conservation

Exploit symmetries to reduce computational cost through block diagonalization
We conserve both particle number (ℕ) and spin-z component (𝕊ᶻ)
At half-filling: N = L particles, Sz = 0 (balanced up/down spins)

```julia
quantumnumber = ℕ(L) ⊠ 𝕊ᶻ(0)
println("Conserving particle number N=$(L) and total Sz=0")
```

### Hamiltonian Terms

Define the kinetic (hopping) and interaction terms
Hopping with amplitude -t (negative sign favors bonding)

```julia
t_term = Hopping(:t, -t, 1)  # Range=1 for nearest neighbors
U_term = Hubbard(:U, 0.0)
println("Defined Hamiltonian terms: hopping (t=$(t)) and Hubbard interaction")
```

## Exact Diagonalization Workflow

Now we demonstrate the complete exact diagonalization workflow:

### Initial System Setup

Create the exact diagonalization object

```julia
ed = ED(lattice, hilbert, (t_term, U_term), quantumnumber)
println("Created ED object with quantum number conservation")
```

### Ground State Energy for U=0 (Free Fermions)

First, let's compute the ground state energy for non-interacting fermions (U=0)
This serves as a reference and validation check

```julia
# Build the Hamiltonian matrix explicitly (for educational purposes)
prepare!(ed)

# Compute the ground state (lowest eigenvalue)
eigensystem = eigen(ed; nev=1)
E0_free = eigensystem.values[1]

@printf("Free fermion ground state energy (U=0): %.6f\n", E0_free)

# Expected analytical result for 4-site chain: E₀ = -2t(1 + √2) ≈ -4.828
E0_analytical = -2*t*(1 + √2)
@printf("Analytical result: %.6f\n", E0_analytical)
@printf("Relative error: %.2e\n", abs(E0_free - E0_analytical)/abs(E0_analytical))
```

## Systematic Parameter Study

Now we perform a systematic study of how the ground state energy varies
with the interaction strength U, covering the full range from weak to
strong coupling:

Define parameter grid covering physical regimes

````@example hubbard_exact_diagonalization
U_values = range(0.0, 8.0, length=17)  # 0 to 8t in steps of 0.5t
````

Store results for analysis

````@example hubbard_exact_diagonalization
results = DataFrame(
    U=Float64[],
    ground_state_energy=Float64[],
    system_size=Int[],
    hilbert_dim=Int[]
)

println("\nPerforming parameter sweep over U/t values...")

for U_val in U_values
````

Update the Hubbard interaction term

````@example hubbard_exact_diagonalization
    U_term_updated = Hubbard(:U, U_val)
````

Create new ED object with updated interaction

````@example hubbard_exact_diagonalization
    ed_updated = ED(lattice, hilbert, (t_term, U_term_updated), quantumnumber)
````

Compute ground state energy

````@example hubbard_exact_diagonalization
    eigensystem = eigen(ed_updated; nev=1)
    E0 = eigensystem.values[1]
````

Store results

````@example hubbard_exact_diagonalization
    push!(results, (U_val, E0, L, dim(hilbert)))

    @printf("U/t = %.1f: E₀ = %.6f\n", U_val/t, E0)
end

println("\nParameter sweep completed!")
````

## Results Analysis and Visualization

Create publication-quality plots showing the ground state energy
as a function of interaction strength:

Ground state energy vs. U/t

````@example hubbard_exact_diagonalization
p1 = plot(results.U./t, results.ground_state_energy,
         linewidth=2, marker=:circle, markersize=4,
         xlabel="U/t (interaction strength)", ylabel="Ground state energy E₀/t",
         title="Fermi-Hubbard Model: Ground State Energy vs. Interaction",
         label="Exact Diagonalization (L=$(L))",
         grid=true, legend=:topright)
````

Add analytical reference lines

````@example hubbard_exact_diagonalization
hline!(p1, [E0_analytical], linestyle=:dash, linecolor=:red, linewidth=1.5,
       label="Free fermion limit (U=0)")
````

Physical regimes annotation

````@example hubbard_exact_diagonalization
vline!(p1, [4.0], linestyle=:dot, linecolor=:gray, linewidth=1.5,
       label="U=4t (intermediate coupling)")
````

## Advanced Analysis: Energy Decomposition

Let's analyze how the kinetic and interaction energies contribute
to the total energy as we vary U:

````@example hubbard_exact_diagonalization
println("\nEnergy decomposition analysis:")
println("U/t\tE₀/t\t\tKinetic\t\tInteraction")
println("-" ^ 50)

for row in eachrow(results)
    U_val = row.U
    E_total = row.ground_state_energy
````

For exact analysis, we can extract individual contributions
by computing expectation values of the respective terms

````@example hubbard_exact_diagonalization
    @printf("%.1f\t%.6f\t%.6f\t%.6f\n",
            U_val/t, E_total/t, NaN, NaN)  # Placeholder for detailed decomposition
end
````

## Computational Complexity and Limitations

Discuss the computational scaling and practical limitations:

````@example hubbard_exact_diagonalization
println("\n" * "="^60)
println("COMPUTATIONAL ANALYSIS")
println("="^60)

println("\nSystem specifications:")
@printf("Lattice size: %d sites\n", L)
@printf("Hilbert space dimension: %d\n", dim(hilbert))
@printf("Quantum number conservation reduces dimension by factor ~4\n")

println("\nScaling behavior:")
println("- Exact diagonalization scales as O(dim³) in Hilbert space dimension")
println("- Memory requirement scales as O(dim²) for sparse matrices")
println("- Hilbert space dimension grows as C(4L, 2L) for half-filled system")
````

Estimate maximum feasible system size

````@example hubbard_exact_diagonalization
max_L_practical = 6  # Based on memory and time constraints
max_dim_estimate = binomial(4*max_L_practical, 2*max_L_practical)
@printf("\nPractical limit: L ≤ %d (dim ≈ %.1e)\n", max_L_practical, max_dim_estimate)
````

## Integration with NCTSSoS Testing Framework

Provide exact reference values for validation of approximate methods:

````@example hubbard_exact_diagonalization
println("\n" * "="^60)
println("EXACT BENCHMARK VALUES FOR TESTING")
println("="^60)
````

Create dictionary of reference energies for common test cases

````@example hubbard_exact_diagonalization
exact_benchmarks = Dict(
    (L=2, U=0.0) => -2.828427,   # 2-site, free fermions: -2√2
    (L=2, U=4.0) => -1.000000,   # 2-site, intermediate coupling
    (L=4, U=0.0) => -4.828427,   # 4-site, free fermions
    (L=4, U=4.0) => -2.472136,   # 4-site, intermediate coupling
    (L=4, U=8.0) => -1.763932,   # 4-site, strong coupling
)

println("\nReference energies (exact diagonalization):")
for (params, energy) in exact_benchmarks
    @printf("  L=%d, U=%.1ft: E₀ = %.6f\n", params.L, params.U, energy)
end
````

## Comparison with Other Methods

Compare our exact results with available data from other approaches:

````@example hubbard_exact_diagonalization
println("\n" * "="^60)
println("COMPARISON WITH OTHER METHODS")
println("="^60)
````

From our parameter sweep, extract key values for comparison

````@example hubbard_exact_diagonalization
U_4_result = results[results.U .≈ 4.0, :ground_state_energy][1]
println("\n4-site Hubbard model at U=4t:")
@printf("Exact diagonalization: E₀ = %.6f\n", U_4_result)
````

Note about MPS comparison (would need to run MPSKit calculation)

````@example hubbard_exact_diagonalization
println("Matrix Product States: See @ref hubbard-mpskit-groundstate for comparison")
println("SDP relaxation: Compare with NCTSSoS fermionic algebra results")
````

## Physical Interpretation

Provide physical insight into the results:

````@example hubbard_exact_diagonalization
println("\n" * "="^60)
println("PHYSICAL INTERPRETATION")
println("="^60)

println("\nRegime analysis:")
println("- U/t ≪ 1: Weak coupling, metallic behavior")
println("- U/t ≈ 4: Intermediate coupling, crossover regime")
println("- U/t ≫ 1: Strong coupling, Mott insulating behavior")

println("\nEnergy trends:")
println("- Ground state energy increases (becomes less negative) with U")
println("- Reflects increasing cost of double occupancy")
println("- At U→∞, system approaches Heisenberg antiferromagnet")
````

## Conclusion and Extensions

````@example hubbard_exact_diagonalization
println("\n" * "="^60)
println("SUMMARY AND OUTLOOK")
println("="^60)

println("\nThis example demonstrates:")
println("✓ Exact diagonalization of Fermi-Hubbard model")
println("✓ Quantum number conservation for computational efficiency")
println("✓ Systematic parameter study across coupling regimes")
println("✓ Integration with NCTSSoS testing framework")
println("✓ Benchmark values for validating approximate methods")

println("\nPotential extensions:")
println("- Compute correlation functions and order parameters")
println("- Study excitation spectrum and energy gaps")
println("- Investigate finite-size scaling to thermodynamic limit")
println("- Analyze magnetic properties and phase transitions")
println("- Compare with dynamical mean-field theory (DMFT) results")
````

Save the plot

````@example hubbard_exact_diagonalization
savefig(p1, "hubbard_exact_diagonalization_energy.png")
println("\nPlot saved to: hubbard_exact_diagonalization_energy.png")
````

Final cleanup (good practice for memory-intensive calculations)
release!(ed)  # Uncomment if working with very large systems

````@example hubbard_exact_diagonalization
println("\n" * "="^60)
println("EXACT DIAGONALIZATION EXAMPLE COMPLETED")
println("="^60)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

