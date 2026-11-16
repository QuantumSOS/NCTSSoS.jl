```@meta
EditURL = "../literate/pauli_algebra_interface.jl"
```

# [Simplified Quantum Spin Models with Pauli Algebra](@id pauli-algebra-interface)

`NCTSSoS.jl` provides a convenient interface for working with common quantum algebras,
eliminating the need to manually specify commutation relations and constraints.
This tutorial demonstrates the `pauli_variables` function for quantum spin systems,
which significantly simplifies the problem setup for [polynomial optimization](@ref polynomials)
problems in quantum many-body physics [wang2024Certifying](@cite).

## The Problem: Manual Constraint Specification

If we used free variables for quantum spins, we would have to specify
every Pauli relation. This is error-prone and obscures the physics.
At each site, the six ordered products of distinct axes must have the
correct signs, and each axis must square to the identity. Operators on
different sites must commute.

### Heisenberg XXX Model

The Heisenberg XXX Hamiltonian for a 1D chain with periodic boundary conditions is:

```math
H = \frac{1}{4} \sum_{i=1}^{N} \left( \sigma_i^x \sigma_{i+1}^x + \sigma_i^y \sigma_{i+1}^y + \sigma_i^z \sigma_{i+1}^z \right)
```

where $\sigma_i^{x,y,z}$ are the Pauli operators at site $i$.
Six ordered-product relations at each of six sites would require
**36 equations**, before the square and cross-site relations. For example,
$\sigma_i^x\sigma_i^y=i\sigma_i^z$ but
$\sigma_i^y\sigma_i^x=-i\sigma_i^z$. A wrong sign changes the problem.

## The Solution: Typed Algebra Variables (Recommended API)

The `pauli_variables` function encapsulates all these constraints
automatically through Julia's type system, allowing you to focus on the physics
rather than the algebra.

### Simplified Approach with `pauli_variables`

The same problem can be solved much more concisely:

!!! note "SDP Solver"
    This example uses [Mosek](https://www.mosek.com/) through `MosekBackend`.
    For open-source alternatives, see [Optimizers](@ref overview-of-optimizers).

````julia
using NCTSSoS, Mosek
N = 6  # Number of spins in the chain

registry, (σx, σy, σz) = pauli_variables(1:N)
````

````
(Registry(Pauli, 6 sites, 18 variables), ([X1, X2, X3, X4, X5, X6], [Y1, Y2, Y3, Y4, Y5, Y6], [Z1, Z2, Z3, Z4, Z5, Z6]))
````

Check the relations directly; there is no constraint list to maintain.

````julia
@assert σx[1] * σy[1] == im * σz[1]
@assert σy[1] * σx[1] == -im * σz[1]
@assert σx[1]^2 == one(σx[1])
@assert σx[1] * σy[2] == σy[2] * σx[1]
````

Construct the Hamiltonian (same formula, cleaner variables)

````julia
ham = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)]
          for op in [σx, σy, σz] for i in 1:N)
````

````
(0.25) Z1 Z2 + (0.25) Z2 Z3 + (0.25) Z3 Z4 + (0.25) Z4 Z5 + (0.25) Z1 Z6 + (0.25) Z5 Z6 + (0.25) X1 X2 + (0.25) Y1 Y2 + … (18 terms)
````

Create the optimization problem - constraints are handled automatically by the algebra type!

````julia
problem = Problem(ham)
````

````
Problem: minimize (0.25) Z1 Z2 + (0.25) Z2 Z3 + (0.25) Z3 Z4 + (0.25) Z4 Z5 + (0.25) Z1 Z6 + (0.25) Z5 Z6 + (0.25) X1 X2 + (0.25) Y1 Y2 + … (18 terms); 0 operator constraints
````

The algebra supplies the relations, so no manual constraint equations are needed.

## Computing Ground State Energy

Now we can solve for the ground state energy lower bound using the
`solve` function. We configure it to use a second-order
[moment relaxation](@ref moment-sohs-hierarchy) [wang2024Certifying](@cite).

````julia
relaxation = Relaxation(order=2)
solver = MosekBackend()

res = solve(problem, relaxation; solver)
energy_per_site = res.bound / N
````

````
-0.4671292756829241
````

The result provides a certified lower bound on the ground state energy per site.
For the 6-site XXX Heisenberg chain, this yields approximately **-0.467129**,
which matches the exact value to high precision [wang2024Certifying](@cite).

## Advantages of Typed Algebra Variables

The `pauli_variables` interface provides several key benefits:

1. **Automatic constraint generation**: All Pauli commutation relations are
   encoded by the `Pauli` algebra without manual specification.

2. **Error prevention**: Eliminates typos and sign errors in constraint equations.

3. **Code clarity**: The physics intent is immediately clear from the function name.

4. **Scalability**: Works seamlessly for any system size without code modification.

5. **Consistency**: Ensures the same algebraic structure across different problems.

6. **Extensibility**: Can still add custom constraints when needed for specific
   physical scenarios.

## Other Algebra Types

Similar functions exist for other quantum algebras:

- `fermionic_variables`: Fermionic operators with anticommutation relations
- `bosonic_variables`: Bosonic operators with commutation relations
- `site_variables(Projector(), nsites, nletters)`: Projectors (P² = P)
- `site_variables(Unipotent(), nsites, nletters)`: Self-adjoint operators with U² = I
- `site_variables(Free(), nsites, nletters)`: Free self-adjoint operators within each site

## NCTSSoS Steps

This interface extends naturally to other quantum systems:

- For systems with more complex geometries, see the [2D lattice examples](@ref ground-state-energy)
- For correlation function bounds, see [certifying ground state properties](@ref certify-property)
- For non-local correlations, explore [Bell inequalities](@ref bell-inequalities)

The typed algebra approach demonstrates how `NCTSSoS.jl` bridges the gap
between physical intuition and mathematical formalism, making quantum many-body
optimization more accessible and reliable.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

