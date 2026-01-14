# Simplified Quantum Spin Models with Pauli Algebra

`NCTSSoS.jl` provides a convenient interface for working with common quantum algebras,
eliminating the need to manually specify commutation relations and constraints.
This tutorial demonstrates the `create_pauli_variables` function for quantum spin systems,
which significantly simplifies the problem setup for polynomial optimization
problems in quantum many-body physics [wang2024Certifying](@cite).

## The Problem: Manual Constraint Specification (Legacy API - Deprecated)

> **Deprecated API**
>
> The following example demonstrates the **legacy API** using `@ncpolyvar` and manual
> constraint specification. This API is deprecated and will be removed in a future version.
> Please use the new typed algebra API demonstrated in the next section.

When working with quantum spin systems, the old approach required manually
defining all Pauli operator commutation relations. This was tedious, error-prone,
and obscured the physics of the problem. Let's see this with a concrete example.

### Traditional Approach: Heisenberg XXX Model (Deprecated)

The Heisenberg XXX Hamiltonian for a 1D chain with periodic boundary conditions is:

```math
H = \frac{1}{4} \sum_{i=1}^{N} \left( \sigma_i^x \sigma_{i+1}^x + \sigma_i^y \sigma_{i+1}^y + \sigma_i^z \sigma_{i+1}^z \right)
```

where $\sigma_i^{x,y,z}$ are the Pauli operators at site $i$. In the legacy API,
solving for the ground state energy using polynomial
optimization required:

```julia
# DEPRECATED - Do not use in new code
using NCTSSoS, MosekTools
N = 6  # Number of spins in the chain

# Step 1: Declare non-commutative variables for Pauli operators
@ncpolyvar x[1:N] y[1:N] z[1:N]  # Requires DynamicPolynomials

# Step 2: Construct the Hamiltonian
ham = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)] for op in [x, y, z] for i in 1:N)

# Step 3: Manually specify all 36 Pauli commutation relations (!)
eq_cons = reduce(vcat, [
    [x[i] * y[i] - im * z[i],   # sigma_x * sigma_y = i*sigma_z
     y[i] * x[i] + im * z[i],   # sigma_y * sigma_x = -i*sigma_z
     ...  # and so on for all 6 relations at each site
    ]
    for i in 1:N
])

# Step 4: Create the optimization problem with all constraints
pop_old = cpolyopt(ham;
    eq_constraints=eq_cons,                     # Pauli commutation relations
    comm_gps=[[x[i], y[i], z[i]] for i in 1:N], # Operators on different sites commute
    is_unipotent=true                           # Pauli operators square to identity
)
```

This approach required **36 constraint equations** for just 6 spins! As the system
size grew, this became increasingly cumbersome and error-prone.

## The Solution: Typed Algebra Variables (Recommended API)

The `create_pauli_variables` function encapsulates all these constraints
automatically through Julia's type system, allowing you to focus on the physics
rather than the algebra.

### Simplified Approach with `create_pauli_variables`

The same problem can be solved much more concisely:

````julia
using NCTSSoS, MosekTools
N = 6  # Number of spins in the chain

registry, (σx, σy, σz) = create_pauli_variables(1:N)
````

Construct the Hamiltonian (same formula, cleaner variables)

````julia
ham = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)]
          for op in [σx, σy, σz] for i in 1:N)
````

Create the optimization problem - constraints are handled automatically by the algebra type!

````julia
pop = polyopt(ham, registry)
````

The new interface is much cleaner: just 4 lines instead of 40+, and eliminates
the possibility of typos in constraint equations.

## Computing Ground State Energy

Now we can solve for the ground state energy lower bound using the
`cs_nctssos` solver. We configure it to use a second-order
moment relaxation [wang2024Certifying](@cite).

````julia
solver_config = SolverConfig(
    optimizer=Mosek.Optimizer,  # SDP solver backend
    order=2                     # Relaxation order (higher = tighter bound)
)

res = cs_nctssos(pop, solver_config)
energy_per_site = res.objective / N
````

The result provides a certified lower bound on the ground state energy per site.
For the 6-site XXX Heisenberg chain, this yields approximately **-0.467129**,
which matches the exact value to high precision [wang2024Certifying](@cite).

## Advantages of Typed Algebra Variables

The `create_pauli_variables` interface provides several key benefits:

1. **Automatic constraint generation**: All Pauli commutation relations are
   encoded correctly through the `PauliAlgebra` type without manual specification.

2. **Error prevention**: Eliminates typos and sign errors in constraint equations.

3. **Code clarity**: The physics intent is immediately clear from the function name.

4. **Scalability**: Works seamlessly for any system size without code modification.

5. **Consistency**: Ensures the same algebraic structure across different problems.

6. **Extensibility**: Can still add custom constraints when needed for specific
   physical scenarios.

## Other Algebra Types

Similar functions exist for other quantum algebras:

- `create_fermionic_variables`: Fermionic operators with anticommutation relations
- `create_bosonic_variables`: Bosonic operators with commutation relations
- `create_projector_variables`: Projector operators (P² = P)
- `create_unipotent_variables`: Unipotent operators (U² = I)
- `create_noncommutative_variables`: Generic non-commutative variables

## Next Steps

This interface extends naturally to other quantum systems:

- For systems with more complex geometries, see the 2D lattice examples
- For correlation function bounds, see certifying ground state properties
- For non-local correlations, explore Bell inequalities

The typed algebra approach demonstrates how `NCTSSoS.jl` bridges the gap
between physical intuition and mathematical formalism, making quantum many-body
optimization more accessible and reliable.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

