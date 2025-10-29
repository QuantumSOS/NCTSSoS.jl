# # [Simplified Quantum Spin Models with Pauli Algebra Interface](@id pauli-algebra-interface)

# `NCTSSoS.jl` provides a convenient interface for working with common quantum algebras,
# eliminating the need to manually specify commutation relations and constraints.
# This tutorial demonstrates the [`pauli_algebra`](@ref) constructor for quantum spin systems,
# which significantly simplifies the problem setup for [polynomial optimization](@ref polynomials)
# problems in quantum many-body physics [wang2024Certifying](@cite).

# ## The Problem: Manual Constraint Specification

# When working with quantum spin systems, the standard approach requires manually
# defining all Pauli operator commutation relations. This is tedious, error-prone,
# and obscures the physics of the problem. Let's see this with a concrete example.

# ### Traditional Approach: Heisenberg XXX Model

# The Heisenberg XXX Hamiltonian for a 1D chain with periodic boundary conditions is:
#
# ```math
# H = \frac{1}{4} \sum_{i=1}^{N} \left( \sigma_i^x \sigma_{i+1}^x + \sigma_i^y \sigma_{i+1}^y + \sigma_i^z \sigma_{i+1}^z \right)
# ```
#
# where $\sigma_i^{x,y,z}$ are the Pauli operators at site $i$. To solve for the
# [ground state energy](@ref ground-state-energy) using polynomial optimization,
# we traditionally need to:

using NCTSSoS, MosekTools
N = 6  # Number of spins in the chain

# **Step 1**: Declare non-commutative variables for Pauli operators
@ncpolyvar x[1:N] y[1:N] z[1:N]

# **Step 2**: Construct the Hamiltonian
ham = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)] for op in [x, y, z] for i in 1:N)

# **Step 3**: Manually specify all Pauli commutation relations
# For each site i, we need to encode:
# [σ^x_i, σ^y_i] = iσ^z_i, [σ^y_i, σ^z_i] = iσ^x_i, [σ^z_i, σ^x_i] = iσ^y_i
eq_cons = reduce(vcat, [
    [x[i] * y[i] - im * z[i],   # σ^x_i σ^y_i = iσ^z_i
     y[i] * x[i] + im * z[i],   # σ^y_i σ^x_i = -iσ^z_i
     y[i] * z[i] - im * x[i],   # σ^y_i σ^z_i = iσ^x_i
     z[i] * y[i] + im * x[i],   # σ^z_i σ^y_i = -iσ^x_i
     z[i] * x[i] - im * y[i],   # σ^z_i σ^x_i = iσ^y_i
     x[i] * z[i] + im * y[i]]   # σ^x_i σ^z_i = -iσ^y_i
    for i in 1:N
])

# **Step 4**: Create the optimization problem with all constraints
pop_old = cpolyopt(ham;
    eq_constraints=eq_cons,                     # Pauli commutation relations
    comm_gps=[[x[i], y[i], z[i]] for i in 1:N], # Operators on different sites commute
    is_unipotent=true                           # Pauli operators square to identity
)

# This approach requires **36 constraint equations** for just 6 spins! As the system
# size grows, this becomes increasingly cumbersome and error-prone.

# ## The Solution: Algebra Constructors

# The [`pauli_algebra`](@ref) constructor encapsulates all these constraints
# automatically, allowing you to focus on the physics rather than the algebra.

# ### Simplified Approach with `pauli_algebra`

# The same problem can be solved much more concisely:

sys = pauli_algebra(N)
x_new, y_new, z_new = sys.variables

# Construct the Hamiltonian (same formula, different variables)
ham_new = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)]
              for op in [x_new, y_new, z_new] for i in 1:N)

# Create the optimization problem - constraints are automatic!
pop_new = cpolyopt(ham_new, sys)

# Both approaches produce identical optimization problems, but the new interface is
# **10 lines shorter** and eliminates the possibility of typos in constraint equations.

# ## Computing Ground State Energy

# Now we can solve for the ground state energy lower bound using the
# [`cs_nctssos`](@ref) solver. We configure it to use a second-order
# [moment relaxation](@ref moment-sohs-hierarchy) [wang2024Certifying](@cite).

solver_config = SolverConfig(
    optimizer=Mosek.Optimizer,  # SDP solver backend
    order=2                     # Relaxation order (higher = tighter bound)
)

res = cs_nctssos(pop_new, solver_config)
energy_per_site = res.objective / N

# The result provides a certified lower bound on the ground state energy per site.
# For the 6-site XXX Heisenberg chain, this yields approximately **-0.467129**,
# which matches the exact value to high precision [wang2024Certifying](@cite).


# ## Advantages of Algebra Constructors

# The [`pauli_algebra`](@ref) interface provides several key benefits:

# 1. **Automatic constraint generation**: All Pauli commutation relations are
#    encoded correctly without manual specification.
#
# 2. **Error prevention**: Eliminates typos and sign errors in constraint equations.
#
# 3. **Code clarity**: The physics intent is immediately clear from `pauli_algebra(N)`.
#
# 4. **Scalability**: Works seamlessly for any system size without code modification.
#
# 5. **Consistency**: Ensures the same algebraic structure across different problems.
#
# 6. **Extensibility**: Can still add custom constraints when needed for specific
#    physical scenarios.

# ## Next Steps

# This interface extends naturally to other quantum systems:
#
# - For systems with more complex geometries, see the [2D lattice examples](@ref ground-state-energy)
# - For correlation function bounds, see [certifying ground state properties](@ref certify-property)
# - For non-local correlations, explore [Bell inequalities](@ref bell-inequalities)
#
# The algebra constructor approach demonstrates how `NCTSSoS.jl` bridges the gap
# between physical intuition and mathematical formalism, making quantum many-body
# optimization more accessible and reliable.
