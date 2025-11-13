# # [Simplified Fermionic Systems with Fermionic Algebra Interface](@id fermionic-algebra-interface)

# `NCTSSoS.jl` provides a convenient interface for working with fermionic quantum systems,
# eliminating the need to manually specify anti-commutation relations and nilpotency constraints.
# This tutorial demonstrates the [`fermionic_algebra`](@ref) constructor for fermionic many-body systems,
# which significantly simplifies the problem setup for [polynomial optimization](@ref polynomials)
# problems involving fermions (electrons, quarks, etc.).

# ## The Problem: Manual Constraint Specification

# When working with fermionic systems, the standard approach requires manually
# defining all fermionic operator anti-commutation relations and nilpotency constraints.
# This is tedious, error-prone, and obscures the physics of the problem.

# ### Traditional Approach: Free Fermion System

# Consider a simple system of N fermionic modes with creation operators $c^\dagger_i$
# and annihilation operators $c_i$. These operators must satisfy:
#
# ```math
# \begin{aligned}
# \{c_i, c_j\} &= 0, \quad \{c^\dagger_i, c^\dagger_j\} = 0 \\
# \{c_i, c^\dagger_j\} &= \delta_{ij} \\
# c_i^2 &= 0, \quad (c^\dagger_i)^2 = 0
# \end{aligned}
# ```
#
# where $\{A, B\} = AB + BA$ denotes the anti-commutator. To study the
# ground state energy of a free fermion Hamiltonian using polynomial optimization,
# we traditionally need to:

using NCTSSoS, MosekTools
N = 2  # Number of fermionic modes

# **Step 1**: Declare non-commutative variables for fermionic operators
@ncpolyvar c[1:N] c_dag[1:N]

# **Step 2**: Construct a simple Hamiltonian (single-particle energies)
# H = ε₁ n₁ + ε₂ n₂, where nᵢ = c†ᵢ cᵢ is the number operator
ε = [-1.0, -0.5]  # Single-particle energies
# Symmetrized for SDP compatibility
ham_manual = sum(ComplexF64(ε[i]/2) * (c_dag[i]*c[i] + c[i]*c_dag[i]) for i in 1:N)

# **Step 3**: Manually specify all anti-commutation relations and nilpotency
eq_cons_manual = []

# Anti-commutation: {cᵢ, cⱼ} = 0
for i in 1:N, j in i:N
    push!(eq_cons_manual, ComplexF64(1.0) * (c[i]*c[j] + c[j]*c[i]))
end

# Anti-commutation: {c†ᵢ, c†ⱼ} = 0
for i in 1:N, j in i:N
    push!(eq_cons_manual, ComplexF64(1.0) * (c_dag[i]*c_dag[j] + c_dag[j]*c_dag[i]))
end

# Canonical anti-commutation: {cᵢ, c†ⱼ} = δᵢⱼ
for i in 1:N, j in 1:N
    if i == j
        push!(eq_cons_manual, ComplexF64(1.0) * (c[i]*c_dag[i] + c_dag[i]*c[i]) - 1)
    else
        push!(eq_cons_manual, ComplexF64(1.0) * (c[i]*c_dag[j] + c_dag[j]*c[i]))
    end
end

# Nilpotent constraints: cᵢ² = 0, (c†ᵢ)² = 0
for i in 1:N
    push!(eq_cons_manual, ComplexF64(1.0) * c[i]*c[i])
    push!(eq_cons_manual, ComplexF64(1.0) * c_dag[i]*c_dag[i])
end

# **Step 4**: Create optimization problem with all constraints
pop_old = cpolyopt(ham_manual;
    eq_constraints=eq_cons_manual,
    comm_gps=[vcat(collect(c), collect(c_dag))],  # All ops in one group
    is_unipotent=false,
    is_projective=false
)

# This approach requires **14 constraint equations** for just 2 modes! The number
# of constraints scales as O(N²), making this approach impractical for larger systems.

# ## The Solution: Algebra Constructors

# The [`fermionic_algebra`](@ref) constructor encapsulates all these constraints
# automatically, allowing you to focus on the physics rather than the algebraic details.

# ### Simplified Approach with `fermionic_algebra`

# The same problem can be solved much more concisely:

sys = fermionic_algebra(N)
c_new, c_dag_new = sys.variables

# Construct the Hamiltonian (symmetrized for SDP)
ham_new = sum(ComplexF64(ε[i]/2) * (c_dag_new[i]*c_new[i] + c_new[i]*c_dag_new[i])
              for i in 1:N)

# Create the optimization problem - constraints are automatic!
pop_new = cpolyopt(ham_new, sys)

# Both approaches produce identical optimization problems, but the new interface is
# **dramatically shorter** and eliminates the possibility of errors in constraint specification.

# ## Computing Ground State Energy

# Let's solve for the ground state energy of the free fermion system.
# The expected ground state (vacuum) has energy 0, as both modes are empty.

solver_config = SolverConfig(
    optimizer=Mosek.Optimizer,  # SDP solver backend
    order=1                     # Relaxation order (1 is sufficient for free fermions)
)

res = cs_nctssos(pop_new, solver_config)

# The result should be approximately 0, corresponding to the vacuum state.
# For this free fermion system with negative single-particle energies,
# the ground state is the vacuum (no particles).
println("Ground state energy: ", res.objective)

# ## Adding Physical Constraints

# One key advantage of the algebra interface is the ease of adding custom constraints.
# For example, let's find the ground state energy with exactly 1 particle:

# Total particle number operator (symmetrized)
n1 = ComplexF64(0.5) * (c_dag_new[1]*c_new[1] + c_new[1]*c_dag_new[1])
n2 = ComplexF64(0.5) * (c_dag_new[2]*c_new[2] + c_new[2]*c_dag_new[2])
N_total = n1 + n2

# Constraint: N_total = 1
particle_constraint = N_total - 1

# Create problem with additional constraint
pop_constrained = cpolyopt(ham_new, sys; eq_constraints=[particle_constraint])

# Solve with particle number constraint
solver_config_higher = SolverConfig(
    optimizer=Mosek.Optimizer,
    order=2  # Higher order for particle number constraint
)
res_constrained = cs_nctssos(pop_constrained, solver_config_higher)

# The result should be approximately ε₁ = -1.0 (particle in lowest energy mode).
# The single particle occupies the lowest energy state.
println("Constrained ground state energy: ", res_constrained.objective)
println("Expected energy (ε₁): ", ε[1])

# ## Constraint Count Scaling

# The fermionic algebra constructor automatically generates constraints that scale
# as O(N²) for N modes:
#
# | N modes | Constraints | Formula     |
# |---------|-------------|-------------|
# | 1       | 5           | 2N² + 3N    |
# | 2       | 14          | 2N² + 3N    |
# | 3       | 27          | 2N² + 3N    |
# | 5       | 65          | 2N² + 3N    |
# | 10      | 230         | 2N² + 3N    |
#
# While this scaling means larger systems have many constraints, the algebra
# constructor ensures they are all correct, and modern SDP solvers handle
# O(100) constraints efficiently.

# ## Advantages of Fermionic Algebra Constructor

# The [`fermionic_algebra`](@ref) interface provides several key benefits:

# 1. **Automatic constraint generation**: All fermionic anti-commutation relations
#    and nilpotency constraints are encoded correctly.
#
# 2. **Error prevention**: Eliminates typos and sign errors in the 2N² + 3N constraints.
#
# 3. **Code clarity**: The physics intent is immediately clear from `fermionic_algebra(N)`.
#
# 4. **Scalability**: Works seamlessly for any number of modes without code modification.
#
# 5. **Consistency**: Ensures the same algebraic structure across different problems.
#
# 6. **Constraint-based approach**: All properties enforced via polynomial equality
#    constraints, making the solver's job transparent.
#
# 7. **Easy to extend**: Custom physical constraints (particle number, symmetries)
#    can be added simply via the `eq_constraints` parameter.

# ## Physical Applications

# The fermionic algebra interface is well-suited for studying:
#
# - **Free fermion systems**: Non-interacting electrons with arbitrary single-particle
#   Hamiltonians
#
# - **Fermi-Hubbard model**: Interacting electrons on lattices, central to
#   condensed matter physics
#
# - **Quantum chemistry**: Electronic structure problems with fermionic operators
#
# - **Superconductivity**: BCS theory and pairing phenomena
#
# - **Fermionic tensor networks**: Optimization problems in fermionic DMRG and PEPS

# ## Implementation Notes

# ### Symmetrization for SDP Solvers

# Note that fermionic number operators like $c^\dagger_i c_i$ are not symmetric
# polynomials (in the algebraic sense). For use as objectives in `cpolyopt`,
# we symmetrize them as:
#
# ```math
# n_i^{\text{symm}} = \frac{1}{2}(c^\dagger_i c_i + c_i c^\dagger_i) = \frac{1}{2}(\{c_i, c^\dagger_i\} + [c^\dagger_i, c_i])
# ```
#
# This symmetrized form equals the standard number operator plus a constant,
# and is compatible with SDP solvers while preserving the physics.

# ### Constraint Encoding

# The [`fermionic_algebra`](@ref) constructor uses a pure constraint-based approach:
#
# - **No modifications** to `FastPolynomials`, `PolyOpt`, or `ComplexPolyOpt`
# - All fermionic properties encoded as polynomial equality constraints
# - Nilpotent constraints ($c_i^2 = 0$) are technically redundant (implied by
#   anti-commutation) but included explicitly for numerical stability
# - SDP solver enforces all constraints during optimization

# ## Next Steps

# This interface extends naturally to more complex fermionic systems:
#
# - For interacting fermions, construct Hubbard-type Hamiltonians with the
#   fermionic operators
# - For systems with spin, use separate fermionic modes for spin-up and spin-down
# - For quantum chemistry, map molecular orbitals to fermionic modes
#
# The algebra constructor approach makes fermionic quantum many-body optimization
# more accessible and reliable, bridging physical intuition with rigorous
# mathematical formalism.
