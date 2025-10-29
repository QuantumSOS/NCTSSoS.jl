# # Simplified Quantum Spin Models with Pauli Algebra Interface

# NCTSSoS provides a convenient interface for working with common quantum algebras,
# eliminating the need to manually specify commutation relations and constraints.
# This tutorial demonstrates the `pauli_algebra` constructor for quantum spin systems.

# ## Basic Usage: Heisenberg XXX Model

# The traditional approach requires manually defining Pauli commutation relations:

using NCTSSoS, MosekTools
N = 6
@ncpolyvar x[1:N] y[1:N] z[1:N]

# Construct the XXX Heisenberg Hamiltonian
ham = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)] for op in [x, y, z] for i in 1:N)

# **Old way**: Manual constraint specification
eq_cons = reduce(vcat, [
    [x[i] * y[i] - im * z[i],
     y[i] * x[i] + im * z[i],
     y[i] * z[i] - im * x[i],
     z[i] * y[i] + im * x[i],
     z[i] * x[i] - im * y[i],
     x[i] * z[i] + im * y[i]]
    for i in 1:N
])

pop_old = cpolyopt(ham;
    eq_constraints=eq_cons,
    comm_gps=[[x[i], y[i], z[i]] for i in 1:N],
    is_unipotent=true
)

# **New way**: Using the algebra constructor
sys = pauli_algebra(N)
x_new, y_new, z_new = sys.variables

# Construct the same Hamiltonian with new variables
ham_new = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)]
              for op in [x_new, y_new, z_new] for i in 1:N)

pop_new = cpolyopt(ham_new, sys)

# Both approaches produce identical optimization problems, but the new interface is
# more concise and less error-prone.

# ## Solving the Optimization Problem

solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=2)
res = cs_nctssos(pop_new, solver_config)
energy_per_site = res.objective / N

# The result matches the known ground state energy per site of approximately -0.467129
# for the 6-site XXX Heisenberg model.

# ## Combining Algebra with Custom Constraints

# ## Summary

# The `pauli_algebra` interface provides:
# 1. **Automatic constraint generation**: No need to manually write commutation relations
# 2. **Error prevention**: Correct algebra structure guaranteed
# 3. **Code clarity**: Physics intent is clear from `pauli_algebra(N)`
# 4. **Flexibility**: Can still add custom constraints when needed
# 5. **Consistency**: Same interface works for different system sizes
