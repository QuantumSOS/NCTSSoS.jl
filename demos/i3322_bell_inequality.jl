# Demo: I₃₃₂₂ Bell inequality — maximum quantum violation
#
# Persona: knows polynomial optimization (has used ncpol2sdpa), new to NCTSSoS.
# Tests: API translation, UnipotentAlgebra grouping, sparsity comparison.
#
# Problem:
#   Maximize I = A₁(B₁+B₂) + A₂(B₁+B₂) + A₃(B₁-B₂) - A₁ - 2B₁ - B₂
#   over dichotomic observables Aᵢ²=Bⱼ²=I, [Aᵢ,Bⱼ]=0.
#
# Reference: the classical bound is 0, the quantum bound is ≈ 0.25096.

using NCTSSoS
using Clarabel

# --- Variables ---
# Separate groups: Alice (3 measurements) and Bob (2 measurements)
# → different-group variables commute, encoding [Aᵢ,Bⱼ]=0
registry, (A, B) = create_unipotent_variables([("A", 1:3), ("B", 1:2)])

# --- Bell functional ---
I_bell = A[1]*(B[1]+B[2]) + A[2]*(B[1]+B[2]) + A[3]*(B[1]-B[2]) - A[1] - 2*B[1] - B[2]

# Minimize -I to get the maximum quantum violation
pop = polyopt(-I_bell, registry)

# --- Solve: dense (no sparsity) at NPA level 2 ---
config_dense = SolverConfig(optimizer=Clarabel.Optimizer, order=2)
result_dense = cs_nctssos(pop, config_dense)
bound_dense = -result_dense.objective
println("Dense bound (order=2):  ", bound_dense)

# --- Solve: with correlative + term sparsity ---
config_sparse = SolverConfig(optimizer=Clarabel.Optimizer, order=2,
                             cs_algo=MF(), ts_algo=MMD())
result_sparse = cs_nctssos(pop, config_sparse)
bound_sparse = -result_sparse.objective
println("Sparse bound (order=2): ", bound_sparse)

# --- Compare ---
gap = abs(bound_dense - bound_sparse)
println("Gap (dense - sparse):   ", gap)
println("Expected quantum bound: ≈ 0.25096")
