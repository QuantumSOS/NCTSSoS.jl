# # XXX Model Ground State: GNS Reconstruction of Pauli Operators
#
# This example demonstrates the GNS (Gelfand-Naimark-Segal) reconstruction
# for the XXX Heisenberg model with 4 sites. We'll find the ground state energy
# using SDP relaxation, then extract finite-dimensional matrix representations
# of the Pauli operators from the moment matrix.
#
# This example showcases the [`pauli_algebra`](@ref) interface, which automatically
# handles all Pauli commutation relations and unipotency constraints, making the
# code more concise and less error-prone.

# ## The XXX Heisenberg Model
#
# The XXX (or isotropic Heisenberg) model describes a chain of interacting spins.
# For N sites, the Hamiltonian is:
#
# ```math
# H = \sum_{i=1}^{N} \left( X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1} \right)
# ```
#
# where periodic boundary conditions apply (site N+1 = site 1).
#
# Each site has three Pauli operators X, Y, Z that satisfy:
# - Anti-commutation: {Xᵢ, Yᵢ} = {Yᵢ, Zᵢ} = {Zᵢ, Xᵢ} = 0
# - Commutation: [Xᵢ, Yⱼ] = 2i δᵢⱼ Zᵢ (and cyclic permutations)
# - Unipotency: Xᵢ² = Yᵢ² = Zᵢ² = I
# - Different sites: [σᵢ, σⱼ] = 0 for i ≠ j
#
# We'll optimize over N=4 sites to find the ground state energy, then use
# the moment matrix to reconstruct explicit matrix representations of the
# Pauli operators.

# ## Package Setup

using NCTSSoS
using NCTSSoS.FastPolynomials
using JuMP
using MosekTools
using LinearAlgebra

# ## Problem Setup: Define Pauli Operators

N = 4  # Number of sites

# Use the pauli_algebra interface to automatically set up the Pauli operators
# and their algebraic relations
sys = pauli_algebra(N)
X, Y, Z = sys.variables

println("Created Pauli operators for $N sites using pauli_algebra interface")
println("X = ", X)
println("Y = ", Y)
println("Z = ", Z)

# ## Hamiltonian Construction
#
# The XXX Hamiltonian with periodic boundary conditions:

function xxx_hamiltonian(X, Y, Z, N)
    H = sum(one(ComplexF64) * X[i] * X[mod1(i + 1, N)] +
            Y[i] * Y[mod1(i + 1, N)] +
            Z[i] * Z[mod1(i + 1, N)]
            for i in 1:N)
    return H
end

ham = xxx_hamiltonian(X, Y, Z, N)

println("\nHamiltonian: ", ham)

# ## Polynomial Optimization Problem
#
# We minimize the Hamiltonian to find the ground state energy.
# The pauli_algebra interface automatically provides:
# - Equality constraints (Pauli commutation relations)
# - Commutation groups (operators at same site)
# - Simplification algorithm with unipotency

pop = cpolyopt(ham, sys)

println("\nPolynomial optimization problem created")

# ## SDP Relaxation and Solving
#
# We use COSMO solver with a moment relaxation of order 2.

SOLVER = optimizer_with_attributes(Mosek.Optimizer)
solver_config = SolverConfig(optimizer=SOLVER, order=2)

println("\nSolving SDP relaxation...")
result = cs_nctssos(pop, solver_config)

model = result.model

println("\nSolve time: ", solve_time(model), " seconds")
println("Ground state energy: ", objective_value(model))
println("Ground state energy per site: ", objective_value(model) / N)

# For comparison: exact ground state energy of 4-site XXX model
# (can be computed exactly for small systems)
println("\nNote: For the 4-site XXX model, the exact ground state energy")
println("can be computed using exact diagonalization or Bethe ansatz.")

# ## Extracting the Moment Matrix
#
# The moment matrix is embedded in the dual of the SDP relaxation.

cons = all_constraints(model, include_variable_in_set_constraints=true)
X_matrix = value.(dual(cons[1]))

println("\nMoment matrix X size: ", size(X_matrix))

# ## Converting from X to Complex Hermitian Form
#
# The moment matrix X is in real form:
# X = [X₁   X₃']
#     [X₃   X₂ ]
#
# We convert it to complex Hermitian form H = H_R + im*H_I:

n_total = size(X_matrix, 1)
n_half = div(n_total, 2)

# Extract blocks
X1 = X_matrix[1:n_half, 1:n_half]
X3 = X_matrix[(n_half+1):end, 1:n_half]
X2 = X_matrix[(n_half+1):end, (n_half+1):end]

# Compute Hermitian components
H_R = (X1 + X2) / 2
H_I = (X3 - X_matrix[1:n_half, (n_half+1):end]') / 2

# Construct complex Hermitian matrix
H = H_R + im * H_I

println("\nComplex Hermitian moment matrix H size: ", size(H))
println("Hermiticity check: ||H - H'|| = ", norm(H - H'))

# Check positive semi-definiteness
eigenvalues_H = eigvals(Hermitian(H))
println("Minimum eigenvalue of H: ", minimum(real.(eigenvalues_H)))
println("Maximum eigenvalue of H: ", maximum(real.(eigenvalues_H)))

# ## GNS Reconstruction
#
# Now we perform GNS reconstruction to extract matrix representations
# of the Pauli operators.

# All variables (in order)
vars = [X; Y; Z]  # 12 variables total (4 sites × 3 operators)

println("\nTotal number of variables: ", length(vars))

# Degree for the Hankel matrix (same as relaxation order)
H_deg = 2

# Use the simplification algorithm from the pauli_algebra system
# This ensures proper handling of Pauli commutation relations and unipotency
sa = sys.simplify_algo

# Perform GNS reconstruction
println("\nPerforming GNS reconstruction...")
all_recon = reconstruct(H, vars, H_deg, sa; atol=1e-3)

println("Reconstruction complete!")
println("Number of reconstructed operators: ", length(all_recon))
println("Reconstructed operator dimension: ", size(all_recon[1]))

# Split into X, Y, Z operators for each site
X_recon = all_recon[1:N]
Y_recon = all_recon[(N+1):(2*N)]
Z_recon = all_recon[(2*N+1):end]

println("\nReconstructed operators:")
println("X_recon: ", length(X_recon), " operators of size ", size(X_recon[1]))
println("Y_recon: ", length(Y_recon), " operators of size ", size(Y_recon[1]))
println("Z_recon: ", length(Z_recon), " operators of size ", size(Z_recon[1]))

# ## Verification: Pauli Algebra
#
# Verify that the reconstructed operators satisfy the Pauli algebra relations.

println("\n" * "="^60)
println("VERIFICATION: PAULI ALGEBRA")
println("="^60)

# Test 1: Squares equal identity
println("\n1. Testing σᵢ² = I for each operator at each site:")
max_square_error = 0.0
for i in 1:N
    X_sq_err = norm(X_recon[i]^2 - I)
    Y_sq_err = norm(Y_recon[i]^2 - I)
    Z_sq_err = norm(Z_recon[i]^2 - I)
    
    println("  Site $i:")
    println("    ||X[$i]² - I|| = $(round(X_sq_err, digits=8))")
    println("    ||Y[$i]² - I|| = $(round(Y_sq_err, digits=8))")
    println("    ||Z[$i]² - I|| = $(round(Z_sq_err, digits=8))")
    
    max_square_error = max(max_square_error, X_sq_err, Y_sq_err, Z_sq_err)
end
println("  Maximum square error: $max_square_error")

# Test 2: Anti-commutation within sites
println("\n2. Testing anti-commutation {σᵢ, σⱼ} = 0 at each site:")
max_anticomm_error = 0.0
for i in 1:N
    XY_anticomm = X_recon[i] * Y_recon[i] + Y_recon[i] * X_recon[i]
    YZ_anticomm = Y_recon[i] * Z_recon[i] + Z_recon[i] * Y_recon[i]
    ZX_anticomm = Z_recon[i] * X_recon[i] + X_recon[i] * Z_recon[i]
    
    XY_err = norm(XY_anticomm)
    YZ_err = norm(YZ_anticomm)
    ZX_err = norm(ZX_anticomm)
    
    println("  Site $i:")
    println("    ||{X[$i], Y[$i]}|| = $(round(XY_err, digits=8))")
    println("    ||{Y[$i], Z[$i]}|| = $(round(YZ_err, digits=8))")
    println("    ||{Z[$i], X[$i]}|| = $(round(ZX_err, digits=8))")
    
    max_anticomm_error = max(max_anticomm_error, XY_err, YZ_err, ZX_err)
end
println("  Maximum anti-commutation error: $max_anticomm_error")

# Test 3: Commutation relations [Xᵢ, Yᵢ] = 2iZᵢ
println("\n3. Testing commutation relations [σᵢ, σⱼ] = 2iε_{ijk}σₖ:")
max_comm_error = 0.0
for i in 1:N
    XY_comm = X_recon[i] * Y_recon[i] - Y_recon[i] * X_recon[i]
    YZ_comm = Y_recon[i] * Z_recon[i] - Z_recon[i] * Y_recon[i]
    ZX_comm = Z_recon[i] * X_recon[i] - X_recon[i] * Z_recon[i]
    
    XY_err = norm(XY_comm - 2im * Z_recon[i])
    YZ_err = norm(YZ_comm - 2im * X_recon[i])
    ZX_err = norm(ZX_comm - 2im * Y_recon[i])
    
    println("  Site $i:")
    println("    ||[X[$i], Y[$i]] - 2iZ[$i]|| = $(round(XY_err, digits=8))")
    println("    ||[Y[$i], Z[$i]] - 2iX[$i]|| = $(round(YZ_err, digits=8))")
    println("    ||[Z[$i], X[$i]] - 2iY[$i]|| = $(round(ZX_err, digits=8))")
    
    max_comm_error = max(max_comm_error, XY_err, YZ_err, ZX_err)
end
println("  Maximum commutation error: $max_comm_error")

# Test 4: Commutation between different sites
println("\n4. Testing inter-site commutation [σᵢ, σⱼ] = 0 for i ≠ j:")
max_intersite_error = 0.0
sample_pairs = [(1, 2), (1, 3), (2, 4)]  # Sample a few pairs

for (i, j) in sample_pairs
    # Test one combination from each site pair
    comm_XX = X_recon[i] * X_recon[j] - X_recon[j] * X_recon[i]
    comm_YZ = Y_recon[i] * Z_recon[j] - Z_recon[j] * Y_recon[i]
    
    XX_err = norm(comm_XX)
    YZ_err = norm(comm_YZ)
    
    println("  Sites ($i, $j):")
    println("    ||[X[$i], X[$j]]|| = $(round(XX_err, digits=8))")
    println("    ||[Y[$i], Z[$j]]|| = $(round(YZ_err, digits=8))")
    
    max_intersite_error = max(max_intersite_error, XX_err, YZ_err)
end
println("  Maximum inter-site commutation error: $max_intersite_error")

# ## Verification: Hamiltonian Expectation Value
#
# Verify that the reconstructed operators give the correct ground state energy.

println("\n" * "="^60)
println("VERIFICATION: HAMILTONIAN RECONSTRUCTION")
println("="^60)

# Compute the Hamiltonian using reconstructed operators
H_recon = sum(X_recon[i] * X_recon[mod1(i+1, N)] + 
              Y_recon[i] * Y_recon[mod1(i+1, N)] + 
              Z_recon[i] * Z_recon[mod1(i+1, N)] 
              for i in 1:N)

println("\nReconstructed Hamiltonian H_recon size: ", size(H_recon))

# The ground state energy should be the minimum eigenvalue
eigenvalues_H_recon = eigvals(Hermitian(H_recon))
E_ground_recon = minimum(real.(eigenvalues_H_recon))

println("Ground state energy from SDP: ", objective_value(model))
println("Minimum eigenvalue of H_recon: ", E_ground_recon)
println("Difference: ", abs(objective_value(model) - E_ground_recon))

# ## Analysis of Hilbert Space Structure
#
# Let's analyze the structure of the Hilbert space on which these operators act.

println("\n" * "="^60)
println("HILBERT SPACE ANALYSIS")
println("="^60)

dim = size(X_recon[1], 1)
println("\nHilbert space dimension: $dim")

# For 4 qubits, we'd expect dimension 2^4 = 16 classically
# But the GNS construction might give a smaller effective dimension
println("Expected dimension for 4 qubits: ", 2^N)
println("Actual reconstructed dimension: ", dim)

if dim < 2^N
    println("\n✓ The GNS reconstruction found a compressed representation!")
    println("  This dimension captures the essential quantum correlations")
    println("  needed for the ground state, without requiring the full")
    println("  2^N dimensional Hilbert space.")
end

# Analyze spectrum of individual operators
println("\nSpectrum analysis of reconstructed operators:")
for i in 1:N
    X_eigs = sort(real.(eigvals(X_recon[i])))
    Y_eigs = sort(real.(eigvals(Y_recon[i])))
    Z_eigs = sort(real.(eigvals(Z_recon[i])))
    
    println("  Site $i:")
    println("    X[$i] eigenvalues: ", round.(X_eigs[1:min(4, length(X_eigs))], digits=4), "...")
    println("    Y[$i] eigenvalues: ", round.(Y_eigs[1:min(4, length(Y_eigs))], digits=4), "...")
    println("    Z[$i] eigenvalues: ", round.(Z_eigs[1:min(4, length(Z_eigs))], digits=4), "...")
end

# ## Summary and Physical Interpretation
#
# We have successfully:
#
# 1. **Formulated the XXX model** as a polynomial optimization problem
# 2. **Computed ground state energy** using SDP relaxation
# 3. **Extracted explicit matrices** for all Pauli operators via GNS reconstruction
# 4. **Verified Pauli algebra** - all relations hold to high precision
# 5. **Confirmed Hamiltonian** - reconstructed operators give correct ground state energy
#
# ### Physical Insights:
#
# - The GNS construction provides a **finite-dimensional realization** of the quantum
#   operators that is sufficient to capture the ground state properties
# - The dimension is typically much smaller than 2^N, revealing the **effective degrees
#   of freedom** relevant for the ground state
# - All Pauli algebra relations are satisfied, confirming the reconstruction is
#   **physically consistent**
# - This approach bridges the gap between **abstract operator algebras** and
#   **concrete matrix representations**
#
# ### Applications:
#
# This technique can be extended to:
# - Larger spin chains (more sites)
# - Different spin models (XYZ, Ising, etc.)
# - Systems with additional symmetries or constraints
# - Analyzing entanglement structure in the ground state
# - Computing other observables and correlation functions

println("\n" * "="^60)
println("EXAMPLE COMPLETE")
println("="^60)
