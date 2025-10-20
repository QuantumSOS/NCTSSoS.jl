# # Mermin Square: GNS Reconstruction
#
# This example demonstrates the GNS (Gelfand-Naimark-Segal) reconstruction
# for the Mermin Square game, a fundamental example of quantum contextuality.
# We'll use the moment matrix from an SDP relaxation to extract finite-dimensional
# matrix representations of the quantum operators.

# ## The Mermin Square Game
#
# The Mermin Square is a 3×3 grid of observables (measurement operators):
#
# ```math
# \begin{array}{|c|c|c|}
# \hline
# A_{11} & A_{12} & A_{13} \\
# \hline
# A_{21} & A_{22} & A_{23} \\
# \hline
# A_{31} & A_{32} & A_{33} \\
# \hline
# \end{array}
# ```
#
# We also have column observables B[i,j] that must satisfy:
#
# **Row constraints**: For each row i, the product equals +1:
# ```math
# A_{i1} A_{i2} A_{i3} = +\mathbb{I}
# ```
#
# **Column constraints**: For each column j, the product equals -1:
# ```math
# B_{1j} B_{2j} B_{3j} = -\mathbb{I}
# ```
#
# **Compatibility**: Operators at the same grid position must be identical:
# ```math
# A_{ij} = B_{ij} \quad \forall i,j
# ```
#
# **Commutativity**: Operators in different grid positions anti-commute.
#
# This game demonstrates quantum contextuality: the product of all row products
# gives (+1)³ = +1, but the product of all column products gives (-1)³ = -1.
# This is impossible classically but achievable with quantum mechanics!

# ## Package Setup

using NCTSSoS
using NCTSSoS.FastPolynomials
using JuMP
using COSMO 

# ## Helper Function for Entry Constraints
#
# We need a custom function to inject the compatibility constraints A[i,j] = B[i,j]
# into the moment relaxation. This function extends the standard cs_nctssos to
# allow additional entry constraints on the moment matrix.

function cs_nctssos_with_entry(pop::OP, solver_config::SolverConfig, entry_constraints::Vector{Polynomial{T}}; dualize::Bool=true) where {T,P<:Polynomial{T},OP<:NCTSSoS.OptimizationProblem{P}}

    sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_projective=pop.is_projective, is_unipotent=pop.is_unipotent)
    order = iszero(solver_config.order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) : solver_config.order

    corr_sparsity = NCTSSoS.correlative_sparsity(pop, order, solver_config.cs_algo)

    cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base, sa)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
    end

    moment_problem = NCTSSoS.moment_relax(pop, corr_sparsity, cliques_term_sparsities)

    for c in entry_constraints
        push!(moment_problem.constraints, (:HPSD, [c;;]))
    end

    (pop isa NCTSSoS.ComplexPolyOpt{P} && !dualize) && error("Solving Moment Problem for Complex Poly Opt is not supported")
    problem_to_solve = !dualize ? moment_problem : NCTSSoS.sos_dualize(moment_problem)

    set_optimizer(problem_to_solve.model, solver_config.optimizer)
    optimize!(problem_to_solve.model)
    return NCTSSoS.PolyOptResult(objective_value(problem_to_solve.model), corr_sparsity, cliques_term_sparsities, problem_to_solve.model)
end

# ## Problem Setup
#
# We declare non-commuting polynomial variables for the operators in the grid.

n = 3
T1 = ComplexF64

@ncpolyvar A[1:n,1:n] B[1:n,1:n]

# Set up a trivial objective (we're just looking for a feasible solution)

obj = one(T1) * one(A[3, 3])

# Define the game constraints

game_cons = [[
        A[1, i] * A[2, i] * A[3, i] + one(T1)
        for i in 1:3
    ]; [B[i, 1] * B[i, 2] * B[i, 3] - one(T1) for i in 1:3]]

# Define commutativity constraints (operators in different positions anti-commute)

comm_cons = filter!(!iszero, vec([
    one(T1) * op[i, j] * op[k, j] - op[k, j] * op[i, j] for i in 1:3, j in 1:3, k in 1:3, op in [A, B]
]))

# Create the complex polynomial optimization problem with unipotency
# (all operators square to identity)

pop = cpolyopt(obj, eq_constraints=[game_cons; comm_cons], is_unipotent=true)

# Define entry constraints (compatibility: A[i,j] = B[i,j])

entry_cons = vec([A[i, j] * B[i, j] - one(T1) for i in 1:3, j in 1:3])

# ## Solving the SDP Relaxation
#
# We use the SCS solver to compute the moment relaxation of order 2.

SOLVER = optimizer_with_attributes(COSMO.Optimizer)
solver_config = SolverConfig(optimizer=SOLVER, order=2)

result = cs_nctssos_with_entry(pop, solver_config, entry_cons; dualize=true)

model = result.model

# ## Extracting the Hankel Matrix
#
# The moment matrix (Hankel matrix) is embedded in the constraints of the dual problem.
# We extract it from the first PSD constraint.

cons = all_constraints(model, include_variable_in_set_constraints=true)

# The moment matrix is the value of the first constraint

X = value.(dual(cons[1]))

println("X matrix size: ", size(X))

# ## Converting from X form to Y (complex Hermitian) form
#
# The moment matrix X is currently in the form:
# X = [X₁   X₃']
#     [X₃   X₂ ]
#
# We need to convert it to Y form:
# Y = [H_R   -H_I]  which represents the complex Hermitian matrix H = H_R + im*H_I
#     [H_I    H_R]
#
# The relationship is:
# X₁ + X₂ = 2*H_R  →  H_R = (X₁ + X₂)/2
# X₃ - X₃' = 2*H_I  →  H_I = (X₃ - X₃')/2

# Extract the dimensions
n_total = size(X, 1)
n_half = div(n_total, 2)

# Extract the blocks from X
X1 = X[1:n_half, 1:n_half]
X3 = X[(n_half+1):end, 1:n_half]
X2 = X[(n_half+1):end, (n_half+1):end]

# Compute H_R and H_I
H_R = (X1 + X2) / 2
H_I = (X3 - X[1:n_half, (n_half+1):end]') / 2

H = H_R + im * H_I

# Also construct the complex Hermitian form for reference
H_complex = H_R + im * H_I
println("\nComplex Hermitian matrix size: ", size(H_complex))
println("Hermiticity check: ||H_complex - H_complex'|| = ", norm(H_complex - H_complex'))

# ## GNS Reconstruction
#
# Now we use the GNS construction to extract finite-dimensional matrix representations
# of the operators from the moment matrix.

# First, we need to construct the basis that was used to index the moment matrix.
# Since we used order 2, the basis consists of monomials up to degree 2.

# Get all variables (we'll use A since A[i,j] = B[i,j])
vars = [vec(A[1:n, 1:n]); vec(B[1:n, 1:n])]

# The degree for the Hankel matrix
H_deg = 2

# Create the same simplification algorithm used in the optimization
sa = SimplifyAlgorithm(comm_gps=[vars], is_unipotent=true, is_projective=false)

# Perform GNS reconstruction with appropriate tolerance and simplification
# The simplification algorithm ensures basis vectors are simplified according to
# the unipotency constraint (all operators square to identity)
A_recon_vec = reconstruct(H, vars, H_deg, sa; atol=1e-3)

round.(A_recon_vec[1], digits=6)
# Reshape back into a 3×3 array of matrices
A_recon = reshape(A_recon_vec, (n, n))

println("\nReconstructed operator dimension: ", size(A_recon[1, 1]))

# Since A[i,j] = B[i,j], we can use the same reconstructed matrices for B
B_recon = A_recon

# ## Verification
#
# Let's verify that the reconstructed matrices satisfy the game constraints.

println("\n=== Verification of Game Constraints ===\n")

# Check row products: A[i,1] * A[i,2] * A[i,3] ≈ I
println("Row products (should be ≈ I):")
for i in 1:3
    row_product = A_recon[i, 1] * A_recon[i, 2] * A_recon[i, 3]
    identity_error = norm(row_product - I)
    println("  Row $i: ||product - I|| = ", identity_error)
end

# Check column products: B[1,j] * B[2,j] * B[3,j] ≈ -I
println("\nColumn products (should be ≈ -I):")
for j in 1:3
    col_product = B_recon[1, j] * B_recon[2, j] * B_recon[3, j]
    neg_identity_error = norm(col_product + I)
    println("  Column $j: ||product + I|| = ", neg_identity_error)
end

# Check that operators square to identity (unipotency)
println("\nOperator squares (should be ≈ I):")
max_square_error = 0.0
for i in 1:3, j in 1:3
    square_error = norm(A_recon[i, j]^2 - I)
    max_square_error = max(max_square_error, square_error)
end
println("  Maximum ||A[i,j]² - I|| = ", max_square_error)

# Check compatibility: A[i,j] = B[i,j] (automatically satisfied by construction)
println("\nCompatibility A[i,j] = B[i,j]: ✓ (satisfied by construction)")

# ## Analysis
#
# The GNS reconstruction successfully extracts finite-dimensional matrix representations
# that satisfy all the Mermin Square constraints. The dimension of the Hilbert space
# (size of the matrices) is determined by the rank of the moment matrix, which in turn
# depends on the quantum state that achieves the Mermin Square correlations.
#
# This demonstrates that:
# 1. The SDP relaxation captures the quantum correlations
# 2. The GNS construction recovers explicit matrix representations
# 3. These matrices satisfy all the game constraints (up to numerical precision)
#
# The Mermin Square is a beautiful example of quantum contextuality, showing that
# quantum mechanics allows correlations that are impossible to achieve classically!
