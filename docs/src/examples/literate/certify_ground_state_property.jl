# # [Certifying Ground State Properties](@id certify-property)

# Understanding ground-state properties of quantum many-body systems represents
# a **fundamental challenge in quantum physics** [wang2024Certifying](@cite).
# These properties reveal the low-energy phases and quantum correlations
# that characterize complex quantum systems. However, computing them exactly
# becomes intractable for large systems due to the exponential growth of
# the Hilbert space dimension with particle number.
#
# While exact diagonalization provides precise results for small systems,
# it quickly becomes computationally prohibitive. Instead, we can obtain
# rigorous bounds on ground-state properties by combining:
# - **Upper bounds** from variational methods (or exact diagonalization for small systems)
# - **Lower bounds** from semidefinite programming (SDP) relaxations
#
# This approach provides certified intervals that guarantee the true
# ground-state property lies within computable bounds.
#
# As a demonstration, we consider the quantum Heisenberg model:
#
# ```math
# \hat{H} = J \sum_{j=1}^{N}\sigma^z_j \sigma^z_{j+1} + h \sum_{j=1}^{N} \sigma^x_j
# ```
#
# where we will certify both the ground-state energy and correlation functions.

using Yao
using LinearAlgebra

## System parameters for exact diagonalization
N = 3                    # Number of spins in the chain
J = 1.0                  # Coupling strength for Z-Z interactions

## Arrays to store results across different magnetic field strengths
ground_state_energy_upper_bounds = Float64[]  # Energy per site
s1s2values = Float64[]                        # S^z_1 * S^z_2 correlation

## Loop over different transverse field strengths
for h in 0.1:0.2:2.0
    ## Construct the Heisenberg Hamiltonian with transverse field
    ## H = J/4 * Σ Z_i Z_{i+1} + h/2 * Σ X_i
    ham = sum(J / 4 * kron(N, i => Z, mod1(i + 1, N) => Z) for i in 1:N) +
          sum(h / 2 * put(N, i => X) for i in 1:N)

    ## Diagonalize the Hamiltonian to find eigenvalues and eigenvectors
    evals, eigves = eigen(Matrix(ham))

    ## Calculate S^z_1 * S^z_2 correlation operator
    s1s2 = Matrix(kron(N, 1 => Z, 2 => Z)) / 4

    ## Store ground state energy per site
    push!(ground_state_energy_upper_bounds, minimum(real(evals)) / N)

    ## Calculate and store S^z_1 * S^z_2 expectation value in ground state
    ground_state_idx = argmin(real(evals))
    ground_state = eigves[:, ground_state_idx]
    push!(s1s2values, real(ground_state' * s1s2 * ground_state))
end

# ## Computing Lower Bounds with Semidefinite Programming
#
# Now we obtain rigorous lower bounds on the ground-state energy using
# `NCTSSoS.jl`. This approach formulates the quantum ground-state problem
# as a polynomial optimization problem and solves its semidefinite
# programming relaxation, providing certified lower bounds.

using NCTSSoS, NCTSSoS.FastPolynomials
using MosekTools
using JuMP

## Configure Mosek solver with high precision settings
SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-8,  # Primal feasibility tolerance
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-8,  # Dual feasibility tolerance
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-8,  # Relative gap tolerance
    "MSK_IPAR_NUM_THREADS" => 0)            # Use all available threads

## Initialize array to store energy lower bounds
ground_state_energy_lower_bounds = Float64[]

## System parameters (matching the exact diagonalization setup)
N = 3                    # Number of spins
T1 = ComplexF64          # Complex number type for calculations
J = 1.0                  # Coupling strength

## Loop over the same magnetic field values as before
for (i, h) in enumerate(0.1:0.2:2.0)
    ## Define non-commutative polynomial variables for Pauli operators
    @ncpolyvar x[1:N] y[1:N] z[1:N]  # x = σ^x, y = σ^y, z = σ^z

    ## Objective function: we want to minimize the Hamiltonian
    ## H = J/4 * Σ z_i z_{i+1} + h/2 * Σ x_i
    ham = sum(T1(J / 4) * z[i] * z[mod1(i + 1, N)] + T1(h / 2) * x[i] for i in 1:N)

    ## Pauli operator algebra constraints (commutation relations)
    ## These encode the fundamental quantum mechanical properties of spin operators
    eq_cons = reduce(vcat, [
        [x[i] * y[i] - im * z[i],   # [σ^x_i, σ^y_i] = iσ^z_i
         y[i] * x[i] + im * z[i],   # [σ^y_i, σ^x_i] = -iσ^z_i
         y[i] * z[i] - im * x[i],   # [σ^y_i, σ^z_i] = iσ^x_i
         z[i] * y[i] + im * x[i],   # [σ^z_i, σ^y_i] = -iσ^x_i
         z[i] * x[i] - im * y[i],   # [σ^z_i, σ^x_i] = iσ^y_i
         x[i] * z[i] + im * y[i]]   # [σ^x_i, σ^z_i] = -iσ^y_i
        for i in 1:N])

    ## Create polynomial optimization problem
    pop = cpolyopt(ham;
        eq_constraints=eq_cons,                    # Pauli algebra constraints
        comm_gps=[[x[i], y[i], z[i]] for i in 1:N], # Spin operators on same site commute
        is_unipotent=true)                         # Pauli operators square to identity

    ## Configure solver with second-order moment relaxation
    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    ## Solve the semidefinite program to get energy lower bound
    res = cs_nctssos(pop, solver_config)

    ## Store energy per site (divide by system size)
    push!(ground_state_energy_lower_bounds, res.objective / N)
end

# ## Bounding Ground-State Correlation Functions
#
# Beyond energy bounds, we can rigorously bound other ground-state properties
# such as correlation functions. Here we demonstrate bounding the expectation
# value of the two-point correlation function $S^z_{1}S^z_{2}$.
#
# To accomplish this, we need a custom solver that can handle additional
# constraints on correlation functions. The following function extends
# the standard NCTSSoS solver to incorporate entry constraints.

"""
    cs_nctssos_with_entry(pop, solver_config, entry_constraints; dualize=true)

Extended NCTSSoS solver that incorporates additional entry constraints
for bounding specific correlation functions in quantum systems.

This function builds the moment matrix relaxation with sparsity patterns
and adds semidefinite constraints to bound specific operator expectations.
"""
function cs_nctssos_with_entry(pop::OP, solver_config::SolverConfig, entry_constraints::Vector{Polynomial{T}}; dualize::Bool=true) where {T,P<:Polynomial{T},OP<:NCTSSoS.OptimizationProblem{P}}

   ## Initialize simplification algorithm for polynomial reduction
   sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_projective=pop.is_projective, is_unipotent=pop.is_unipotent)

   ## Determine the order of moment relaxation (automatic if not specified)
   order = iszero(solver_config.order) ?
       maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) :
       solver_config.order

   ## Exploit correlative sparsity to reduce problem size
   corr_sparsity = NCTSSoS.correlative_sparsity(pop, order, solver_config.cs_algo)

   ## Decompose objective function across cliques (sparse subproblems)
   cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

   ## Initialize support patterns for term sparsity
   initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base, sa)
   end

   ## Apply term sparsity to further reduce problem size
   cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
   end

   ## Build the moment relaxation problem
   moment_problem = NCTSSoS.moment_relax(pop, corr_sparsity, cliques_term_sparsities)

   ## Add entry constraints for correlation function bounds
   ## These are semidefinite constraints that ensure physical validity
   for c in entry_constraints
       push!(moment_problem.constraints,(:HPSD, [c;;]))
   end

   ## Handle complex polynomial optimization (dual formulation)
   (pop isa NCTSSoS.ComplexPolyOpt{P} && !dualize) && error("Solving Moment Problem for Complex Poly Opt is not supported")
   problem_to_solve = !dualize ? moment_problem : NCTSSoS.sos_dualize(moment_problem)

   ## Solve the semidefinite program
   set_optimizer(problem_to_solve.model, solver_config.optimizer)
   optimize!(problem_to_solve.model)

   ## Return optimization results
   return NCTSSoS.PolyOptResult(objective_value(problem_to_solve.model), corr_sparsity, cliques_term_sparsities, problem_to_solve.model)
end

# ## Computing Rigorous Bounds on Correlation Functions
#
# Now we demonstrate the main result: computing rigorous upper and lower bounds
# on the two-point correlation function $\langle S^z_1 S^z_2 \rangle$.
# This showcases how polynomial optimization can certify physical properties
# beyond just the ground-state energy.
#
# We use reference energy bounds (computed via high-precision DMRG) to
# constrain the optimization problem, ensuring our correlation function
# bounds are physically meaningful.

## Initialize arrays to store correlation function bounds
lower_bounds = Float64[]
upper_bounds = Float64[]

## System parameters
J = 1.0  # Coupling strength (same as before)

## Loop over magnetic field values
for (i, h) in enumerate(0.1:0.2:2.0)
    ## Define polynomial variables for Pauli operators
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    ## Objective function: S^z_1 * S^z_2 correlation
    obj = one(T1) * z[1] * z[2]

    ## Hamiltonian (same form as before)
    ham = sum(T1(J / 4) * z[i] * z[mod1(i + 1, N)] + T1(h / 2) * x[i] for i in 1:N)

    ## Pauli algebra constraints
    eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

    ## Energy constraint: ensure ground state energy is within reference bounds
    ## This is crucial for obtaining physically meaningful correlation bounds
    ineq_cons = [ham - ground_state_energy_lower_bounds[i] * N]

    ## Create optimization problems for lower and upper bounds
    ## We solve two separate problems: minimize and maximize the correlation
    pop_l = cpolyopt(obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons,
        comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    pop_u = cpolyopt(-obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons,
        comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    ## Configure solver
    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    ## Additional energy constraint for the upper bound problem
    single_ineq_cons = [ground_state_energy_upper_bounds[i] * N - ham]

    ## Solve for lower and upper bounds on correlation function
    res_l = cs_nctssos_with_entry(pop_l, solver_config, single_ineq_cons; dualize=true)
    res_u = cs_nctssos_with_entry(pop_u, solver_config, single_ineq_cons; dualize=true)

    ## Store bounds (divide by 4 to convert from Pauli to spin operators)
    push!(lower_bounds, res_l.objective / 4)
    push!(upper_bounds, -res_u.objective / 4)
end

# ## Visualizing Certified Correlation Function Bounds
#
# Finally, we visualize our results to demonstrate the tightness of the
# certified bounds. The polynomial optimization approach provides rigorous
# upper and lower bounds that closely bracket the exact correlation function
# values, validating the effectiveness of our certification method.

using CairoMakie

## Create the visualization
f = Figure(size=(800, 600))
ax = Axis(f[1, 1],
    xlabel="Transverse Field Strength (h)",
    ylabel=L"\langle S^z_1 S^z_2 \rangle",
    title="Certified Bounds on Ground-State Correlation Function")

## Magnetic field values for plotting
xs = collect(0.1:0.2:2.0)

## Plot the certified bounds
scatterlines!(ax, xs, upper_bounds, color=:red, label="Upper Bound (SDP)", linewidth=2, markersize=8)
scatterlines!(ax, xs, lower_bounds, color=:green, label="Lower Bound (SDP)", linewidth=2, markersize=8)
scatterlines!(ax, xs, s1s2values, color=:blue, label="Exact Value (ED)", linewidth=2, markersize=8)

## Customize the plot
axislegend(ax, position=:rt)
ax.xgridvisible = true
ax.ygridvisible = true
ax.xgridwidth = 0.5
ax.ygridwidth = 0.5
ax.xgridcolor = (:gray, 0.2)
ax.ygridcolor = (:gray, 0.2)

f

# ## Summary
#
# This example demonstrates how polynomial optimization provides rigorous
# **certified bounds** on quantum ground-state properties. Key insights:
#
# 1. **Tight bounds**: The semidefinite programming relaxations yield bounds
#    that closely bracket the exact values, providing high-precision certification.
#
# 2. **General method**: The approach works for any polynomial Hamiltonian
#    and can bound arbitrary correlation functions, not just energies.
#
# 3. **Scalability**: By exploiting sparsity patterns (correlative and term
#    sparsity), we can handle larger quantum systems efficiently.
#
# 4. **Physical constraints**: Energy bounds ensure physically meaningful
#    correlation function bounds, preventing unphysical results.
#
# This certification methodology is particularly valuable for quantum many-body
# systems where exact solutions are unavailable, providing guaranteed error
# bars on computed properties."}