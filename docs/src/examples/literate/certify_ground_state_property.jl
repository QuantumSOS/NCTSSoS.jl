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

@show ground_state_energy_upper_bounds
@show s1s2values

# ## Computing Lower Bounds with Semidefinite Programming
#
# Now we obtain rigorous lower bounds on the ground-state energy using
# `NCTSSoS.jl`. This approach formulates the quantum ground-state problem
# as a polynomial optimization problem and solves its semidefinite
# programming relaxation, providing certified lower bounds.
#
# The NCTSSoS API provides type-safe Pauli algebra with automatic
# simplification. No explicit commutation constraints are needed - the
# `PauliAlgebra` type handles the Pauli operator relations automatically.

using NCTSSoS
using MosekTools
using JuMP

const MOI = NCTSSoS.MOI

## Configure Mosek solver with high precision settings
SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    MOI.Silent() => true,
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
    ## Create Pauli algebra variables
    ## PauliAlgebra automatically handles:
    ##   - σ² = I (unipotent property)
    ##   - Cyclic product rules: σx·σy = iσz, σy·σz = iσx, σz·σx = iσy
    ##   - Commutation relations between different sites
    registry, (x, y, z) = create_pauli_variables(1:N)  # x = σx, y = σy, z = σz

    ## Objective function: we want to minimize the Hamiltonian
    ## H = J/4 * Σ z_i z_{i+1} + h/2 * Σ x_i
    ham = sum(T1(J / 4) * z[i] * z[mod1(i + 1, N)] + T1(h / 2) * x[i] for i in 1:N)

    ## Create polynomial optimization problem
    ## No explicit algebra constraints needed - PauliAlgebra handles them automatically
    pop = polyopt(ham, registry)

    ## Configure solver with second-order moment relaxation
    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    ## Solve the semidefinite program to get energy lower bound
    res = cs_nctssos(pop, solver_config)

    ## Store energy per site (divide by system size)
    push!(ground_state_energy_lower_bounds, res.objective / N)
end

@show ground_state_energy_lower_bounds

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
function _check_solver_status_allow_slow(model)
    status = termination_status(model)
    status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) &&
        return status

    primal = primal_status(model)
    dual = dual_status(model)

    if status in (MOI.ITERATION_LIMIT, MOI.SLOW_PROGRESS) &&
       (primal in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT) ||
        dual in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT))
        return status
    end

    error("Solver failed: termination=$status, primal=$primal, dual=$dual")
end

function cs_nctssos_with_entry(
    pop::OP,
    solver_config::SolverConfig,
    entry_constraints::Vector{P};
    dualize::Bool=true
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, OP<:NCTSSoS.OptimizationProblem{A,P}}

   ## Compute sparsity structure (correlative + term sparsity)
   sparsity = NCTSSoS.compute_sparsity(pop, solver_config)

   ## Build the moment relaxation problem
   moment_problem = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)

   ## Add entry constraints for correlation function bounds
   ## These are semidefinite constraints that ensure physical validity
   for c in entry_constraints
       push!(moment_problem.constraints,(:HPSD, [c;;]))
   end

   if dualize
       sos_problem = NCTSSoS.sos_dualize(moment_problem)
       set_optimizer(sos_problem.model, solver_config.optimizer)
       optimize!(sos_problem.model)
       _check_solver_status_allow_slow(sos_problem.model)
       return NCTSSoS.PolyOptResult(
           objective_value(sos_problem.model),
           sparsity,
           sos_problem.model,
           moment_problem.n_unique_moment_matrix_elements
       )
   else
       result = NCTSSoS.solve_moment_problem(moment_problem, solver_config.optimizer)
       _check_solver_status_allow_slow(result.model)
       return NCTSSoS.PolyOptResult(result.objective, sparsity, result.model, result.n_unique_elements)
   end
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
    ## Create Pauli algebra variables
    registry, (x, y, z) = create_pauli_variables(1:N)

    ## Objective function: S^z_1 * S^z_2 correlation
    obj = one(T1) * z[1] * z[2]

    ## Hamiltonian (same form as before)
    ham = sum(T1(J / 4) * z[j] * z[mod1(j + 1, N)] + T1(h / 2) * x[j] for j in 1:N)

    ## Energy constraint: ensure ground state energy is within reference bounds
    ## This is crucial for obtaining physically meaningful correlation bounds
    ineq_cons = [ham - ground_state_energy_lower_bounds[i] * N]

    ## Create optimization problems for lower and upper bounds
    ## We solve two separate problems: minimize and maximize the correlation
    pop_l = polyopt(obj, registry; ineq_constraints=ineq_cons)
    pop_u = polyopt(-obj, registry; ineq_constraints=ineq_cons)

    ## Configure solver
    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    ## Additional energy constraint for the upper bound problem
    single_ineq_cons = [ground_state_energy_upper_bounds[i] * N - ham]

    ## Solve for lower and upper bounds on correlation function
    res_l = cs_nctssos_with_entry(pop_l, solver_config, single_ineq_cons; dualize=false)
    res_u = cs_nctssos_with_entry(pop_u, solver_config, single_ineq_cons; dualize=false)

    ## Store bounds (divide by 4 to convert from Pauli to spin operators)
    push!(lower_bounds, res_l.objective / 4)
    push!(upper_bounds, -res_u.objective / 4)
end

@show lower_bounds
@show upper_bounds

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
# 5. **Type-safe algebra**: The `PauliAlgebra` type automatically handles
#    Pauli operator relations (σ² = I, commutation rules) without requiring
#    explicit constraints.
#
# This certification methodology is particularly valuable for quantum many-body
# systems where exact solutions are unavailable, providing guaranteed error
# bars on computed properties.
