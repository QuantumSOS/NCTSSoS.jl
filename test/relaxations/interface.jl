# High-Level API Interface Tests
# Tests the user-facing API: polyopt, cs_nctssos, cs_nctssos_higher
# Includes:
#   - PolyOpt constructor tests (Polynomial and NCStatePolynomial)
#   - Basic dualization tests
#
# Note: Problem-specific optimization tests are in problems/ subdirectory.

using Test, NCTSSoS, JuMP

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

# PolyOpt Constructor Tests

@testset "PolyOpt Constructor" begin
    @testset "Pauli Algebra Example" begin
        N = 1
        reg, (sx, sy, sz) = create_pauli_variables(1:N)

        ham = sum(ComplexF64(1 / 2) * op[1] for op in [sx, sy, sz])

        pop = polyopt(ham, reg)
        @test pop.objective == ham
        @test isempty(pop.eq_constraints)
    end

    @testset "NonCommutative Algebra Basic" begin
        nvars = 10
        ncons = 3
        reg, (x,) = create_noncommutative_variables([("x", 1:nvars)])

        objective = sum(x .^ 2)
        constraints = [sum(Float64(i) .* x) for i = 1:ncons]

        @testset "Unconstrained" begin
            pop = polyopt(objective, reg)

            @test isempty(pop.eq_constraints)
            @test isempty(pop.ineq_constraints)
            @test length(pop.registry) == nvars
        end

        @testset "Constrained Optimization Problem" begin
            pop = polyopt(objective, reg; ineq_constraints=constraints)

            @test pop.ineq_constraints == constraints
            @test isempty(pop.eq_constraints)

            # Add an additional non-zero constraint
            extra_constraint = 1.0 * x[1]
            all_constraints = [constraints; extra_constraint]
            pop = polyopt(objective, reg; ineq_constraints=all_constraints)
            @test length(pop.ineq_constraints) == length(all_constraints)

            pop = polyopt(
                objective,
                reg;
                eq_constraints=constraints[2:2:end],
                ineq_constraints=constraints[1:2:end],
            )

            @test length(pop.eq_constraints) == 1
            @test length(pop.ineq_constraints) == 2
        end
    end

    @testset "Unipotent Algebra" begin
        reg, (x,) = create_unipotent_variables([("x", 1:5)])

        objective = sum(x .^ 2)
        pop = polyopt(objective, reg)
    end

    @testset "Projector Algebra" begin
        reg, (P,) = create_projector_variables([("P", 1:5)])

        objective = 1.0 * sum(P .^ 2)
        pop = polyopt(objective, reg)
    end
end

@testset "PolyOpt Constructor (NCStatePolynomial)" begin
    @testset "Example 7.2.1" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])

        # ς(Monomial) → StateWord, 1.0 * StateWord → StatePolynomial
        sp1 = sum([1.0, 1.0] .* map(a -> prod(ς.(a)), [[x[1] * y[2]], [y[2] * x[1]]]))
        sp2 = sum([1.0, -1.0] .* map(a -> prod(ς.(a)), [[x[1] * y[1]], [x[2] * y[2]]]))
        sp = sum([sp1 * sp1, sp2 * sp2])

        sp1_sq = sum(
            [1.0, 1.0, 1.0, 1.0] .* map(
                a -> prod(ς.(a)),
                [
                    [x[1] * y[2], x[1] * y[2]],
                    [y[2] * x[1], y[2] * x[1]],
                    [x[1] * y[2], y[2] * x[1]],
                    [y[2] * x[1], x[1] * y[2]],
                ],
            ),
        )
        sp2_sq = sum(
            [1.0, -1.0, -1.0, 1.0] .* map(
                a -> prod(ς.(a)),
                [
                    [x[1] * y[1], x[1] * y[1]],
                    [x[1] * y[1], x[2] * y[2]],
                    [x[2] * y[2], x[1] * y[1]],
                    [x[2] * y[2], x[2] * y[2]],
                ],
            ),
        )
        true_obj = sum([sp1_sq, sp2_sq])

        # StatePolynomial * NormalMonomial → NCStatePolynomial
        pop = polyopt(sp * one(NormalMonomial{UnipotentAlgebra,UInt8}), reg)
        @test pop.objective == true_obj * one(NormalMonomial{UnipotentAlgebra,UInt8})
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
    end

    @testset "Example 7.2.2" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])

        cov(a, b) = 1.0 * ς(x[a] * y[b]) - ς(x[a]) * ς(y[b])
        sp =
            cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) +
            cov(3, 1) - cov(3, 2)

        pop = polyopt(sp * one(NormalMonomial{UnipotentAlgebra,UInt8}), reg)
        true_obj = sum(
            [
                1.0,
                -1.0,
                1.0,
                -1.0,
                1.0,
                -1.0,
                1.0,
                -1.0,
                1.0,
                -1.0,
                -1.0,
                1.0,
                1.0,
                -1.0,
                -1.0,
                1.0,
            ] .* map(
                a -> prod(ς.(a)),
                ([
                    [x[1] * y[1]],
                    [x[1], y[1]],
                    [x[1] * y[2]],
                    [x[1], y[2]],
                    [x[1] * y[3]],
                    [x[1], y[3]],
                    [x[2] * y[1]],
                    [x[2], y[1]],
                    [x[2] * y[2]],
                    [x[2], y[2]],
                    [x[2] * y[3]],
                    [x[2], y[3]],
                    [x[3] * y[1]],
                    [x[3], y[1]],
                    [x[3] * y[2]],
                    [x[3], y[2]],
                ]),
            ),
        )
        @test pop.objective == true_obj * one(NormalMonomial{UnipotentAlgebra,UInt8})
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
    end

    @testset "Example 8.1.2" begin
        reg, (A, B) = create_unipotent_variables([("A", 1:3), ("B", 1:3)])

        J1 = 0.5 * (ς(A[1]) + ς(A[2]) + ς(A[3]) + ς(B[1] * A[1]) + ς(B[1] * A[2]))

        J2 =
            0.5 * (ς(A[1]) + ς(A[2]) - ς(A[3]) + ς(B[2] * A[1]) - ς(B[2] * A[2])) +
            0.5 * (
                ς(A[1] * B[3] * A[1]) - ς(A[1] * B[3] * A[2]) - ς(A[2] * B[3] * A[1]) +
                ς(A[2] * B[3] * A[2])
            )


        L = 4.0 + ς(A[1]) + ς(A[2])

        sp = sum([
            2.0 * J1 * J2,
            2.0 * J1 * L,
            2.0 * J2 * L,
            -1.0 * J1 * J1,
            -1.0 * J2 * J2,
            -1.0 * L * L,
        ])
        pop = polyopt(sp * one(NormalMonomial{UnipotentAlgebra,UInt8}), reg)
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
    end
end

# Basic Dualization Tests

@testset "Dualization" begin
    @testset "Naive Example" begin
        N = 1
        registry, (sx, sy, sz) = create_pauli_variables(1:N)

        ham = sum(ComplexF64(1 / 2) * op[1] for op in [sx, sy, sz])

        pop = polyopt(ham, registry)

        solver_config = SolverConfig(optimizer=SOLVER, order=1)

        # Both dualize=true and dualize=false now work for complex (Pauli) algebra
        res_mom = cs_nctssos(pop, solver_config; dualize=false)
        res_sos = cs_nctssos(pop, solver_config; dualize=true)
        # Both should give the same result
        @test res_mom.objective ≈ res_sos.objective atol = 1e-6
        @test res_sos.objective ≈ -0.8660254037844387 atol = 1e-6
    end

    @testset "Trivial Example" begin
        n = 2
        true_min = 3.0
        registry, (x,) = create_noncommutative_variables([("x", 1:n)])

        f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min

        pop = polyopt(f, registry)
        order = 2

        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order
        )

        result = cs_nctssos(pop, solver_config; dualize=true)

        @test isapprox(result.objective, true_min, atol=1e-6)
    end

    # This test requires high precision solver - COSMO gives Inf for one method
    if USE_LOCAL
        @testset "With Constraints" begin
            n = 2
            true_min = 3.0
            registry, (x,) = create_noncommutative_variables([("x", 1:n)])

            f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min
            r = -10.0
            g1 = r - x[1]
            g2 = r - x[2]
            g3 = x[1] - r
            g4 = x[2] - r

            pop = polyopt(f, registry; ineq_constraints=[g1, g2, g3, g4])
            order = 2

            solver_config = SolverConfig(
                optimizer=SOLVER,
                order=order
            )

            result_mom = cs_nctssos(pop, solver_config; dualize=false)
            result_sos = cs_nctssos(pop, solver_config; dualize=true)

            @test isapprox(result_mom.objective, result_sos.objective, atol=1e-3)
        end
    end
end

# Unit Tests for Extracted Helper Functions

@testset "compute_relaxation_order" begin
    @testset "Auto-compute from polynomial degree" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Degree 2 polynomial → order 1
        obj_deg2 = 1.0*x[1]^2 + 1.0*x[2]
        pop = polyopt(obj_deg2, reg)
        @test NCTSSoS.compute_relaxation_order(pop, 0) == 1

        # Degree 4 polynomial → order 2
        obj_deg4 = 1.0*x[1]^4 + 1.0*x[2]^2
        pop = polyopt(obj_deg4, reg)
        @test NCTSSoS.compute_relaxation_order(pop, 0) == 2

        # Degree 3 polynomial → order 2 (ceil(3/2))
        obj_deg3 = 1.0*x[1]^3 + 1.0*x[2]
        pop = polyopt(obj_deg3, reg)
        @test NCTSSoS.compute_relaxation_order(pop, 0) == 2
    end

    @testset "User-specified order takes precedence" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        obj = 1.0*x[1]^2  # degree 2 → auto would be 1
        pop = polyopt(obj, reg)

        @test NCTSSoS.compute_relaxation_order(pop, 3) == 3
        @test NCTSSoS.compute_relaxation_order(pop, 5) == 5
    end

    @testset "Constraints affect order" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        obj = 1.0 * x[1]  # degree 1
        constraint = 1.0*x[1]^4 + 1.0*x[2]^4  # degree 4
        pop = polyopt(obj, reg; ineq_constraints=[constraint])

        # Should use max degree across all polynomials
        @test NCTSSoS.compute_relaxation_order(pop, 0) == 2
    end

    @testset "Trivial polynomial defaults to order 1" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        obj = 5.0 * one(x[1])  # constant, degree 0
        pop = polyopt(obj, reg)

        # ceil(0/2) = 0, but we use max(1, ...) to default to 1 for trivial problems
        @test NCTSSoS.compute_relaxation_order(pop, 0) == 1
    end
end

@testset "project_to_clique" begin
    @testset "Polynomial projection" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:4)])

        # poly = x1*x2 + x2*x3 + x3*x4
        poly = 1.0*x[1]*x[2] + 1.0*x[2]*x[3] + 1.0*x[3]*x[4]

        # Get actual variable indices from the monomials
        idx1 = collect(variable_indices(x[1]))[1]
        idx2 = collect(variable_indices(x[2]))[1]
        idx3 = collect(variable_indices(x[3]))[1]
        idx4 = collect(variable_indices(x[4]))[1]

        # Clique {x1,x2} should keep only x1*x2
        proj_12 = NCTSSoS.project_to_clique(poly, [idx1, idx2])
        @test proj_12 == 1.0*x[1]*x[2]

        # Clique {x2,x3} should keep only x2*x3
        proj_23 = NCTSSoS.project_to_clique(poly, [idx2, idx3])
        @test proj_23 == 1.0*x[2]*x[3]

        # Clique {x1,x2,x3} should keep x1*x2 + x2*x3
        proj_123 = NCTSSoS.project_to_clique(poly, [idx1, idx2, idx3])
        @test proj_123 == 1.0*x[1]*x[2] + 1.0*x[2]*x[3]
    end

    @testset "Empty projection returns zero" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        poly = 1.0*x[1]*x[2]

        # Get variable index for x3
        idx3 = collect(variable_indices(x[3]))[1]

        # Clique {x3} has no overlap with variables in poly (x1*x2)
        proj = NCTSSoS.project_to_clique(poly, [idx3])
        @test iszero(proj)
    end

    @testset "Full projection returns original" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        poly = 1.0*x[1]*x[2] + 2.0*x[2]*x[3]

        # Get all variable indices
        idx1 = collect(variable_indices(x[1]))[1]
        idx2 = collect(variable_indices(x[2]))[1]
        idx3 = collect(variable_indices(x[3]))[1]

        # Clique containing all variables
        proj = NCTSSoS.project_to_clique(poly, [idx1, idx2, idx3])
        @test proj == poly
    end

    @testset "NCStatePolynomial projection" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])

        # Create state polynomial with terms on different variable sets
        sp = 1.0*ς(x[1])*ς(y[1]) + 2.0*ς(x[2])*ς(y[2]) + 3.0*ς(x[3])*ς(y[3])
        ncstp = sp * one(NormalMonomial{UnipotentAlgebra,UInt8})

        # Get variable indices for x1,x2,y1,y2
        idx_x1 = collect(variable_indices(x[1]))[1]
        idx_x2 = collect(variable_indices(x[2]))[1]
        idx_y1 = collect(variable_indices(y[1]))[1]
        idx_y2 = collect(variable_indices(y[2]))[1]

        proj = NCTSSoS.project_to_clique(ncstp, [idx_x1, idx_x2, idx_y1, idx_y2])

        # Should have 2 terms (x1*y1 and x2*y2)
        @test length(NCTSSoS.monomials(proj)) == 2
    end
end

@testset "_check_solver_status" begin
    using JuMP

    @testset "Acceptable statuses don't throw" begin
        # Create a simple feasible model
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        obj = x[1]^2 + x[2]^2 + 1.0
        pop = polyopt(obj, reg)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)

        # This should solve successfully and not throw
        result = cs_nctssos(pop, solver_config)
        @test result.objective ≈ 1.0 atol=1e-4
    end

    @testset "Status constants are defined" begin
        # Verify the acceptable statuses set exists and contains expected values
        @test MOI.OPTIMAL ∈ NCTSSoS._ACCEPTABLE_STATUSES
        @test MOI.ALMOST_OPTIMAL ∈ NCTSSoS._ACCEPTABLE_STATUSES
        @test MOI.LOCALLY_SOLVED ∈ NCTSSoS._ACCEPTABLE_STATUSES

        # These should NOT be acceptable
        @test MOI.INFEASIBLE ∉ NCTSSoS._ACCEPTABLE_STATUSES
        @test MOI.DUAL_INFEASIBLE ∉ NCTSSoS._ACCEPTABLE_STATUSES
        @test MOI.NUMERICAL_ERROR ∉ NCTSSoS._ACCEPTABLE_STATUSES
    end
end

@testset "compute_sparsity" begin
    @testset "Returns SparsityResult with correct fields" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        obj = 1.0*x[1]^2 + 1.0*x[2]^2 + 1.0*x[3]^2
        pop = polyopt(obj, reg)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)

        sparsity = compute_sparsity(pop, solver_config)

        # Check struct type
        @test sparsity isa SparsityResult

        # Check all fields are populated
        @test !isempty(sparsity.corr_sparsity.cliques)
        @test !isempty(sparsity.initial_activated_supps)
        @test !isempty(sparsity.cliques_term_sparsities)

        # Number of cliques should match
        n_cliques = length(sparsity.corr_sparsity.cliques)
        @test length(sparsity.initial_activated_supps) == n_cliques
        @test length(sparsity.cliques_term_sparsities) == n_cliques
    end

    @testset "Can inspect sparsity before solving" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:4)])

        # Create a sparse problem: x1*x2 + x3*x4 (two disconnected cliques)
        obj = 1.0*x[1]*x[2] + 1.0*x[3]*x[4]
        pop = polyopt(obj, reg)
        solver_config = SolverConfig(optimizer=SOLVER, order=1, cs_algo=MF())

        # Compute sparsity without solving
        sparsity = compute_sparsity(pop, solver_config)

        # With MF elimination, should detect 2 cliques
        @test length(sparsity.corr_sparsity.cliques) >= 1

        # initial_activated_supps should be accessible
        for (i, supp) in enumerate(sparsity.initial_activated_supps)
            @test supp isa Vector{<:NormalMonomial}
        end
    end

    @testset "Sparsity matches what cs_nctssos uses" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        obj = x[1]^2 + x[2]^2 + 1.0
        pop = polyopt(obj, reg)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)

        # Get sparsity directly
        sparsity = compute_sparsity(pop, solver_config)

        # Solve and get result
        result = cs_nctssos(pop, solver_config)

        # Sparsity in result should match
        @test result.sparsity.corr_sparsity.cliques == sparsity.corr_sparsity.cliques
        @test length(result.sparsity.cliques_term_sparsities) == length(sparsity.cliques_term_sparsities)
    end

    @testset "NCStatePolynomial PolyOpt returns SparsityResult with StateType" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])

        sp = 1.0*ς(x[1])*ς(y[1]) + 1.0*ς(x[2])*ς(y[2])
        pop = polyopt(sp * one(NormalMonomial{UnipotentAlgebra,UInt8}), reg)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)

        sparsity = compute_sparsity(pop, solver_config)

        @test sparsity isa SparsityResult
        # Verify it's a state polynomial sparsity (ST != Nothing)
        @test sparsity isa SparsityResult{<:AlgebraType, <:Integer, <:Any, <:Any, <:NCTSSoS.StateType}
        @test !isempty(sparsity.initial_activated_supps)
    end
end
