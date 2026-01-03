# =============================================================================
# High-Level API Interface Tests
# =============================================================================
# Tests the user-facing API: polyopt, cs_nctssos, cs_nctssos_higher
# Includes:
#   - PolyOpt constructor tests
#   - StatePolyOpt constructor tests
#   - Basic dualization tests
#
# Note: Problem-specific optimization tests are in problems/ subdirectory.
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "standalone_setup.jl"))

# =============================================================================
# PolyOpt Constructor Tests
# =============================================================================

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

@testset "StatePolyOpt Constructor" begin
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

        # StatePolynomial * Monomial → NCStatePolynomial
        pop = polyopt(sp * one(Monomial{UnipotentAlgebra,UInt8}), reg)
        @test pop.objective == true_obj * one(Monomial{UnipotentAlgebra,UInt8})
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
    end

    @testset "Example 7.2.2" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])

        cov(a, b) = 1.0 * ς(x[a] * y[b]) - ς(x[a]) * ς(y[b])
        sp =
            cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) +
            cov(3, 1) - cov(3, 2)

        pop = polyopt(sp * one(Monomial{UnipotentAlgebra,UInt8}), reg)
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
        @test pop.objective == true_obj * one(Monomial{UnipotentAlgebra,UInt8})
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
        pop = polyopt(sp * one(Monomial{UnipotentAlgebra,UInt8}), reg)
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
    end
end

# =============================================================================
# Basic Dualization Tests
# =============================================================================

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
