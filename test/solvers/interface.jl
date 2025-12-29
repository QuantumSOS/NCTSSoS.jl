# High-Level API Interface Tests
# ===============================
# Tests the user-facing API: polyopt, cs_nctssos, cs_nctssos_higher
# Includes:
#   - PolyOpt constructor tests (formerly pop.jl)
#   - StatePolyOpt constructor tests
#   - Solver integration tests

using Test, NCTSSoS

# =============================================================================
# PolyOpt Constructor Tests (merged from pop.jl)
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
# Solver Integration Tests
# =============================================================================

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

if LOCAL_TESTING
    @testset "1D Transverse Field Ising Model" begin
        N = 3
        registry, (sx, sy, sz) = create_pauli_variables(1:N)

        J = 1.0
        h = 2.0
        for (periodic, true_ans) in zip((true, false), (-1.0175918, -1.0104160))
            ham = sum(-complex(J / 4) * sz[i] * sz[mod1(i + 1, N)] for i in 1:(periodic ? N : N - 1)) + sum(-h / 2 * sx[i] for i in 1:N)

            pop = polyopt(ham, registry)

            solver_config = SolverConfig(optimizer=SOLVER, order=2)

            res = cs_nctssos(pop, solver_config)
            @test res.objective / N ≈ true_ans atol = 1e-6
        end
    end

    @testset "1D Heisenberg Chain" begin
        N = 6
        registry, (sx, sy, sz) = create_pauli_variables(1:N)

        ham = sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [sx, sy, sz] for i in 1:N)

        pop = polyopt(ham, registry)

        solver_config = SolverConfig(optimizer=SOLVER, order=2)

        res = cs_nctssos(pop, solver_config)

        @test res.objective / N ≈ -0.467129 atol = 1e-6
    end

    @testset "I_3322 Example with Sparsity" begin
        # Use projector algebra (P² = P)
        registry, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]

        pop = polyopt(-f, registry)

        for (cs_algo, ts_algo, ans) in zip([NoElimination(), MF(), MF()],
            [NoElimination(), MMD(), MaximalElimination()],
            [-0.2508755573198166, -0.9999999892255513, -0.3507010331201541])
            solver_config = SolverConfig(optimizer=SOLVER; order=3, cs_algo=cs_algo, ts_algo=ts_algo)
            result = cs_nctssos(pop, solver_config)
            @test isapprox(result.objective, ans; atol=1e-5)
        end
    end
end

@testset "Majumdar Gosh Model" begin
    num_sites = 6
    J1_interactions =
        unique!([tuple(sort([i, mod1(i + 1, num_sites)])...) for i = 1:num_sites])
    J2_interactions =
        unique!([tuple(sort([i, mod1(i + 2, num_sites)])...) for i = 1:num_sites])

    J1 = 2.0
    J2 = 1.0

    true_ans = -num_sites / 4 * 6

    ij2idx_dict = Dict(
        zip(
            [(i, j) for i in 1:num_sites, j in 1:num_sites if j > i],
            1:(num_sites*(num_sites-1)÷2),
        ),
    )

    # Use projector algebra (P² = P)
    registry, (hij,) = create_projector_variables([("h", 1:(num_sites*(num_sites-1)÷2))])

    objective = (
        sum([J1 * hij[ij2idx_dict[(i, j)]] for (i, j) in J1_interactions]) + sum([J2 * hij[ij2idx_dict[(i, j)]] for (i, j) in J2_interactions])
    )

    gs = unique!([
        (
            hij[ij2idx_dict[tuple(sort([i, j])...)]] *
            hij[ij2idx_dict[tuple(sort([j, k])...)]] +
            hij[ij2idx_dict[tuple(sort([j, k])...)]] *
            hij[ij2idx_dict[tuple(sort([i, j])...)]] -
            0.5 * (
                hij[ij2idx_dict[tuple(sort([i, j])...)]] +
                hij[ij2idx_dict[tuple(sort([j, k])...)]] -
                hij[ij2idx_dict[tuple(sort([i, k])...)]]
            )
        ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
        (i != j && j != k && i != k)
    ])

    pop = polyopt(-objective, registry; eq_constraints=gs)

    solver_config = SolverConfig(optimizer=SOLVER; order=1)

    result = cs_nctssos(pop, solver_config)

    @test isapprox(result.objective, true_ans; atol=1e-4)
end

@testset "Problem Creation Interface" begin
    n = 2
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    pop = polyopt(f, registry; ineq_constraints=[g], eq_constraints=[h1])

    solver_config = SolverConfig(
        optimizer=SOLVER;
        order=2,
        cs_algo=MF(),
        ts_algo=MMD(),
    )

    result = cs_nctssos(pop, solver_config)
    @test isapprox(result.objective, -1.0; atol=1e-4)

    result_higher = cs_nctssos_higher(pop, result, solver_config)
    @test isapprox(result.objective, result_higher.objective; atol=1e-4)
end

@testset "README Example Unconstrained" begin
    registry, (x,) = create_noncommutative_variables([("x", 1:3)])
    f =
        1.0 +
        x[1]^4 +
        x[2]^4 +
        x[3]^4 +
        x[1] * x[2] +
        x[2] * x[1] +
        x[2] * x[3] +
        x[3] * x[2]

    pop = polyopt(f, registry)

    solver_config_dense = SolverConfig(optimizer=SOLVER)

    result_dense = cs_nctssos(pop, solver_config_dense)

    result_cs =
        cs_nctssos(pop, SolverConfig(optimizer=SOLVER; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(
        pop,
        SolverConfig(optimizer=SOLVER; cs_algo=MF(), ts_algo=MMD()),
    )

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)
end

@testset "README Example Constrained" begin
    registry, (x,) = create_noncommutative_variables([("x", 1:2)])
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0

    pop = polyopt(f, registry; ineq_constraints=[g], eq_constraints=[h1])

    result_dense = cs_nctssos(pop, SolverConfig(optimizer=SOLVER))

    result_cs =
        cs_nctssos(pop, SolverConfig(optimizer=SOLVER; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(
        pop,
        SolverConfig(optimizer=SOLVER; cs_algo=MF(), ts_algo=MMD()),
    )

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)

    result_cs_ts_higher = cs_nctssos_higher(
        pop,
        result_cs_ts,
        SolverConfig(optimizer=SOLVER; cs_algo=MF(), ts_algo=MMD()),
    )

    @test isapprox(result_dense.objective, result_cs_ts_higher.objective, atol=1e-4)
end
