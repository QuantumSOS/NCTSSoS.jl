# PolyOpt Constructor tests for state polynomial optimization problems.
#
# Moved from `test/relaxations/interface.jl` to keep state polynomial tests isolated
# under `test/state_poly/`.

using Test, NCTSSoS

@testset "PolyOpt Constructor (NCStatePolynomial)" begin
    @testset "StatePolynomial input is promoted to NCStatePolynomial" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])

        objective_sp = -1.0 * ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
        eq_sp = -1.0 + ς(x[1])
        ineq_sp = 2.0 + ς(y[1])

        pop_state = polyopt(objective_sp, reg; eq_constraints=[eq_sp], ineq_constraints=[ineq_sp])
        pop_mixed = polyopt(
            objective_sp,
            reg;
            eq_constraints=[eq_sp * one(typeof(x[1]))],
            ineq_constraints=[ineq_sp * one(typeof(x[1]))],
        )
        pop_mixed_vector = polyopt(
            objective_sp,
            reg;
            eq_constraints=Any[eq_sp, eq_sp * one(typeof(x[1]))],
            ineq_constraints=Any[ineq_sp, ineq_sp * one(typeof(x[1]))],
        )
        pop_nc = polyopt(
            objective_sp * one(typeof(x[1])),
            reg;
            eq_constraints=[eq_sp * one(typeof(x[1]))],
            ineq_constraints=[ineq_sp * one(typeof(x[1]))],
        )

        @test pop_state.objective isa NCTSSoS.NCStatePolynomial{Float64,Arbitrary,UnipotentAlgebra,UInt8}
        @test pop_state.objective == pop_nc.objective
        @test pop_state.eq_constraints == pop_nc.eq_constraints
        @test pop_state.ineq_constraints == pop_nc.ineq_constraints
        @test pop_mixed.objective == pop_nc.objective
        @test pop_mixed.eq_constraints == pop_nc.eq_constraints
        @test pop_mixed.ineq_constraints == pop_nc.ineq_constraints
        @test pop_mixed_vector.objective == pop_nc.objective
        @test pop_mixed_vector.eq_constraints == pop_nc.eq_constraints
        @test pop_mixed_vector.ineq_constraints == pop_nc.ineq_constraints
        @test all(isone(ncsw.nc_word) for ncsw in monomials(pop_state.objective))

        err_bad_constraint = try
            polyopt(objective_sp, reg; eq_constraints=Any[1.0])
            nothing
        catch e
            e
        end
        @test err_bad_constraint isa ArgumentError
        @test occursin("State-polynomial constraints must be", sprint(showerror, err_bad_constraint))

        # NCStatePolynomial with non-identity NC words must be rejected
        bad_ncsw = NCTSSoS.NCStateWord(ς(x[1]), x[2])
        bad_obj = NCTSSoS.NCStatePolynomial([1.0], [bad_ncsw])
        @test_throws ArgumentError polyopt(bad_obj, reg)
        err_bad_nc = try
            polyopt(bad_obj, reg)
            nothing
        catch e
            e
        end
        @test occursin("non-identity NC word", sprint(showerror, err_bad_nc))
        @test occursin("expect", sprint(showerror, err_bad_nc))

        # Non-identity NC words in constraints are allowed (localizing matrices)
        good_obj = objective_sp * one(typeof(x[1]))
        bad_eq = NCTSSoS.NCStatePolynomial([1.0], [bad_ncsw])
        @test polyopt(good_obj, reg; eq_constraints=[bad_eq]) isa NCTSSoS.PolyOpt
        @test polyopt(good_obj, reg; ineq_constraints=[bad_eq]) isa NCTSSoS.PolyOpt

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        result_state = cs_nctssos(pop_state, config)
        result_nc = cs_nctssos(pop_nc, config)

        @test result_state.objective ≈ result_nc.objective atol = 1e-5
        @test result_state.n_unique_moment_matrix_elements == result_nc.n_unique_moment_matrix_elements
        @test result_state.moment_matrix_sizes == result_nc.moment_matrix_sizes
    end

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

        pop = polyopt(sp, reg)
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

        pop = polyopt(sp, reg)
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
        pop = polyopt(sp, reg)
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
    end
end

