# PolyOpt Constructor tests for state polynomial optimization problems.
#
# Moved from `test/relaxations/interface.jl` to keep state polynomial tests isolated
# under `test/state_poly/`.

using Test, NCTSSoS

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

