# test/problems/nc_polynomial/nc_example2.jl
# Tests: NC Example 2 - Constrained NC polynomial (2 variables)
#
# Standard benchmark from NCTSSOS paper.
# Expected optimal value: -1.0

using Test, NCTSSoS, JuMP

@testset "NC Example 2 (constrained)" begin
    n = 2
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])

    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0

    pop = polyopt(f, reg; eq_constraints=[h1], ineq_constraints=[g])

    @testset "Dense (Moment)" begin
        config = SolverConfig(optimizer=SOLVER, order=2)
        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ -1.0 atol = 1e-6
    end

    @testset "Dense (SOS)" begin
        config = SolverConfig(optimizer=SOLVER, order=2)
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ -1.0 atol = 1e-6
    end

    @testset "Term Sparsity (Moment)" begin
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MF(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ -1.0 atol = 1e-6
    end

    @testset "Term Sparsity (SOS)" begin
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ -1.0 atol = 1e-6
    end

    @testset "cs_nctssos_higher" begin
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MF(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)
        result_higher = cs_nctssos_higher(pop, result, config)
        @test result.objective ≈ result_higher.objective atol = 1e-4
    end
end
