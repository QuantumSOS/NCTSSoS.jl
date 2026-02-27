# test/problems/nc_polynomial/nc_example1.jl
# Tests: NC Example 1 - Unconstrained NC polynomial (3 variables)
#
# Standard benchmark from NCTSSOS paper.

using Test, NCTSSoS, JuMP

# Expectations in test/data/expectations/nc_example1.json

@testset "NC Example 1 (unconstrained)" begin
    n = 3
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])

    f = 1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
        2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
        x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
        9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
        54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

    pop = polyopt(f, reg)

    @testset "Dense (order=2)" begin
        oracle = expectations_oracle("expectations/nc_example1.json", "Dense_d2")
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-4
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Dense (SOS)" begin
        oracle = expectations_oracle("expectations/nc_example1.json", "SOS_d2")
        config = SolverConfig(optimizer=SOLVER, order=2)
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ oracle.opt atol = 1e-6
    end

    @testset "Term Sparsity MMD (order=2)" begin
        oracle = expectations_oracle("expectations/nc_example1.json", "TS_d2")
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-4
        @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end
