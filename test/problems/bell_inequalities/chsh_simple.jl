# test/problems/bell_inequalities/chsh_simple.jl
# Tests: CHSH Bell inequality - basic sparsity configurations (order=1)
#
# Coverage: Dense, Correlative Sparsity (MF), Term Sparsity (MMD)
# Expected optimal value: -2sqrt(2) ≈ -2.8284 (quantum bound)

using Test, NCTSSoS, JuMP

# Expectations in test/data/expectations/chsh_simple.json

function create_chsh_problem()
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    pop = polyopt(-f, reg)
    return pop, reg
end

@testset "CHSH Simple (order=1)" begin

    @testset "Dense" begin
        oracle = expectations_oracle("expectations/chsh_simple.json", "Dense_d1")
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Correlative Sparsity (MF)" begin
        oracle = expectations_oracle("expectations/chsh_simple.json", "CS_d1")
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Term Sparsity (MMD)" begin
        oracle = expectations_oracle("expectations/chsh_simple.json", "TS_d1")
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

end
