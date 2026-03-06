using Test, NCTSSoS, JuMP

if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

function wang_magron_example_3_4()
    reg, (x,) = create_noncommutative_variables([("x", 1:3)])
    f = 1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
        2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
        x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
        9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
        54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]
    return polyopt(f, reg)
end

function wang_magron_example_3_8()
    reg, (x,) = create_noncommutative_variables([("x", 1:2)])
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h = x[1] * x[2] + x[2] * x[1] - 2.0
    return polyopt(f, reg; eq_constraints=[h], ineq_constraints=[g])
end

function wang_magron_example_5_3()
    reg, (x,) = create_noncommutative_variables([("x", 1:2)])
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    objective = (1.0 * tr(f)) * one(typeof(x[1]))
    g = (1.0 * tr(4.0 - x[1]^2 - x[2]^2)) * one(typeof(x[1]))
    h = (1.0 * tr(x[1] * x[2] + x[2] * x[1] - 2.0)) * one(typeof(x[1]))
    return polyopt(objective, reg; eq_constraints=[h], ineq_constraints=[g])
end

@testset "Wang-Magron (2021) Paper Examples" begin
    @testset "Example 3.4" begin
        pop = wang_magron_example_3_4()

        oracle_dense = expectations_oracle(
            "expectations/wang_magron_examples.json",
            "Example_3_4_Dense_d2"
        )
        dense_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        dense_result = cs_nctssos(pop, dense_config; dualize=false)
        @test dense_result.objective ≈ oracle_dense.opt atol = 1e-4
        @test flatten_sizes(dense_result.moment_matrix_sizes) == oracle_dense.sides
        @test dense_result.n_unique_moment_matrix_elements == oracle_dense.nuniq

        oracle_ts = expectations_oracle(
            "expectations/wang_magron_examples.json",
            "Example_3_4_TS_d2"
        )
        ts_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        ts_result = cs_nctssos(pop, ts_config; dualize=false)
        @test ts_result.objective ≈ oracle_ts.opt atol = 1e-4
        @test sort(flatten_sizes(ts_result.moment_matrix_sizes)) == sort(oracle_ts.sides)
        @test ts_result.n_unique_moment_matrix_elements == oracle_ts.nuniq
    end

    @testset "Example 3.8" begin
        pop = wang_magron_example_3_8()

        oracle_dense = expectations_oracle(
            "expectations/wang_magron_examples.json",
            "Example_3_8_Dense_d2"
        )
        dense_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        dense_result = cs_nctssos(pop, dense_config; dualize=false)
        @test dense_result.objective ≈ oracle_dense.opt atol = 1e-6
        @test flatten_sizes(dense_result.moment_matrix_sizes) == oracle_dense.sides
        @test dense_result.n_unique_moment_matrix_elements == oracle_dense.nuniq

        oracle_ts = expectations_oracle(
            "expectations/wang_magron_examples.json",
            "Example_3_8_TS_d2"
        )
        ts_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        ts_result = cs_nctssos(pop, ts_config; dualize=false)
        @test ts_result.objective ≈ oracle_ts.opt atol = 1e-6
        @test sort(flatten_sizes(ts_result.moment_matrix_sizes)) == sort(oracle_ts.sides)
        @test ts_result.n_unique_moment_matrix_elements == oracle_ts.nuniq
    end

    @testset "Example 5.3" begin
        pop = wang_magron_example_5_3()

        oracle_dense = expectations_oracle(
            "expectations/wang_magron_examples.json",
            "Example_5_3_Dense_d2"
        )
        dense_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        dense_result = cs_nctssos(pop, dense_config; dualize=false)
        @test dense_result.objective ≈ oracle_dense.opt atol = 1e-6
        @test flatten_sizes(dense_result.moment_matrix_sizes) == oracle_dense.sides
        @test dense_result.n_unique_moment_matrix_elements == oracle_dense.nuniq

        oracle_ts = expectations_oracle(
            "expectations/wang_magron_examples.json",
            "Example_5_3_TS_d2"
        )
        ts_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        ts_result = cs_nctssos(pop, ts_config; dualize=false)
        @test ts_result.objective ≈ oracle_ts.opt atol = 1e-6
        @test ts_result.n_unique_moment_matrix_elements == oracle_ts.nuniq
    end
end
