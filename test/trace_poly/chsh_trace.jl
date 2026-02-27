# test/problems/bell_inequalities/chsh_trace.jl
# Tests: CHSH Bell inequality - Trace Polynomial formulation
#
# Uses trace polynomial formulation with NCTSSoS.tr().

using Test, NCTSSoS, JuMP

# Expectations in test/data/expectations/chsh_trace.json

@testset "CHSH Trace Polynomial" begin
    reg, (vars,) = create_unipotent_variables([("v", 1:4)])
    x = vars[1:2]
    y = vars[3:4]

    @testset "Variable Declaration" begin
        # Tracial Bell formulation uses the transpose trick: Bob observables appear as Bᵀ,
        # so we must NOT impose commutation between Alice and Bob variables.
        @test x[1] * y[1] != y[1] * x[1]
        @test x[2] * y[2] != y[2] * x[2]
    end

    p = -1.0 * NCTSSoS.tr(x[1] * y[1]) - NCTSSoS.tr(x[1] * y[2]) -
        NCTSSoS.tr(x[2] * y[1]) + NCTSSoS.tr(x[2] * y[2])
    tpop = polyopt(p * one(typeof(x[1])), reg)

    @testset "Dense" begin
        oracle = expectations_oracle("expectations/chsh_trace.json", "Dense")
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(tpop, config)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Term Sparsity (MMD)" begin
        oracle = expectations_oracle("expectations/chsh_trace.json", "TS")
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(tpop, config)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end
