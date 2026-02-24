# test/problems/bell_inequalities/chsh_trace.jl
# Tests: CHSH Bell inequality - Trace Polynomial formulation
#
# Uses trace polynomial formulation with NCTSSoS.tr().

using Test, NCTSSoS, JuMP

# Oracle values from NCTSSOS
const CHSH_TRACE_ORACLES = (
    Dense = (opt=-2.828427124746234, sides=[9], nuniq=21),
    TS    = (opt=-2.8284271247321175, nuniq=10),  # sides vary by implementation
)

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
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(tpop, config)
        @test result.objective ≈ CHSH_TRACE_ORACLES.Dense.opt atol = 1e-5
        @test flatten_sizes(result.moment_matrix_sizes) == CHSH_TRACE_ORACLES.Dense.sides
        @test result.n_unique_moment_matrix_elements == CHSH_TRACE_ORACLES.Dense.nuniq
    end

    @testset "Term Sparsity (MMD)" begin
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(tpop, config)
        @test result.objective ≈ CHSH_TRACE_ORACLES.TS.opt atol = 1e-5
        @test result.n_unique_moment_matrix_elements == CHSH_TRACE_ORACLES.TS.nuniq
    end
end
