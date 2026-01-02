# =============================================================================
# Sparsity Algorithm Variants Tests
# =============================================================================
# Tests different term sparsity algorithms on trace polynomial optimization.
# Moved from relaxations/sparsity.jl as this tests actual problem solving.
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "..", "setup.jl"))

@testset "Sparsity Algorithm Variants" begin
    @testset "Term Sparsity Algorithms" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]

        p = -1.0 * NCTSSoS.tr(x[1] * y[1]) - NCTSSoS.tr(x[1] * y[2]) - NCTSSoS.tr(x[2] * y[1]) + NCTSSoS.tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)

        expected = -2.8284

        for (name, algo) in [
            ("NoElimination", NoElimination()),
            ("MMD", MMD()),
            ("MaximalElimination", MaximalElimination())
        ]
            config = SolverConfig(optimizer=SOLVER, order=1, cs_algo=NoElimination(), ts_algo=algo)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ expected atol = 1e-4
        end
    end
end
