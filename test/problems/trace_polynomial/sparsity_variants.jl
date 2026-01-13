# Sparsity Algorithm Variants Tests
# Tests different term sparsity algorithms on trace polynomial optimization.
# Uses CHSH trace polynomial (expected: -2√2 ≈ -2.8284).

using Test, NCTSSoS, JuMP

# Solver: use Mosek if available, otherwise error
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

# Expected: -2√2 ≈ -2.8284271247462307
const EXPECTED_CHSH_TRACE = -2.8284271247462307

@testset "Sparsity Algorithm Variants" begin
    @testset "Term Sparsity Algorithms" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]

        p = -1.0 * NCTSSoS.tr(x[1] * y[1]) - NCTSSoS.tr(x[1] * y[2]) - NCTSSoS.tr(x[2] * y[1]) + NCTSSoS.tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)

        for (name, algo) in [
            ("NoElimination", NoElimination()),
            ("MMD", MMD()),
            ("MaximalElimination", MaximalElimination())
        ]
            @testset "$name" begin
                config = SolverConfig(optimizer=SOLVER, order=1, cs_algo=NoElimination(), ts_algo=algo)
                result = cs_nctssos(tpop, config)
                @test result.objective ≈ EXPECTED_CHSH_TRACE atol = 1e-6
            end
        end
    end
end
