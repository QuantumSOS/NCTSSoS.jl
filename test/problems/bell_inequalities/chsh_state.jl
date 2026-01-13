# test/problems/bell_inequalities/chsh_state.jl
# Tests: CHSH Bell inequality - State Polynomial formulation
#
# Uses state polynomial formulation with ς() state function.

using Test, NCTSSoS, JuMP

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

# Oracle values from NCTSSOS
const CHSH_STATE_ORACLES = (
    Dense = (opt=-2.828427124746234, sides=[9], nuniq=21),
    TS    = (opt=-2.8284271247321175, nuniq=10),  # sides vary by implementation
)

if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

@testset "CHSH State Polynomial" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp = -ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
    spop = polyopt(sp * one(typeof(x[1])), reg)

    @testset "Dense" begin
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(spop, config)
        @test result.objective ≈ CHSH_STATE_ORACLES.Dense.opt atol = 1e-5
        @test flatten_sizes(result.moment_matrix_sizes) == CHSH_STATE_ORACLES.Dense.sides
        @test result.n_unique_moment_matrix_elements == CHSH_STATE_ORACLES.Dense.nuniq
    end

    @testset "Term Sparsity (MMD)" begin
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(spop, config)
        @test result.objective ≈ CHSH_STATE_ORACLES.TS.opt atol = 1e-5
        @test result.n_unique_moment_matrix_elements == CHSH_STATE_ORACLES.TS.nuniq
    end
end
