# =============================================================================
# test/problems/bell_inequalities/chsh_simple.jl
# =============================================================================
# Tests: CHSH Bell inequality - basic sparsity configurations (order=1)
# Dependencies: SOLVER
# Requires --local: no
#
# Coverage: Dense, Correlative Sparsity (MF), Term Sparsity (MMD)
# Expected optimal value: -2sqrt(2) ≈ -2.8284 (quantum bound)
# =============================================================================

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
const CHSH_SIMPLE_ORACLES = (
    Dense_d1 = (opt=-2.8284271321623193, sides=[5], nuniq=11),
    CS_d1    = (opt=-2.8284271247170496, sides=[4, 4], nuniq=10),
    TS_d1    = (opt=-2.8284271247321175, sides=[3, 3, 1], nuniq=6),
)

flatten_sizes(sizes) = reduce(vcat, sizes)

function create_chsh_problem()
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    pop = polyopt(-f, reg)
    return pop, reg
end

@testset "CHSH Simple (order=1)" begin

    @testset "Dense" begin
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ CHSH_SIMPLE_ORACLES.Dense_d1.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == CHSH_SIMPLE_ORACLES.Dense_d1.sides
        @test result.n_unique_moment_matrix_elements == CHSH_SIMPLE_ORACLES.Dense_d1.nuniq
    end

    @testset "Correlative Sparsity (MF)" begin
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ CHSH_SIMPLE_ORACLES.CS_d1.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == CHSH_SIMPLE_ORACLES.CS_d1.sides
        @test result.n_unique_moment_matrix_elements == CHSH_SIMPLE_ORACLES.CS_d1.nuniq
    end

    @testset "Term Sparsity (MMD)" begin
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ CHSH_SIMPLE_ORACLES.TS_d1.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == CHSH_SIMPLE_ORACLES.TS_d1.sides
        @test result.n_unique_moment_matrix_elements == CHSH_SIMPLE_ORACLES.TS_d1.nuniq
    end

end
