# =============================================================================
# test/problems/nc_polynomial/nc_example1.jl
# =============================================================================
# Tests: NC Example 1 - Unconstrained NC polynomial (3 variables)
# Dependencies: SOLVER
# Requires --local: no
#
# Standard benchmark from NCTSSOS paper.
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
const NC_EXAMPLE1_ORACLES = (
    Dense_d2 = (opt=-3.936210849199504e-10, sides=[13], nuniq=73),
    TS_d2    = (opt=-0.0035512846020968616, sides=[1, 1, 2, 2, 2, 2, 3, 3, 3, 4], nuniq=21),
)

flatten_sizes(sizes) = reduce(vcat, sizes)

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
        oracle = NC_EXAMPLE1_ORACLES.Dense_d2
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Dense (SOS)" begin
        config = SolverConfig(optimizer=SOLVER, order=2)
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ 0.0 atol = 1e-6
    end

    @testset "Term Sparsity MMD (order=2)" begin
        oracle = NC_EXAMPLE1_ORACLES.TS_d2
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
