# =============================================================================
# test/problems/nc_polynomial/nc_readme.jl
# =============================================================================
# Tests: README documentation examples
# Dependencies: SOLVER
# Requires --local: no
#
# Validates that examples from README.md work correctly.
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

@testset "README Examples" begin

    @testset "Unconstrained" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        f = 1.0 + x[1]^4 + x[2]^4 + x[3]^4 +
            x[1] * x[2] + x[2] * x[1] + x[2] * x[3] + x[3] * x[2]

        pop = polyopt(f, reg)

        result_dense = cs_nctssos(pop, SolverConfig(optimizer=SOLVER))
        result_cs = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, cs_algo=MF()))
        result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, cs_algo=MF(), ts_algo=MMD()))

        @test result_dense.objective ≈ result_cs.objective atol = 1e-6
        @test result_cs.objective ≈ result_cs_ts.objective atol = 1e-6
    end

    @testset "Constrained" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
        g = 4.0 - x[1]^2 - x[2]^2
        h1 = x[1] * x[2] + x[2] * x[1] - 2.0

        pop = polyopt(f, reg; ineq_constraints=[g], eq_constraints=[h1])

        result_dense = cs_nctssos(pop, SolverConfig(optimizer=SOLVER))
        result_cs = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, cs_algo=MF()))
        result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, cs_algo=MF(), ts_algo=MMD()))

        @test result_dense.objective ≈ result_cs.objective atol = 1e-6
        @test result_cs.objective ≈ result_cs_ts.objective atol = 1e-6

        result_cs_ts_higher = cs_nctssos_higher(
            pop,
            result_cs_ts,
            SolverConfig(optimizer=SOLVER, cs_algo=MF(), ts_algo=MMD())
        )
        @test result_dense.objective ≈ result_cs_ts_higher.objective atol = 1e-6
    end

end
