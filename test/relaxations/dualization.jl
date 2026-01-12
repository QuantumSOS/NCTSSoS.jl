# =============================================================================
# test/relaxations/dualization.jl
# =============================================================================
# Tests: SOS vs Moment method equivalence (dualization)
# Dependencies: SOLVER
# Requires --local: no
#
# Verifies that primal (Moment) and dual (SOS) formulations give the same
# optimal value, validating the dualization implementation.
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

@testset "Dualization (SOS ≈ Moment)" begin

    @testset "CHSH (Unipotent)" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(-f, reg)

        config = SolverConfig(optimizer=SOLVER, order=1)
        result_mom = cs_nctssos(pop, config; dualize=false)
        result_sos = cs_nctssos(pop, config; dualize=true)

        @test result_mom.objective ≈ result_sos.objective atol = 1e-5
    end

    @testset "NC Example 1 (NonCommutative)" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        f = 1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
            2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
            x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
            9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
            54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]
        pop = polyopt(f, reg)

        config = SolverConfig(optimizer=SOLVER, order=2)
        result_mom = cs_nctssos(pop, config; dualize=false)
        result_sos = cs_nctssos(pop, config; dualize=true)

        @test result_mom.objective ≈ result_sos.objective atol = 1e-5
    end

    @testset "NC Example 2 Constrained" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
        g = 4.0 - x[1]^2 - x[2]^2
        h1 = x[1] * x[2] + x[2] * x[1] - 2.0
        pop = polyopt(f, reg; eq_constraints=[h1], ineq_constraints=[g])

        config = SolverConfig(optimizer=SOLVER, order=2)
        result_mom = cs_nctssos(pop, config; dualize=false)
        result_sos = cs_nctssos(pop, config; dualize=true)

        @test result_mom.objective ≈ result_sos.objective atol = 1e-5
    end

end
