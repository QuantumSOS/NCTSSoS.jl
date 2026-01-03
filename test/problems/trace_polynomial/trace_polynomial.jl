# =============================================================================
# Trace Polynomial Optimization Tests (Examples 6.x)
# =============================================================================
# Consolidates trace polynomial tests using the tr() operator:
#   - Example 6.1: Projector algebra with product of traces
#   - Example 6.2.1: Squared trace expressions
#   - Example 6.2.2: Covariance trace inequality
#
# Note: Basic CHSH trace polynomial tests are in chsh.jl
# Results verified against NCTSSOS oracles (trace_poly_oracles.jl)
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__DIR__), "..", "standalone_setup.jl"))

# Load oracle values
include(joinpath(dirname(@__DIR__), "..", "oracles", "results", "trace_poly_oracles.jl"))

@testset "Trace Polynomial Examples (6.x)" begin

    # =========================================================================
    # Example 6.1: Projector Algebra
    # =========================================================================
    @testset "Example 6.1 (Projector Algebra)" begin
        reg, (x,) = create_projector_variables([("x", 1:3)])

        p = (NCTSSoS.tr(x[1] * x[2] * x[3]) + NCTSSoS.tr(x[1] * x[2]) * NCTSSoS.tr(x[3])) * one(typeof(x[1]))
        tpop = polyopt(p, reg)

        oracle_d2 = TRACE_POLY_ORACLES["Trace_6_1_Dense_d2"]
        oracle_d3 = TRACE_POLY_ORACLES["Trace_6_1_Dense_d3"]

        @testset "Order 2" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracle_d2.opt atol = 1e-6
        end

        if USE_LOCAL
            @testset "Order 3" begin
                config = SolverConfig(optimizer=SOLVER, order=3)
                result = cs_nctssos(tpop, config)
                @test result.objective ≈ oracle_d3.opt atol = 1e-6
            end
        end
    end

    # =========================================================================
    # Example 6.2.1: Squared Trace Expressions
    # =========================================================================
    @testset "Example 6.2.1 (Squared Traces)" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]

        p = (1.0 * NCTSSoS.tr(x[1] * y[2]) + NCTSSoS.tr(x[2] * y[1])) *
            (1.0 * NCTSSoS.tr(x[1] * y[2]) + NCTSSoS.tr(x[2] * y[1])) +
            (1.0 * NCTSSoS.tr(x[1] * y[1]) - NCTSSoS.tr(x[2] * y[2])) *
            (1.0 * NCTSSoS.tr(x[1] * y[1]) - NCTSSoS.tr(x[2] * y[2]))

        tpop = polyopt((-1.0 * p) * one(typeof(x[1])), reg)
        oracle = TRACE_POLY_ORACLES["Trace_6_2_1_Dense_d2"]

        # Order=2 gives the correct tight bound of -4.0 (requires Mosek)
        if USE_LOCAL
            @testset "Order 2 (tight bound)" begin
                config = SolverConfig(optimizer=SOLVER, order=2)
                result = cs_nctssos(tpop, config)
                @test result.objective ≈ oracle.opt atol = 1e-5
            end
        end
    end

    # =========================================================================
    # Example 6.2.2: Covariance Trace Inequality
    # =========================================================================
    @testset "Example 6.2.2 (Covariance)" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:6)])
        x = vars[1:3]
        y = vars[4:6]

        cov(i, j) = NCTSSoS.tr(x[i] * y[j]) - NCTSSoS.tr(x[i]) * NCTSSoS.tr(y[j])
        p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) +
                   cov(2, 1) + cov(2, 2) - cov(2, 3) +
                   cov(3, 1) - cov(3, 2))

        tpop = polyopt(p * one(typeof(x[1])), reg)
        oracle_dense = TRACE_POLY_ORACLES["Trace_6_2_2_Dense_d2"]
        oracle_ts = TRACE_POLY_ORACLES["Trace_6_2_2_TS_d2"]

        @testset "Dense" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracle_dense.opt atol = 1e-5
        end

        @testset "Sparse (MF + MMD)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracle_ts.opt atol = 1e-4
        end
    end
end
