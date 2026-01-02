# =============================================================================
# State Polynomial Optimization Tests (Examples 7.2.x)
# =============================================================================
# Consolidates state polynomial tests using the ς() operator:
#   - 7.2.1: Squared expectations (products of state values)
#   - 7.2.2: Covariance Bell inequality
#   - 7.2.3: Mixed state polynomial with squared expectations
#
# Note: Basic CHSH state polynomial tests are in chsh.jl
# Results verified against NCTSSOS oracles (state_poly_extended_oracles.jl)
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__DIR__), "..", "setup.jl"))

# Load oracle values
include(joinpath(dirname(@__DIR__), "..", "oracles", "results", "state_poly_extended_oracles.jl"))

@testset "State Polynomial Examples (7.2.x)" begin

    # =========================================================================
    # Example 7.2.1: Squared Expectations
    # =========================================================================
    # Known limitation - objectives with squared expectations <A><A>
    # At order=3, the relaxation correctly gives the tight bound of -4.0.
    @testset "Example 7.2.1 (Squared Expectations)" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
        sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
        sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

        spop = polyopt(sp * one(typeof(x[1])), reg)
        oracle = STATE_POLY_EXTENDED_ORACLES["State_7_2_1_Dense_d3"]

        @testset "Order 3 (tight bound)" begin
            config = SolverConfig(optimizer=SOLVER, order=3)
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracle.opt atol = 1e-4
        end
    end

    # =========================================================================
    # Example 7.2.2: Covariance Bell Inequality
    # =========================================================================
    @testset "Example 7.2.2 (Covariance)" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
        cov(a, b) = 1.0 * ς(x[a] * y[b]) - 1.0 * ς(x[a]) * ς(y[b])
        sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)

        spop = polyopt(sp * one(typeof(x[1])), reg)
        oracle_dense = STATE_POLY_EXTENDED_ORACLES["State_7_2_2_Dense_d2"]
        oracle_ts = STATE_POLY_EXTENDED_ORACLES["State_7_2_2_TS_d2"]

        @testset "Dense" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracle_dense.opt atol = 1e-2
        end

        @testset "Sparse (MF + MMD)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracle_ts.opt atol = 1e-2
        end
    end

    # =========================================================================
    # Example 7.2.3: Mixed State Polynomial
    # =========================================================================
    # Mixed state polynomial with products of expectations ς(a)*ς(b) and ς(a)²
    # From NCTSSOS example with n=4 variables
    @testset "Example 7.2.3 (Mixed)" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        # Coefficients from NCTSSOS
        sp = -1.0 * ς(x[2]) - 1.0 * ς(y[1]) - 1.0 * ς(y[2]) +    # linear terms
             1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[2] * y[1]) -        # mixed monomial terms
             1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[2]) +        # more mixed
             1.0 * ς(x[1]) * ς(y[1]) + 1.0 * ς(x[2]) * ς(y[1]) +  # products of expectations
             1.0 * ς(x[2]) * ς(y[2]) +                             # more products
             1.0 * ς(x[1]) * ς(x[1]) + 1.0 * ς(y[2]) * ς(y[2])    # squared expectations

        spop = polyopt(sp * one(typeof(x[1])), reg)
        oracle_dense = STATE_POLY_EXTENDED_ORACLES["State_7_2_3_Dense_d2"]
        oracle_ts = STATE_POLY_EXTENDED_ORACLES["State_7_2_3_TS_d2"]

        @testset "Dense (Moment)" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(spop, config; dualize=false)
            @test result.objective ≈ oracle_dense.opt atol = 1e-2
        end

        @testset "Dense (SOS)" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(spop, config; dualize=true)
            @test result.objective ≈ oracle_dense.opt atol = 1e-2
        end

        @testset "Sparse (MMD)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracle_ts.opt atol = 1e-2
        end

        @testset "Moment vs SOS Consistency" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result_mom = cs_nctssos(spop, config; dualize=false)
            result_sos = cs_nctssos(spop, config; dualize=true)
            @test result_mom.objective ≈ result_sos.objective atol = 1e-3
        end
    end
end
