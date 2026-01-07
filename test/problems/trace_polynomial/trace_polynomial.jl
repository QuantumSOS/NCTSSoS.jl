# =============================================================================
# Trace Polynomial Optimization Tests (Examples 6.x)
# =============================================================================
# Consolidates trace polynomial tests using the tr() operator:
#   - Example 6.1: Projector algebra with product of traces
#   - Example 6.2.0: CHSH trace polynomial with term sparsity variants
#   - Example 6.2.1: Squared trace expressions
#   - Example 6.2.2: Covariance trace inequality
#
# Results verified against NCTSSOS.
# =============================================================================

using Test, NCTSSoS, JuMP

# Solver: use Mosek if available, otherwise error
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
    const SOLVER_NAME = :mosek
end

# =============================================================================
# Expected values: mosek (reference) and cosmo (CI)
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = moment matrix block sizes
#   nuniq = unique moment indices (affine constraints)
# =============================================================================
const EXPECTED_TRACE_POLY = (
    mosek = (
        Ex_6_1_Dense_d2   = (opt=-0.04671737845552321, sides=[31], nuniq=81),
        Ex_6_1_Dense_d3   = (opt=-0.031249989780027937, sides=[108], nuniq=395),
        Ex_6_2_0_Dense_d1 = (opt=-2.8284271247462307,),
        Ex_6_2_1_Dense_d2 = (opt=-4.000000007251562, sides=[53], nuniq=222),
        Ex_6_2_2_Dense_d2 = (opt=-4.999999997079296, sides=[115], nuniq=1010),
        Ex_6_2_2_TS_d2    = (opt=-4.999999995608242, nuniq=199),
    ),
    cosmo = (
        Ex_6_1_Dense_d2   = (opt=-0.04671740059436102, sides=[31], nuniq=81),
        Ex_6_1_Dense_d3   = (opt=-0.031250237578871895, sides=[114], nuniq=395),
        Ex_6_2_0_Dense_d1 = (opt=-2.8284271246604606,),
        Ex_6_2_1_Dense_d2 = (opt=-4.000000001047738, sides=[53], nuniq=294),
        Ex_6_2_2_Dense_d2 = (opt=-5.000000449255744, sides=[115], nuniq=1280),
        Ex_6_2_2_TS_d2    = (opt=-4.999999996016442, nuniq=63),
    ),
)

select_oracle(oracles) = getproperty(oracles, SOLVER_NAME)

@testset "Trace Polynomial Examples (6.x)" begin
    oracles = select_oracle(EXPECTED_TRACE_POLY)

    # =========================================================================
    # Example 6.1: Projector Algebra
    # =========================================================================
    @testset "Example 6.1 (Projector Algebra)" begin
        reg, (x,) = create_projector_variables([("x", 1:3)])

        p = (NCTSSoS.tr(x[1] * x[2] * x[3]) + NCTSSoS.tr(x[1] * x[2]) * NCTSSoS.tr(x[3])) * one(typeof(x[1]))
        tpop = polyopt(p, reg)

        @testset "Order 2" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracles.Ex_6_1_Dense_d2.opt atol = 1e-6
        end

        @testset "Order 3" begin
            config = SolverConfig(optimizer=SOLVER, order=3)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracles.Ex_6_1_Dense_d3.opt atol = 1e-6
        end
    end

    # =========================================================================
    # Example 6.2.0: CHSH Trace Polynomial
    # =========================================================================
    @testset "Example 6.2.0 (CHSH)" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]

        p = -1.0 * NCTSSoS.tr(x[1] * y[1]) - NCTSSoS.tr(x[1] * y[2]) -
            NCTSSoS.tr(x[2] * y[1]) + NCTSSoS.tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)

        @testset "Dense" begin
            config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracles.Ex_6_2_0_Dense_d1.opt atol = 1e-6
        end

        @testset "Term Sparsity Algorithms" begin
            for (name, algo) in [
                ("NoElimination", NoElimination()),
                ("MMD", MMD()),
                ("MaximalElimination", MaximalElimination())
            ]
                @testset "$name" begin
                    config = SolverConfig(
                        optimizer=SOLVER, order=1,
                        cs_algo=NoElimination(), ts_algo=algo
                    )
                    result = cs_nctssos(tpop, config)
                    @test result.objective ≈ oracles.Ex_6_2_0_Dense_d1.opt atol = 1e-6
                end
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

        @testset "Order 2 (tight bound)" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracles.Ex_6_2_1_Dense_d2.opt atol = 1e-5
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

        @testset "Dense" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracles.Ex_6_2_2_Dense_d2.opt atol = 1e-5
        end

        # Note: Use TS only (no CS) as NCTSSOS doesn't support correlative sparsity for trace polynomials.
        # Combining MF + MMD produces looser bounds due to smaller per-clique bases.
        @testset "Sparse (TS only)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracles.Ex_6_2_2_TS_d2.opt atol = 1e-4
        end
    end
end
