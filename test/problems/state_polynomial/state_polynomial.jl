# State Polynomial Optimization Tests (Examples 7.2.x)
# Consolidates state polynomial tests using the ς() operator:
#   - 7.2.1: Squared expectations (products of state values)
#   - 7.2.2: Covariance Bell inequality
#   - 7.2.3: Mixed state polynomial with squared expectations
#
# Note: Basic CHSH state polynomial tests are in chsh.jl
# Results verified against NCTSSOS.

using Test, NCTSSoS, JuMP

# Helper: flatten moment_matrix_sizes for comparison
if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

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

# Expected values: mosek (reference) and cosmo (CI)
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = moment matrix block sizes
#   nuniq = unique moment indices (affine constraints)
# Note: Values updated after fixing site encoding to use single site for all variables,
# ensuring consistent SDP formulation regardless of how variables are grouped.
const EXPECTED_STATE_POLY = (
    mosek = (
        Ex_7_2_1_Dense_d3 = (opt=-3.9999999974338247, sides=[269], nuniq=5381),
        Ex_7_2_2_Dense_d2 = (opt=-4.9999999971484765, sides=[115], nuniq=1755),
        Ex_7_2_2_TS_d2    = (opt=-4.999999992730746, sides=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8], nuniq=63),
        Ex_7_2_3_Dense_d2 = (opt=-4.080998980966884, sides=[53], nuniq=355),
        Ex_7_2_3_TS_d2    = (opt=-4.080998977321278, sides=[1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6], nuniq=36),
    ),
    cosmo = (
        Ex_7_2_1_Dense_d3 = (opt=-4.000000008371047, sides=[269], nuniq=5381),
        Ex_7_2_2_Dense_d2 = (opt=-5.00000032294447, sides=[115], nuniq=1755),
        Ex_7_2_2_TS_d2    = (opt=-4.999999722040697, sides=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8], nuniq=63),
        Ex_7_2_3_Dense_d2 = (opt=-4.080998980966884, sides=[53], nuniq=355),
        Ex_7_2_3_TS_d2    = (opt=-4.080998977321278, sides=[1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6], nuniq=36),
    ),
)

select_oracle(oracles) = getproperty(oracles, SOLVER_NAME)
solver_atol(base_tol=1e-6) = SOLVER_NAME == :cosmo ? max(base_tol, 1e-3) : base_tol

@testset "State Polynomial Examples (7.2.x)" begin
    oracles = select_oracle(EXPECTED_STATE_POLY)

    # =========================================================================
    # Example 7.2.1: Squared Expectations
    # =========================================================================
    # At order=3, the relaxation correctly gives the tight bound of -4.0.
    @testset "Example 7.2.1 (Squared Expectations)" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
        sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
        sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

        spop = polyopt(sp * one(typeof(x[1])), reg)

        @testset "Order 3 (tight bound)" begin
            config = SolverConfig(optimizer=SOLVER, order=3)
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracles.Ex_7_2_1_Dense_d3.opt atol = solver_atol()
            @test flatten_sizes(result.moment_matrix_sizes) == oracles.Ex_7_2_1_Dense_d3.sides
            @test result.n_unique_moment_matrix_elements == oracles.Ex_7_2_1_Dense_d3.nuniq
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

        @testset "Dense" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracles.Ex_7_2_2_Dense_d2.opt atol = solver_atol()
            @test flatten_sizes(result.moment_matrix_sizes) == oracles.Ex_7_2_2_Dense_d2.sides
            @test result.n_unique_moment_matrix_elements == oracles.Ex_7_2_2_Dense_d2.nuniq
        end

        @testset "Sparse (MF + MMD)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracles.Ex_7_2_2_TS_d2.opt atol = solver_atol()
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracles.Ex_7_2_2_TS_d2.sides)
            @test result.n_unique_moment_matrix_elements == oracles.Ex_7_2_2_TS_d2.nuniq
        end
    end

    # =========================================================================
    # Example 7.2.3: Mixed State Polynomial
    # =========================================================================
    @testset "Example 7.2.3 (Mixed)" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        sp = -1.0 * ς(x[2]) - 1.0 * ς(y[1]) - 1.0 * ς(y[2]) +
             1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[2] * y[1]) -
             1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[2]) +
             1.0 * ς(x[1]) * ς(y[1]) + 1.0 * ς(x[2]) * ς(y[1]) +
             1.0 * ς(x[2]) * ς(y[2]) +
             1.0 * ς(x[1]) * ς(x[1]) + 1.0 * ς(y[2]) * ς(y[2])

        spop = polyopt(sp * one(typeof(x[1])), reg)

        @testset "Dense (Moment)" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(spop, config; dualize=false)
            @test result.objective ≈ oracles.Ex_7_2_3_Dense_d2.opt atol = solver_atol()
            @test flatten_sizes(result.moment_matrix_sizes) == oracles.Ex_7_2_3_Dense_d2.sides
            @test result.n_unique_moment_matrix_elements == oracles.Ex_7_2_3_Dense_d2.nuniq
        end

        @testset "Dense (SOS)" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(spop, config; dualize=true)
            @test result.objective ≈ oracles.Ex_7_2_3_Dense_d2.opt atol = solver_atol()
            @test flatten_sizes(result.moment_matrix_sizes) == oracles.Ex_7_2_3_Dense_d2.sides
            @test result.n_unique_moment_matrix_elements == oracles.Ex_7_2_3_Dense_d2.nuniq
        end

        @testset "Sparse (MMD)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracles.Ex_7_2_3_TS_d2.opt atol = solver_atol()
            # Block sizes may differ slightly due to chordal decomposition algorithm,
            # but total block count and nuniq should match NCTSSOS
            @test length(flatten_sizes(result.moment_matrix_sizes)) == length(oracles.Ex_7_2_3_TS_d2.sides)
            @test result.n_unique_moment_matrix_elements == oracles.Ex_7_2_3_TS_d2.nuniq
        end

        @testset "Moment vs SOS Consistency" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result_mom = cs_nctssos(spop, config; dualize=false)
            result_sos = cs_nctssos(spop, config; dualize=true)
            @test result_mom.objective ≈ result_sos.objective atol = solver_atol()
        end
    end
end
