# =============================================================================
# I_3322 Bell Inequality Tests
# =============================================================================
# Tests I_3322 Bell inequality optimization with various sparsity configurations.
# Results verified against NCTSSOS.
#
# The I_3322 inequality uses 3 measurements per party with 2 outcomes each.
# Expected optimal value: ≈ -0.2509 (quantum bound)
#
# NOTE: For Bell inequalities with partition parameter, correlative sparsity (CS)
# is NOT applicable - partition already exploits bipartite structure.
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
# Expected values from NCTSSOS (mosek) and COSMO (CI)
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization → opt ≈ -0.2509)
#   sides = moment matrix block sizes
#   nuniq = unique moment indices (affine constraints)
# =============================================================================
const EXPECTED_I3322 = (
    mosek = (
        Dense_d2 = (opt=-0.25093979763955865, sides=[28], nuniq=154),
        Dense_d3 = (opt=-0.25087575215587177, sides=[88], nuniq=868),
        TS_d3    = (opt=-0.25147090316786014, sides=[13, 12, 12, 12, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2], nuniq=121),
    ),
    cosmo = (
        Dense_d2 = (opt=-0.2509094376715752, sides=[28], nuniq=154),
        Dense_d3 = (opt=-0.2507493348554742, sides=[88], nuniq=868),
        TS_d3    = (opt=-0.25147243754429394, sides=[2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 12, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 2, 4, 4, 4, 6, 9, 12, 13, 12], nuniq=121),
    ),
)

# Helper: flatten moment_matrix_sizes for comparison
flatten_sizes(sizes) = reduce(vcat, sizes)
select_oracle(oracles) = getproperty(oracles, SOLVER_NAME)
solver_atol(base_tol=1e-6) = SOLVER_NAME == :cosmo ? max(base_tol, 1e-3) : base_tol

@testset "I_3322 Bell Inequality" begin
    oracles = select_oracle(EXPECTED_I3322)

    # =========================================================================
    # Problem Setup Helper
    # =========================================================================
    function create_i3322_problem()
        reg, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) +
            x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) -
            x[1] - 2.0 * y[1] - y[2]
        pop = polyopt(-f, reg)
        return pop, reg
    end

    # =========================================================================
    # Dense (No Sparsity) - Order 2
    # =========================================================================
    @testset "Dense (order=2)" begin
        pop, _ = create_i3322_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracles.Dense_d2.opt atol = solver_atol()
        @test flatten_sizes(result.moment_matrix_sizes) == oracles.Dense_d2.sides
        @test result.n_unique_moment_matrix_elements == oracles.Dense_d2.nuniq
    end

    # =========================================================================
    # Dense (No Sparsity) - Order 3 (achieves quantum bound)
    # =========================================================================
    @testset "Dense (order=3)" begin
        pop, _ = create_i3322_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracles.Dense_d3.opt atol = solver_atol()
        @test flatten_sizes(result.moment_matrix_sizes) == oracles.Dense_d3.sides
        @test result.n_unique_moment_matrix_elements == oracles.Dense_d3.nuniq
    end

    # =========================================================================
    # Term Sparsity (order=3)
    # =========================================================================
    @testset "Term Sparsity MMD (order=3)" begin
        pop, _ = create_i3322_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracles.TS_d3.opt atol = solver_atol()
        @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracles.TS_d3.sides)
        @test result.n_unique_moment_matrix_elements == oracles.TS_d3.nuniq
    end
end
