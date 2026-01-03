# =============================================================================
# CHSH Bell Inequality Tests
# =============================================================================
# Tests CHSH Bell inequality optimization with various sparsity configurations.
# Results verified against NCTSSOS.
#
# Expected optimal value: -2sqrt(2) ≈ -2.8284 (quantum bound)
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
end

# =============================================================================
# Expected values from NCTSSOS
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization → opt ≈ -2√2)
#   sides = moment matrix block sizes
#   nuniq = unique moment indices (affine constraints)
# =============================================================================
const EXPECTED_CHSH = (
    # NC Polynomial Formulation
    Dense_d1    = (opt=-2.8284271321623193, sides=[5], nuniq=11),
    CS_d1       = (opt=-2.8284271247170496, sides=[4, 4], nuniq=10),
    TS_d1       = (opt=-2.8284271247321175, sides=[3, 3, 1], nuniq=6),
    CS_TS_d2    = (opt=-3.999999999803662, sides=[3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 2, 2], nuniq=5),
    # State Polynomial Formulation
    State_Dense = (opt=-2.828427124746234, sides=[9], nuniq=21),
    State_TS    = (opt=-2.828427124746232, sides=[9], nuniq=21),
    # Trace Polynomial Formulation
    Trace_Dense = (opt=-2.828427124746234, sides=[9], nuniq=21),
    Trace_TS    = (opt=-2.828427124746232, sides=[9], nuniq=21),
)

# Helper: flatten moment_matrix_sizes for comparison
flatten_sizes(sizes) = reduce(vcat, sizes)

@testset "CHSH Bell Inequality" begin

    # =========================================================================
    # Problem Setup Helper
    # =========================================================================
    function create_chsh_problem()
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        # CHSH: maximize x1*y1 + x1*y2 + x2*y1 - x2*y2
        # We minimize -f, so objective is negated CHSH
        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(-f, reg)
        return pop, reg
    end

    # =========================================================================
    # Dense (No Sparsity)
    # =========================================================================
    @testset "Dense (order=1)" begin
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ EXPECTED_CHSH.Dense_d1.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == EXPECTED_CHSH.Dense_d1.sides
        @test result.n_unique_moment_matrix_elements == EXPECTED_CHSH.Dense_d1.nuniq
    end

    # =========================================================================
    # Correlative Sparsity (MF)
    # =========================================================================
    @testset "Correlative Sparsity MF (order=1)" begin
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ EXPECTED_CHSH.CS_d1.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == EXPECTED_CHSH.CS_d1.sides
        @test result.n_unique_moment_matrix_elements == EXPECTED_CHSH.CS_d1.nuniq
    end

    # =========================================================================
    # Term Sparsity (MMD)
    # =========================================================================
    @testset "Term Sparsity MMD (order=1)" begin
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ EXPECTED_CHSH.TS_d1.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == EXPECTED_CHSH.TS_d1.sides
        @test result.n_unique_moment_matrix_elements == EXPECTED_CHSH.TS_d1.nuniq
    end

    # =========================================================================
    # Combined CS + TS (order=2)
    # WARNING: Combined CS+TS does NOT converge to -2√2 for CHSH!
    # This gives opt=-4.0, a much looser bound than the quantum value -2.828.
    # This is consistent with NCTSSOS behavior - it's an inherent limitation.
    # =========================================================================
    @testset "Combined CS+TS (order=2) - KNOWN LOOSE BOUND" begin
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MF(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ EXPECTED_CHSH.CS_TS_d2.opt atol = 1e-6
        # Block sizes may differ but total count matches
        @test length(flatten_sizes(result.moment_matrix_sizes)) == length(EXPECTED_CHSH.CS_TS_d2.sides)
        @test result.n_unique_moment_matrix_elements == EXPECTED_CHSH.CS_TS_d2.nuniq
    end

    # =========================================================================
    # State Polynomial Formulation
    # =========================================================================
    @testset "State Polynomial" begin
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
            @test result.objective ≈ EXPECTED_CHSH.State_Dense.opt atol = 1e-5
            @test flatten_sizes(result.moment_matrix_sizes) == EXPECTED_CHSH.State_Dense.sides
            @test result.n_unique_moment_matrix_elements == EXPECTED_CHSH.State_Dense.nuniq
        end

        @testset "Term Sparsity MMD" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(spop, config)
            @test result.objective ≈ EXPECTED_CHSH.State_TS.opt atol = 1e-5
            @test flatten_sizes(result.moment_matrix_sizes) == EXPECTED_CHSH.State_TS.sides
            @test result.n_unique_moment_matrix_elements == EXPECTED_CHSH.State_TS.nuniq
        end

        @testset "Direct Moment vs SOS" begin
            config = SolverConfig(optimizer=SOLVER, order=1)
            result_mom = cs_nctssos(spop, config; dualize=false)
            result_sos = cs_nctssos(spop, config; dualize=true)
            @test result_mom.objective ≈ result_sos.objective atol = 1e-5
            @test flatten_sizes(result_mom.moment_matrix_sizes) == EXPECTED_CHSH.State_Dense.sides
            @test flatten_sizes(result_sos.moment_matrix_sizes) == EXPECTED_CHSH.State_Dense.sides
            @test result_mom.n_unique_moment_matrix_elements == EXPECTED_CHSH.State_Dense.nuniq
            @test result_sos.n_unique_moment_matrix_elements == EXPECTED_CHSH.State_Dense.nuniq
        end
    end

    # =========================================================================
    # Trace Polynomial Formulation
    # =========================================================================
    @testset "Trace Polynomial" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]

        p = -1.0 * NCTSSoS.tr(x[1] * y[1]) - NCTSSoS.tr(x[1] * y[2]) -
            NCTSSoS.tr(x[2] * y[1]) + NCTSSoS.tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)

        @testset "Dense" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ EXPECTED_CHSH.Trace_Dense.opt atol = 1e-5
            @test flatten_sizes(result.moment_matrix_sizes) == EXPECTED_CHSH.Trace_Dense.sides
            @test result.n_unique_moment_matrix_elements == EXPECTED_CHSH.Trace_Dense.nuniq
        end

        @testset "Term Sparsity MMD" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ EXPECTED_CHSH.Trace_TS.opt atol = 1e-5
            @test flatten_sizes(result.moment_matrix_sizes) == EXPECTED_CHSH.Trace_TS.sides
            @test result.n_unique_moment_matrix_elements == EXPECTED_CHSH.Trace_TS.nuniq
        end
    end
end
