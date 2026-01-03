# =============================================================================
# CHSH Bell Inequality Tests
# =============================================================================
# Tests CHSH Bell inequality optimization with various sparsity configurations.
# Results verified against NCTSSOS oracles.
#
# Expected optimal value: -2sqrt(2) ≈ -2.8284 (quantum bound)
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__DIR__), "..", "setup.jl"))

# Load oracle values
include(joinpath(dirname(@__DIR__), "..", "oracles", "results", "chsh_oracles.jl"))

# Helper: flatten moment_matrix_sizes for comparison with oracle
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
    # Dense (No Sparsity) - CHSH_Dense_d1
    # =========================================================================
    @testset "Dense (order=1)" begin
        pop, _ = create_chsh_problem()
        oracle = CHSH_ORACLES["CHSH_Dense_d1"]

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    # =========================================================================
    # Correlative Sparsity (MF) - CHSH_CS_d1
    # =========================================================================
    @testset "Correlative Sparsity MF (order=1)" begin
        pop, _ = create_chsh_problem()
        oracle = CHSH_ORACLES["CHSH_CS_d1"]

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    # =========================================================================
    # Term Sparsity (MMD) - CHSH_TS_d1
    # =========================================================================
    @testset "Term Sparsity MMD (order=1)" begin
        pop, _ = create_chsh_problem()
        oracle = CHSH_ORACLES["CHSH_TS_d1"]

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    # =========================================================================
    # Combined CS + TS (order=2) - CHSH_CS_TS_d2
    # =========================================================================
    # WARNING: Combined CS+TS does NOT converge to -2√2 for CHSH!
    # This gives opt=-4.0, a much looser bound than the quantum value -2.828.
    # 
    # Root cause: The cliques share only the identity element, so cross-clique
    # moment constraints are missing. This is a structural limitation of the
    # combined sparsity approach for problems with this specific structure.
    #
    # This is consistent with NCTSSOS behavior and is NOT a bug - it's an
    # inherent limitation of the CS+TS relaxation for certain problem types.
    # =========================================================================
    @testset "Combined CS+TS (order=2) - KNOWN LOOSE BOUND" begin
        pop, _ = create_chsh_problem()
        oracle = CHSH_ORACLES["CHSH_CS_TS_d2"]

        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MF(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)

        # Objective is -4.0, NOT -2.828 (this is expected!)
        @test result.objective ≈ oracle.opt atol = 1e-6
        # Block order may differ: compare sorted block sizes
        @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    # =========================================================================
    # State Polynomial Formulation
    # =========================================================================
    @testset "State Polynomial" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        sp = -ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
        spop = polyopt(sp * one(typeof(x[1])), reg)

        @testset "Dense" begin
            oracle = CHSH_ORACLES["CHSH_State_Dense_d1"]
            config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(spop, config)
            @test result.objective ≈ oracle.opt atol = 1e-5
        end

        @testset "Direct Moment vs SOS" begin
            oracle = CHSH_ORACLES["CHSH_State_Dense_d1"]
            config = SolverConfig(optimizer=SOLVER, order=1)
            result_mom = cs_nctssos(spop, config; dualize=false)
            result_sos = cs_nctssos(spop, config; dualize=true)
            @test result_mom.objective ≈ result_sos.objective atol = 1e-5
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
            oracle = CHSH_ORACLES["CHSH_Trace_Dense_d1"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracle.opt atol = 1e-4
        end

        @testset "Term Sparsity (MaximalElimination)" begin
            oracle = CHSH_ORACLES["CHSH_Trace_TS_d1"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                ts_algo=MaximalElimination()
            )
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ oracle.opt atol = 1e-5
        end
    end
end
