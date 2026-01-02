# =============================================================================
# CHSH Bell Inequality Tests
# =============================================================================
# Consolidates all CHSH-related polynomial optimization tests:
#   - Standard CHSH with moment/SOS methods
#   - Sparsity variants (dense, correlative, term)
#   - State polynomial formulation
#   - Trace polynomial formulation
#
# Expected optimal value: -2sqrt(2) ≈ -2.8284 (quantum bound)
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__DIR__), "..", "setup.jl"))

const CHSH_EXPECTED = -2.8284271321623193

@testset "CHSH Bell Inequality" begin

    # =========================================================================
    # Standard NC Polynomial Formulation
    # =========================================================================
    @testset "Unipotent Variables" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        f = x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(f, reg)

        @testset "Moment Method" begin
            config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(pop, config; dualize=false)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-6
        end

        @testset "SOS Dualization" begin
            config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(pop, config; dualize=true)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-6
        end
    end

    # =========================================================================
    # Sparsity Variants
    # =========================================================================
    @testset "Sparsity Methods" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(-f, reg)

        @testset "Dense (No Sparsity)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-4
        end

        @testset "Correlative Sparsity (MF)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-4
        end

        @testset "Term Sparsity (MMD)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-4
        end

    end

    # =========================================================================
    # State Polynomial Formulation
    # =========================================================================
    @testset "State Polynomial" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        sp = -ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
        spop = polyopt(sp * one(typeof(x[1])), reg)

        @testset "Dense" begin
            config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(spop, config)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-5
        end

        @testset "Direct Moment vs SOS" begin
            config = SolverConfig(optimizer=SOLVER, order=1)
            result_mom = cs_nctssos(spop, config; dualize=false)
            result_sos = cs_nctssos(spop, config; dualize=true)
            @test result_mom.objective ≈ result_sos.objective atol = 1e-5
        end

        @testset "Term Sparsity (MMD)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(spop, config)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-5
        end
    end

    # =========================================================================
    # Trace Polynomial Formulation
    # =========================================================================
    @testset "Trace Polynomial" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]

        p = -1.0 * NCTSSoS.tr(x[1] * y[1]) - NCTSSoS.tr(x[1] * y[2]) - NCTSSoS.tr(x[2] * y[1]) + NCTSSoS.tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)

        @testset "Dense" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-4
        end

        @testset "Term Sparsity (MaximalElimination)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=1,
                ts_algo=MaximalElimination()
            )
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ CHSH_EXPECTED atol = 1e-5
        end
    end
end
