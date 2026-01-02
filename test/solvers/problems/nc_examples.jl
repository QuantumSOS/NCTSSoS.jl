# =============================================================================
# NC Polynomial Examples (from NCTSSOS paper)
# =============================================================================
# Consolidates classical NC polynomial optimization examples:
#   - Example 1: Unconstrained NC polynomial (3 variables)
#   - Example 2: Constrained NC polynomial (2 variables)
#   - Correlative Sparsity: Example with constraints (3 variables)
#   - CS TS n=10: Large-scale sparsity example
#
# These are standard benchmark problems from the NCTSSOS paper.
# Results verified against NCTSSOS oracles.
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__DIR__), "..", "setup.jl"))

# Load oracle values
include(joinpath(dirname(@__DIR__), "..", "oracles", "results", "example1_oracles.jl"))
include(joinpath(dirname(@__DIR__), "..", "oracles", "results", "corr_sparsity_oracles.jl"))
include(joinpath(dirname(@__DIR__), "..", "oracles", "results", "cs_ts_n10_oracles.jl"))

# Helper: flatten moment_matrix_sizes for comparison with oracle
flatten_sizes(sizes) = reduce(vcat, sizes)

@testset "NC Polynomial Examples" begin

    # =========================================================================
    # Example 1: Unconstrained NC polynomial
    # =========================================================================
    # Validated against NCTSSOS oracles: Example1_Dense_d2, Example1_TS_d2
    # =========================================================================
    @testset "Example 1 (unconstrained)" begin
        n = 3
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        f = 1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
            2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
            x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
            9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
            54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

        pop = polyopt(f, reg)

        @testset "Dense (order=2)" begin
            oracle = EXAMPLE1_ORACLES["Example1_Dense_d2"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config; dualize=false)
            # Mosek achieves 3.23e-8 error, COSMO needs 1e-5 tolerance
            @test result.objective ≈ oracle.opt atol = 1e-5
            @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
            @test result.n_unique_moment_matrix_elements == oracle.nuniq
        end

        @testset "Dense (SOS)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2
            )
            result = cs_nctssos(pop, config; dualize=true)
            @test result.objective ≈ 0.0 atol = 1e-6
        end

        @testset "Term Sparsity MMD (order=2)" begin
            oracle = EXAMPLE1_ORACLES["Example1_TS_d2"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config; dualize=false)
            @test result.objective ≈ oracle.opt atol = 1e-4
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
            @test result.n_unique_moment_matrix_elements == oracle.nuniq
        end
    end

    # =========================================================================
    # Example 2: Constrained NC polynomial
    # =========================================================================
    @testset "Example 2 (constrained)" begin
        n = 2
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
        g = 4.0 - x[1]^2 - x[2]^2
        h1 = x[1] * x[2] + x[2] * x[1] - 2.0

        pop = polyopt(f, reg; eq_constraints=[h1], ineq_constraints=[g])

        @testset "Dense (Moment)" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(pop, config; dualize=false)
            @test result.objective ≈ -1.0 atol = 1e-6
        end

        @testset "Dense (SOS)" begin
            config = SolverConfig(optimizer=SOLVER, order=2)
            result = cs_nctssos(pop, config; dualize=true)
            @test result.objective ≈ -1.0 atol = 1e-6
        end

        @testset "Term Sparsity (Moment)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config; dualize=false)
            @test result.objective ≈ -1.0 atol = 1e-6
        end

        @testset "Term Sparsity (SOS)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config; dualize=true)
            @test result.objective ≈ -1.0 atol = 1e-6
        end

        @testset "cs_nctssos_higher" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            result_higher = cs_nctssos_higher(pop, result, config)
            @test result.objective ≈ result_higher.objective atol = 1e-4
        end
    end

    # =========================================================================
    # Correlative Sparsity Example (n=3 with constraints)
    # =========================================================================
    # Validated against NCTSSOS oracles: CorrSparsity_CS_d3, CorrSparsity_TS_d3
    # =========================================================================
    @testset "Correlative Sparsity" begin
        n = 3
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        f = 1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
            2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
            x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
            9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
            54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

        cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0 / 3 for i = 1:n])
        pop = polyopt(f, reg; ineq_constraints=cons)

        @testset "CS MF (order=3)" begin
            oracle = CORR_SPARSITY_ORACLES["CorrSparsity_CS_d3"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config; dualize=false)
            @test result.objective ≈ oracle.opt atol = 1e-5
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
            @test result.n_unique_moment_matrix_elements == oracle.nuniq
        end

        @testset "CS MF (SOS)" begin
            oracle = CORR_SPARSITY_ORACLES["CorrSparsity_CS_d3"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=MF()
            )
            result = cs_nctssos(pop, config; dualize=true)
            @test result.objective ≈ oracle.opt atol = 1e-5
        end

        @testset "TS MMD (order=3)" begin
            oracle = CORR_SPARSITY_ORACLES["CorrSparsity_TS_d3"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=NoElimination(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config; dualize=false)
            result = cs_nctssos_higher(pop, result, config; dualize=false)
            @test result.objective ≈ oracle.opt atol = 1e-5
            # Block structure validated after higher iteration
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
            @test result.n_unique_moment_matrix_elements == oracle.nuniq
        end

        @testset "TS MMD (SOS)" begin
            oracle = CORR_SPARSITY_ORACLES["CorrSparsity_TS_d3"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config; dualize=true)
            @test result.objective ≈ oracle.opt atol = 1e-4
        end
    end

    # =========================================================================
    # CS TS n=10: Large-scale sparsity example (requires Mosek)
    # =========================================================================
    # Validated against NCTSSOS oracle: CS_TS_N10_CS_TS_d3
    # NOTE: Block sizes omitted from oracle (2982 blocks too large to store).
    #       We validate opt, nblocks, and nuniq instead.
    # =========================================================================
    if USE_LOCAL
        @testset "CS TS (n=10)" begin
            n = 10
            reg, (x,) = create_noncommutative_variables([("x", 1:n)])

            # Build polynomial using new API
            f = Polynomial{NonCommutativeAlgebra,UInt8,Float64}(
                Term{Monomial{NonCommutativeAlgebra,UInt8},Float64}[]
            )
            for i = 1:n
                jset = max(1, i - 5):min(n, i + 1)
                jset = setdiff(jset, i)
                f += (2.0 * x[i] + 5.0 * x[i]^3 + 1)^2
                f -= sum([
                    4.0 * x[i] * x[j] +
                    10.0 * x[i]^3 * x[j] +
                    2.0 * x[j] +
                    4.0 * x[i] * x[j]^2 +
                    10.0 * x[i]^3 * x[j]^2 +
                    2.0 * x[j]^2 for j in jset
                ])
                f += sum([
                    1.0 * x[j] * x[k] + 2.0 * x[j]^2 * x[k] + 1.0 * x[j]^2 * x[k]^2
                    for j in jset for k in jset
                ])
            end

            cons = vcat([(1.0 - x[i]^2) for i = 1:n], [(1.0 * x[i] - 1.0 / 3) for i = 1:n])
            pop = polyopt(f, reg; ineq_constraints=cons)

            oracle = CS_TS_N10_ORACLES["CS_TS_N10_CS_TS_d3"]

            @testset "Moment Method (order=3)" begin
                config = SolverConfig(
                    optimizer=SOLVER,
                    order=3,
                    cs_algo=MF(),
                    ts_algo=MMD()
                )
                result = cs_nctssos(pop, config; dualize=false)
                @test result.objective ≈ oracle.opt atol = 1e-3
                @test length(flatten_sizes(result.moment_matrix_sizes)) == oracle.nblocks
                @test result.n_unique_moment_matrix_elements == oracle.nuniq
            end

            @testset "SOS Dualization (order=3)" begin
                config = SolverConfig(
                    optimizer=SOLVER,
                    order=3,
                    cs_algo=MF(),
                    ts_algo=MMD()
                )
                result = cs_nctssos(pop, config; dualize=true)
                @test result.objective ≈ oracle.opt atol = 1e-3
            end
        end
    end

    # =========================================================================
    # README Examples
    # =========================================================================
    @testset "README Examples" begin
        @testset "Unconstrained" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            f = 1.0 + x[1]^4 + x[2]^4 + x[3]^4 +
                x[1] * x[2] + x[2] * x[1] + x[2] * x[3] + x[3] * x[2]

            pop = polyopt(f, reg)

            result_dense = cs_nctssos(pop, SolverConfig(optimizer=SOLVER))
            result_cs = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, cs_algo=MF()))
            result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, cs_algo=MF(), ts_algo=MMD()))

            @test result_dense.objective ≈ result_cs.objective atol = 1e-4
            @test result_cs.objective ≈ result_cs_ts.objective atol = 1e-4
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

            @test result_dense.objective ≈ result_cs.objective atol = 1e-4
            @test result_cs.objective ≈ result_cs_ts.objective atol = 1e-4

            result_cs_ts_higher = cs_nctssos_higher(
                pop,
                result_cs_ts,
                SolverConfig(optimizer=SOLVER, cs_algo=MF(), ts_algo=MMD())
            )
            @test result_dense.objective ≈ result_cs_ts_higher.objective atol = 1e-4
        end
    end
end
