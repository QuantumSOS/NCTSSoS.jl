# NC Polynomial Optimization Benchmark Tests
# ==========================================
# Classical optimization benchmarks adapted for noncommutative polynomials.
# These test sparsity algorithms on structured problems.
#
# Reference: NCTSSOS paper benchmark problems
# All problems are unconstrained NC polynomial minimization.
#
# Oracle values available for Rosenbrock; others use theoretical global minima.
# Tolerances set to 1e-2 or 1e-3 to account for SDP relaxation approximation error.
#
# Benchmarks included:
#   1. Generalized Rosenbrock (n=6)  - Degree 4, Global min = 1
#   2. Broyden Banded (n=4)          - Degree 6, Global min = 0
#   3. Broyden Tridiagonal (n=6)     - Degree 4, Global min = 0
#   4. Chained Singular (n=8)        - Degree 4, Global min = 0
#   5. Chained Wood (n=8)            - Degree 4, Global min ≈ 1

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "..", "standalone_setup.jl"))

# Load oracle values
include(joinpath(dirname(@__FILE__), "..", "..", "oracles", "results", "rosenbrock_oracles.jl"))
include(joinpath(dirname(@__FILE__), "..", "..", "oracles", "results", "benchmarks_oracles.jl"))

# Helper: flatten moment_matrix_sizes for comparison with oracle
flatten_sizes(sizes) = reduce(vcat, sizes)

# All benchmark tests require Mosek for numerical stability
if USE_LOCAL

@testset "NC Polynomial Benchmarks" begin

    # =========================================================================
    # 1. Generalized Rosenbrock (n=6, d=2)
    # =========================================================================
    # f = n + Σᵢ₌₂ⁿ [100xᵢ₋₁⁴ - 200xᵢ₋₁²xᵢ - 2xᵢ + 101xᵢ²]
    # Global minimum: 1 at xᵢ = 0 for all i
    # Degree: 4, so d=2 is sufficient
    @testset "Generalized Rosenbrock (n=6)" begin
        n = 6
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        # Build the Rosenbrock objective
        # f = n + Σᵢ₌₂ⁿ [100xᵢ₋₁⁴ - 200xᵢ₋₁²xᵢ - 2xᵢ + 101xᵢ²]
        f = Float64(n) * one(typeof(x[1]))
        for i in 2:n
            f = f + 100.0 * x[i-1]^4 - 200.0 * x[i-1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
        end

        pop = polyopt(f, reg)
        expected = 1.0

        @testset "NoElimination" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-2
        end

        @testset "MF Correlative Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-2
        end

        @testset "MF + MMD Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-2
        end
    end

    # =========================================================================
    # 2. Broyden Banded (n=4, d=3)
    # =========================================================================
    # f = n + Σᵢ₌₁ⁿ [ 4xᵢ + 4xᵢ² + 10xᵢ³ + 20xᵢ⁴ + 25xᵢ⁶
    #                 + Σⱼ∈Jᵢ (-2xⱼ - 2xⱼ² - 4xᵢxⱼ - 4xᵢxⱼ² - 10xᵢ³xⱼ - 10xᵢ³xⱼ²)
    #                 + Σⱼ,ₖ∈Jᵢ (xⱼxₖ + 2xⱼ²xₖ + xⱼ²xₖ²) ]
    # where Jᵢ = {max(1,i-5), ..., min(n,i+1)} \ {i}
    # Global minimum: 0 at origin
    # Degree: 6, so d=3 is required
    @testset "Broyden Banded (n=4)" begin
        n = 4
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        # Build the Broyden Banded objective
        P = typeof(x[1] + x[2])
        f = Float64(n) * one(P)

        for i in 1:n
            jset = setdiff(max(1, i - 5):min(n, i + 1), i)

            # Single variable terms
            f = f + 4.0 * x[i] + 4.0 * x[i]^2 + 10.0 * x[i]^3 + 20.0 * x[i]^4 + 25.0 * x[i]^6

            # Cross terms with j in Jᵢ
            for j in jset
                f = f - 2.0 * x[j] - 2.0 * x[j]^2 - 4.0 * x[i] * x[j] -
                    4.0 * x[i] * x[j]^2 - 10.0 * x[i]^3 * x[j] - 10.0 * x[i]^3 * x[j]^2
            end

            # Double sum terms
            for j in jset
                for k in jset
                    f = f + x[j] * x[k] + 2.0 * x[j]^2 * x[k] + x[j]^2 * x[k]^2
                end
            end
        end

        pop = polyopt(f, reg)
        expected = 0.0

        @testset "NoElimination" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end

        @testset "MF Correlative Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end

        @testset "MF + MMD Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end
    end

    # =========================================================================
    # 3. Broyden Tridiagonal (n=6, d=2)
    # =========================================================================
    # f = n + Σᵢ₌₁ⁿ gᵢ(x) where:
    # g₁ = 5x₁² + 4x₁⁴ + 4x₂² - 12x₁³ - 12x₁x₂ + 6x₁ + 8x₁²x₂ - 4x₂
    # gᵢ = 5xᵢ² + 4xᵢ⁴ + xᵢ₋₁² + 4xᵢ₊₁² - 12xᵢ³ - 6xᵢ₋₁xᵢ - 12xᵢxᵢ₊₁
    #      + 6xᵢ + 4xᵢ₋₁xᵢ² + 8xᵢ²xᵢ₊₁ + 4xᵢ₋₁xᵢ₊₁ - 2xᵢ₋₁ - 4xᵢ₊₁  (i = 2,...,n-1)
    # gₙ = 5xₙ² + 4xₙ⁴ + xₙ₋₁² - 12xₙ³ - 6xₙ₋₁xₙ + 6xₙ + 4xₙ₋₁xₙ² - 2xₙ₋₁
    # Global minimum: 0 at origin
    # Degree: 4, so d=2 is sufficient
    @testset "Broyden Tridiagonal (n=6)" begin
        n = 6
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        # Build the Broyden Tridiagonal objective
        P = typeof(x[1] + x[2])
        f = Float64(n) * one(P)

        # First variable (i=1)
        f = f + 5.0 * x[1]^2 + 4.0 * x[1]^4 + 4.0 * x[2]^2 - 12.0 * x[1]^3 -
            12.0 * x[1] * x[2] + 6.0 * x[1] + 8.0 * x[1]^2 * x[2] - 4.0 * x[2]

        # Middle variables (i=2,...,n-1)
        for i in 2:(n-1)
            f = f + 5.0 * x[i]^2 + 4.0 * x[i]^4 + x[i-1]^2 + 4.0 * x[i+1]^2 -
                12.0 * x[i]^3 - 6.0 * x[i-1] * x[i] - 12.0 * x[i] * x[i+1] +
                6.0 * x[i] + 4.0 * x[i-1] * x[i]^2 + 8.0 * x[i]^2 * x[i+1] +
                4.0 * x[i-1] * x[i+1] - 2.0 * x[i-1] - 4.0 * x[i+1]
        end

        # Last variable (i=n)
        f = f + 5.0 * x[n]^2 + 4.0 * x[n]^4 + x[n-1]^2 - 12.0 * x[n]^3 -
            6.0 * x[n-1] * x[n] + 6.0 * x[n] + 4.0 * x[n-1] * x[n]^2 - 2.0 * x[n-1]

        pop = polyopt(f, reg)
        expected = 0.0

        @testset "NoElimination" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end

        @testset "MF Correlative Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end

        @testset "MF + MMD Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end
    end

    # =========================================================================
    # 4. Chained Singular (n=8, d=2)
    # =========================================================================
    # f = Σᵢ₌₁,₃,₅,...,ₙ₋₃ [ (xᵢ + 10xᵢ₊₁)² + 5(xᵢ₊₂ - xᵢ₊₃)²
    #                        + (xᵢ₊₁ - 2xᵢ₊₂)⁴ + 10(xᵢ - 10xᵢ₊₃)⁴ ]
    # In NC form: (a+b)² = a² + ab + ba + b², (a+b)⁴ expanded accordingly
    # Global minimum: 0 at origin
    # Degree: 4, so d=2 is sufficient
    @testset "Chained Singular (n=8)" begin
        n = 8
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        P = typeof(x[1] + x[2])
        f = zero(P)

        # Helper for NC square: (a + b)² = a² + ab + ba + b²
        nc_square(a, b) = a^2 + a * b + b * a + b^2

        # Helper for NC fourth power: (a + b)⁴ expanded
        # (a+b)⁴ = ((a+b)²)² but we need to be careful with NC
        # For simplicity, use: (a+b)^4 = a^4 + a^3*b + a^2*b*a + a^2*b^2 + a*b*a^2 + ...
        # Actually, let's directly expand using polynomial arithmetic
        function nc_fourth_power(a, b)
            ab = a + b
            ab2 = ab * ab  # (a+b)²
            return ab2 * ab2  # (a+b)⁴
        end

        for i in 1:2:(n-3)
            # (xᵢ + 10xᵢ₊₁)²
            f = f + nc_square(x[i], 10.0 * x[i+1])

            # 5(xᵢ₊₂ - xᵢ₊₃)²
            f = f + 5.0 * nc_square(x[i+2], -1.0 * x[i+3])

            # (xᵢ₊₁ - 2xᵢ₊₂)⁴
            f = f + nc_fourth_power(x[i+1], -2.0 * x[i+2])

            # 10(xᵢ - 10xᵢ₊₃)⁴
            f = f + 10.0 * nc_fourth_power(x[i], -10.0 * x[i+3])
        end

        pop = polyopt(f, reg)
        expected = 0.0

        @testset "NoElimination" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end

        @testset "MF Correlative Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end

        @testset "MF + MMD Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-3
        end
    end

    # =========================================================================
    # 5. Chained Wood (n=8, d=2)
    # =========================================================================
    # f = (21n - 41) + Σᵢ₌₁,₃,₅,...,ₙ₋₃ gᵢ(x) where each block gᵢ is:
    #   gᵢ = -2xᵢ + xᵢ² + 100xᵢ⁴ - 200xᵢ²xᵢ₊₁ - 40xᵢ₊₁ + 110.1xᵢ₊₁² + 19.8xᵢ₊₁xᵢ₊₃
    #        - 2xᵢ₊₂ + xᵢ₊₂² + 90xᵢ₊₂⁴ - 180xᵢ₊₂²xᵢ₊₃ - 40xᵢ₊₃ + 100.1xᵢ₊₃²
    # Global minimum: ≈1
    # Degree: 4, so d=2 is sufficient
    @testset "Chained Wood (n=8)" begin
        n = 8
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        P = typeof(x[1] + x[2])
        # Constant term
        f = Float64(21 * n - 41) * one(P)

        for i in 1:2:(n-3)
            f = f + (-2.0) * x[i] + x[i]^2 + 100.0 * x[i]^4 - 200.0 * x[i]^2 * x[i+1] +
                (-40.0) * x[i+1] + 110.1 * x[i+1]^2 + 19.8 * x[i+1] * x[i+3] +
                (-2.0) * x[i+2] + x[i+2]^2 + 90.0 * x[i+2]^4 - 180.0 * x[i+2]^2 * x[i+3] +
                (-40.0) * x[i+3] + 100.1 * x[i+3]^2
        end

        pop = polyopt(f, reg)
        expected = 1.0

        @testset "NoElimination" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-2
        end

        @testset "MF Correlative Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-2
        end

        @testset "MF + MMD Sparsity" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ expected atol = 1e-2
        end
    end

    # =========================================================================
    # Sparsity Comparison Test
    # =========================================================================
    # Verifies that sparsity methods produce consistent results
    @testset "Sparsity Method Consistency" begin
        # Use Rosenbrock as a representative test case
        n = 4
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        f = Float64(n) * one(typeof(x[1]))
        for i in 2:n
            f = f + 100.0 * x[i-1]^4 - 200.0 * x[i-1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
        end

        pop = polyopt(f, reg)

        # Collect results from different configurations
        configs = [
            ("Dense", NoElimination(), NoElimination()),
            ("CS only", MF(), NoElimination()),
            ("TS only", NoElimination(), MMD()),
            ("CS + TS", MF(), MMD()),
        ]

        results = Dict{String,Float64}()
        for (name, cs_algo, ts_algo) in configs
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=cs_algo,
                ts_algo=ts_algo
            )
            result = cs_nctssos(pop, config)
            results[name] = result.objective
        end

        # All methods should give similar results (within tolerance)
        baseline = results["Dense"]
        for (name, val) in results
            @test isapprox(val, baseline; atol=0.1)
        end

        # All should be close to the expected minimum
        @test baseline ≈ 1.0 atol = 1e-2
    end

end  # @testset "NC Polynomial Benchmarks"

end  # USE_LOCAL
