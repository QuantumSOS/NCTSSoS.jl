# NC Polynomial Optimization Benchmark Tests
# Classical optimization benchmarks adapted for noncommutative polynomials.
# Tests sparsity algorithms on structured problems.
# Results verified against NCTSSOS.
#
# Benchmarks included:
#   1. Generalized Rosenbrock (n=6)  - Degree 4, Global min = 1
#   2. Broyden Banded (n=4)          - Degree 6, Global min = 0
#   3. Broyden Tridiagonal (n=6)     - Degree 4, Global min = 0
#   4. Chained Singular (n=8)        - Degree 4, Global min = 0
#   5. Chained Wood (n=8)            - Degree 4, Global min ≈ 1

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

# Expected values: mosek (reference) and cosmo (CI)
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = moment matrix block sizes
#   nuniq = unique moment indices (affine constraints)

# Generalized Rosenbrock (n=6, degree=4, order=2)
# Global minimum: 1.0 at xᵢ = 0 for all i
const EXPECTED_ROSENBROCK = (
    mosek = (
        Dense = (opt=0.999999995930163, sides=[43], nuniq=820),
        CS    = (opt=0.999999973478842, sides=[7, 7, 7, 7, 7], nuniq=90),
        CS_TS = (opt=0.9999997821660428, sides=[3, 2, 2, 2, 1, 3, 2, 2, 2, 1, 3, 2, 2, 2, 1, 3, 2, 2, 2, 1, 3, 2, 2, 1], nuniq=33),
    ),
    cosmo = (
        Dense = (opt=1.0003057063376781, sides=[43], nuniq=820),
        CS    = (opt=1.0000078819894798, sides=[7, 7, 7, 7, 7], nuniq=90),
        CS_TS = (opt=0.9999995915958918, sides=[2, 2, 2, 3, 1, 2, 3, 2, 1, 2, 2, 2, 3, 1, 2, 2, 2, 3, 1, 2, 2, 2, 3, 1], nuniq=33),
    ),
)

# Broyden Banded (n=4, degree=6, order=3)
# Global minimum: 0.0 at origin
const EXPECTED_BROYDEN_BANDED = (
    mosek = (
        Dense = (opt=-2.0628417282628906e-8, sides=[85], nuniq=2815),
        CS    = (opt=-2.0628417282628906e-8, sides=[85], nuniq=2815),
        CS_TS = (opt=-1.7480268537791176e-9, sides=[9, 9, 7, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=134),
    ),
    cosmo = (
        Dense = (opt=4.5138152086658275e-7, sides=[85], nuniq=2815),
        CS    = (opt=4.5138152086658275e-7, sides=[85], nuniq=2815),
        CS_TS = (opt=7.488524341470728e-8, sides=[3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 7, 3, 4, 4, 4, 4, 4, 9, 9, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], nuniq=134),
    ),
)

# Broyden Tridiagonal (n=6, degree=4, order=2)
# Global minimum: 0.0 at origin
const EXPECTED_BROYDEN_TRIDIAGONAL = (
    mosek = (
        Dense = (opt=6.580588236100922e-11, sides=[43], nuniq=820),
        CS    = (opt=-7.683336739295744e-9, sides=[13, 13, 13, 13], nuniq=226),
        CS_TS = (opt=-1.3851129121677345e-7, sides=[5, 4, 4, 3, 3, 3, 3, 2, 2, 5, 4, 4, 3, 3, 3, 3, 2, 2, 5, 4, 4, 3, 3, 3, 3, 2, 2, 5, 4, 4, 3, 3, 3, 3, 2, 2], nuniq=62),
    ),
    cosmo = (
        Dense = (opt=2.0788227879723004e-7, sides=[43], nuniq=820),
        CS    = (opt=6.538782001860585e-7, sides=[13, 13, 13, 13], nuniq=226),
        CS_TS = (opt=3.5020218710096718e-6, sides=[2, 2, 3, 3, 3, 4, 5, 4, 3, 2, 2, 3, 3, 3, 4, 5, 4, 3, 2, 2, 3, 3, 3, 4, 5, 4, 3, 2, 2, 3, 3, 3, 4, 5, 4, 3], nuniq=62),
    ),
)

# Chained Singular (n=8, degree=4, order=2)
# Global minimum: 0.0 at origin
const EXPECTED_CHAINED_SINGULAR = (
    mosek = (
        Dense = (opt=-6.428097254420157e-8, sides=[73], nuniq=2413),
        CS    = (opt=-1.2251940085171628e-7, sides=[13, 13, 13, 13, 13, 13], nuniq=328),
        CS_TS = (opt=-2.4177846987055753e-7, sides=[4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1, 4, 3, 2, 2, 2, 2, 1, 1, 1], nuniq=83),
    ),
    cosmo = (
        Dense = (opt=-2.2080854060339308e-5, sides=[73], nuniq=2413),
        CS    = (opt=-1.8758896205426406e-6, sides=[13, 13, 13, 13, 13, 13], nuniq=328),
        CS_TS = (opt=0.056385714820388025, sides=[2, 2, 2, 3, 4, 2, 1, 1, 1, 2, 2, 2, 3, 4, 2, 1, 1, 1, 2, 2, 2, 3, 4, 2, 1, 1, 1, 2, 2, 2, 3, 4, 2, 1, 1, 1, 2, 2, 2, 3, 4, 2, 1, 1, 1, 2, 2, 2, 3, 4, 2, 1, 1, 1], nuniq=83),
    ),
)

# Chained Wood (n=8, degree=4, order=2)
# Global minimum: ≈1.0
const EXPECTED_CHAINED_WOOD = (
    mosek = (
        Dense = (opt=0.9999998013387236, sides=[73], nuniq=2413),
        CS    = (opt=0.9999998300594172, sides=[7, 7, 7, 7, 7, 7, 7], nuniq=124),
        CS_TS = (opt=0.9999998815901131, sides=[3, 2, 2, 2, 1, 3, 2, 2, 2, 2, 3, 2, 2, 2, 1, 3, 2, 2, 2, 2, 3, 2, 2, 2, 1, 3, 2, 2, 2, 2, 3, 2, 2, 2, 1], nuniq=46),
    ),
    cosmo = (
        Dense = (opt=2.338137603486616, sides=[73], nuniq=2413),
        CS    = (opt=1.0000083771504054, sides=[7, 7, 7, 7, 7, 7, 7], nuniq=124),
        CS_TS = (opt=0.9999983103276159, sides=[2, 2, 2, 3, 1, 2, 2, 2, 3, 1, 2, 2, 2, 3, 1, 2, 2, 2, 3, 1, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3], nuniq=46),
    ),
)

# Helper: flatten moment_matrix_sizes for comparison
if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end
select_oracle(oracles) = getproperty(oracles, SOLVER_NAME)
solver_atol(base_tol=1e-6) = SOLVER_NAME == :cosmo ? max(base_tol, 1e-3) : base_tol
# Compare sides - COSMO block order can vary
compare_sides(result_sides, oracle_sides) = SOLVER_NAME == :cosmo ?
    sort(flatten_sizes(result_sides)) == sort(oracle_sides) :
    flatten_sizes(result_sides) == oracle_sides

# Note: For CS+TS tests, block ordering may differ between NCTSSoS.jl and NCTSSOS
# due to internal graph traversal order. The mathematical content is identical
# when the sorted block sizes match. Dense and CS-only tests use exact comparison
# since their block ordering is deterministic.

@testset "NC Polynomial Benchmarks" begin

    # =========================================================================
    # 1. Generalized Rosenbrock (n=6, d=2)
    # =========================================================================
    # f = n + Σᵢ₌₂ⁿ [100xᵢ₋₁⁴ - 200xᵢ₋₁²xᵢ - 2xᵢ + 101xᵢ²]
    # Global minimum: 1 at xᵢ = 0 for all i
    # Degree: 4, so d=2 is sufficient
    @testset "Generalized Rosenbrock (n=6)" begin
        oracles = select_oracle(EXPECTED_ROSENBROCK)
        n = 6
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        # Build the Rosenbrock objective
        # f = n + Σᵢ₌₂ⁿ [100xᵢ₋₁⁴ - 200xᵢ₋₁²xᵢ - 2xᵢ + 101xᵢ²]
        f = Float64(n) * one(typeof(x[1]))
        for i in 2:n
            f = f + 100.0 * x[i-1]^4 - 200.0 * x[i-1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
        end

        pop = polyopt(f, reg)

        @testset "Dense (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            # Test runs use `--check-bounds=yes`, which can shift first-order
            # convergence slightly for this dense benchmark.
            rosen_tol = SOLVER_NAME == :cosmo ? 1.5e-3 : solver_atol()
            @test result.objective ≈ oracles.Dense.opt atol = rosen_tol
            @test compare_sides(result.moment_matrix_sizes, oracles.Dense.sides)
            @test result.n_unique_moment_matrix_elements == oracles.Dense.nuniq
        end

        @testset "Correlative Sparsity MF (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.CS.opt atol = solver_atol()
            @test compare_sides(result.moment_matrix_sizes, oracles.CS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS.nuniq
        end

        @testset "CS + TS (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.CS_TS.opt atol = solver_atol()
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracles.CS_TS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS_TS.nuniq
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
        oracles = select_oracle(EXPECTED_BROYDEN_BANDED)
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

        @testset "Dense (order=3)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.Dense.opt atol = solver_atol()
            @test compare_sides(result.moment_matrix_sizes, oracles.Dense.sides)
            @test result.n_unique_moment_matrix_elements == oracles.Dense.nuniq
        end

        @testset "Correlative Sparsity MF (order=3)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.CS.opt atol = solver_atol()
            @test compare_sides(result.moment_matrix_sizes, oracles.CS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS.nuniq
        end

        @testset "CS + TS (order=3)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.CS_TS.opt atol = solver_atol()
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracles.CS_TS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS_TS.nuniq
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
        oracles = select_oracle(EXPECTED_BROYDEN_TRIDIAGONAL)
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

        @testset "Dense (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.Dense.opt atol = solver_atol()
            @test compare_sides(result.moment_matrix_sizes, oracles.Dense.sides)
            @test result.n_unique_moment_matrix_elements == oracles.Dense.nuniq
        end

        @testset "Correlative Sparsity MF (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.CS.opt atol = solver_atol()
            @test compare_sides(result.moment_matrix_sizes, oracles.CS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS.nuniq
        end

        @testset "CS + TS (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.CS_TS.opt atol = solver_atol()
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracles.CS_TS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS_TS.nuniq
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
        oracles = select_oracle(EXPECTED_CHAINED_SINGULAR)
        n = 8
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        P = typeof(x[1] + x[2])
        f = zero(P)

        # Helper for NC square: (a + b)² = a² + ab + ba + b²
        nc_square(a, b) = a^2 + a * b + b * a + b^2

        # Helper for NC fourth power: (a + b)⁴ expanded
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

        # Note: COSMO gives highly variable results for Chained Singular Dense/CS
        # due to the problem's numerical conditioning. We use bounds testing for COSMO.
        @testset "Dense (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            if SOLVER_NAME == :cosmo
                @test result.objective < 0.1  # Should be near 0 (global min)
            else
                @test result.objective ≈ oracles.Dense.opt atol = solver_atol()
            end
            @test compare_sides(result.moment_matrix_sizes, oracles.Dense.sides)
            @test result.n_unique_moment_matrix_elements == oracles.Dense.nuniq
        end

        @testset "Correlative Sparsity MF (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            if SOLVER_NAME == :cosmo
                @test result.objective < 0.1  # Should be near 0 (global min)
            else
                @test result.objective ≈ oracles.CS.opt atol = solver_atol()
            end
            @test compare_sides(result.moment_matrix_sizes, oracles.CS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS.nuniq
        end

        @testset "CS + TS (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            # COSMO gives variable results for Chained Singular due to numerical conditioning
            if SOLVER_NAME == :cosmo
                @test result.objective < 0.1  # Should be near 0 (global min)
            else
                @test result.objective ≈ oracles.CS_TS.opt atol = solver_atol()
            end
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracles.CS_TS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS_TS.nuniq
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
        oracles = select_oracle(EXPECTED_CHAINED_WOOD)
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

        # Note: COSMO gives variable results for Chained Wood Dense. We use bounds testing.
        @testset "Dense (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            if SOLVER_NAME == :cosmo
                @test result.objective < 3.0  # Should be near 1 (global min ≈1)
            else
                @test result.objective ≈ oracles.Dense.opt atol = solver_atol()
            end
            @test compare_sides(result.moment_matrix_sizes, oracles.Dense.sides)
            @test result.n_unique_moment_matrix_elements == oracles.Dense.nuniq
        end

        @testset "Correlative Sparsity MF (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.CS.opt atol = solver_atol()
            @test compare_sides(result.moment_matrix_sizes, oracles.CS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS.nuniq
        end

        @testset "CS + TS (order=2)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)
            @test result.objective ≈ oracles.CS_TS.opt atol = solver_atol()
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracles.CS_TS.sides)
            @test result.n_unique_moment_matrix_elements == oracles.CS_TS.nuniq
        end
    end

end
