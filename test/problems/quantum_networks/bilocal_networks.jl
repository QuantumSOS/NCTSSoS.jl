# Bilocal Network Tests
# Tests for quantum bilocal network optimization from Example 8.1.1 in the
# state polynomial paper.
#
# Bilocal network scenario:
#        Source 1          Source 2
#           |                  |
#           v                  v
#     Alice --- Bob --- Charlie
#         (x)    (y)      (z)
#
# Variables:
#   - x[1], x[2] : Alice's measurement operators (dichotomic, x² = 1)
#   - y[1], y[2] : Bob's measurement operators
#   - z[1], z[2] : Charlie's measurement operators
#
# The bilocal constraint arises from the independence of the two sources.
# The objective is a state polynomial that tests whether correlations can
# be explained by the bilocal model.
#
# Expected result: Quantum bound of -4
#

using NCTSSoS, Test
using JuMP

# Solver: use Mosek if available, otherwise error
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

# Expected values from NCTSSOS
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization → -4.0 quantum bound)
#   sides = moment matrix block sizes
#   nuniq = unique moment indices (affine constraints)
const EXPECTED_BILOCAL = (
    # Dense: single block of size 637, 15867 constraints
    Dense_d3 = (opt=-3.9999999896685883, sides=[637], nuniq=15867),
    # TS: blocks [16×2, 5×12, 3×24, 1×572], 263 constraints
    TS_d3 = (opt=-3.9999999994438804, sides=[fill(16, 2); fill(5, 12); fill(3, 24); fill(1, 572)], nuniq=263),
)

# Helper: flatten moment_matrix_sizes for comparison
if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

@testset "Bilocal Networks" begin

    # =========================================================================
    # Problem Setup Helper
    # =========================================================================
    function create_bilocal_problem()
        # Create unipotent variables for three parties
        # Index mapping: 1,2 → x; 3,4 → y; 5,6 → z
        reg, (x, y, z) = create_unipotent_variables([("x", 1:2), ("y", 1:2), ("z", 1:2)])

        # Helper to create state expectation from indices
        vars = [x[1], x[2], y[1], y[2], z[1], z[2]]
        function make_expectation(indices)
            m = one(typeof(x[1]))
            for idx in indices
                m = m * vars[idx]
            end
            return ς(m)
        end

        # Support structure from Example 8.1.1
        supp = [
            # Squared terms ς(·)ς(·) with same argument (8 terms)
            [[1, 3, 5], [1, 3, 5]], [[1, 3, 6], [1, 3, 6]], [[2, 3, 5], [2, 3, 5]],
            [[2, 3, 6], [2, 3, 6]], [[1, 4, 5], [1, 4, 5]], [[1, 4, 6], [1, 4, 6]],
            [[2, 4, 5], [2, 4, 5]], [[2, 4, 6], [2, 4, 6]],
            # Cross terms ς(·)ς(·') with different arguments (28 terms)
            [[1, 3, 5], [1, 3, 6]], [[1, 3, 5], [2, 3, 5]], [[1, 3, 5], [2, 3, 6]],
            [[1, 3, 5], [1, 4, 5]], [[1, 3, 5], [1, 4, 6]], [[1, 3, 5], [2, 4, 5]],
            [[1, 3, 5], [2, 4, 6]], [[1, 3, 6], [2, 3, 5]], [[1, 3, 6], [2, 3, 6]],
            [[1, 3, 6], [1, 4, 5]], [[1, 3, 6], [1, 4, 6]], [[1, 3, 6], [2, 4, 5]],
            [[1, 3, 6], [2, 4, 6]], [[2, 3, 5], [2, 3, 6]], [[2, 3, 5], [1, 4, 5]],
            [[2, 3, 5], [1, 4, 6]], [[2, 3, 5], [2, 4, 5]], [[2, 3, 5], [2, 4, 6]],
            [[2, 3, 6], [1, 4, 5]], [[2, 3, 6], [1, 4, 6]], [[2, 3, 6], [2, 4, 5]],
            [[2, 3, 6], [2, 4, 6]], [[1, 4, 5], [1, 4, 6]], [[1, 4, 5], [2, 4, 5]],
            [[1, 4, 5], [2, 4, 6]], [[1, 4, 6], [2, 4, 5]], [[1, 4, 6], [2, 4, 6]],
            [[2, 4, 5], [2, 4, 6]],
            # Linear terms ς(·) (8 terms)
            [[1, 3, 5]], [[1, 3, 6]], [[2, 3, 5]], [[2, 3, 6]],
            [[1, 4, 5]], [[1, 4, 6]], [[2, 4, 5]], [[2, 4, 6]]
        ]

        # Coefficients from the paper
        coe_quad = (1 / 8) .* [
            1, 1, 1, 1, 1, 1, 1, 1,           # Squared terms (8)
            2, 2, 2, -2, 2, 2, -2,             # Cross terms row 1 (7)
            2, 2, -2, 2, 2, -2,                # Cross terms row 2 (6)
            2, -2, 2, 2, -2,                   # Cross terms row 3 (5)
            -2, 2, 2, -2,                      # Cross terms row 4 (4)
            -2, -2, 2, 2, -2, -2               # Cross terms row 5 (6)
        ]
        coe_linear = Float64[-1, -1, -1, -1, -1, 1, 1, -1]

        # Build state polynomial
        sp = coe_quad[1] * make_expectation(supp[1][1]) * make_expectation(supp[1][2])
        for i in 2:36
            term = supp[i]
            sp = sp + coe_quad[i] * make_expectation(term[1]) * make_expectation(term[2])
        end
        for i in 1:8
            term = supp[36+i]
            sp = sp + coe_linear[i] * make_expectation(term[1])
        end

        spop = polyopt(sp * one(typeof(x[1])), reg)
        return spop, reg
    end

    # =========================================================================
    # Dense (No Sparsity)
    # =========================================================================
    @testset "Dense (order=3)" begin
        spop, _ = create_bilocal_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(spop, config)

        @test result.objective ≈ EXPECTED_BILOCAL.Dense_d3.opt atol = 1e-6
        @test_broken flatten_sizes(result.moment_matrix_sizes) == EXPECTED_BILOCAL.Dense_d3.sides
        @test result.n_unique_moment_matrix_elements == EXPECTED_BILOCAL.Dense_d3.nuniq
    end

    # =========================================================================
    # Term Sparsity (MMD)
    # =========================================================================
    @testset "Term Sparsity MMD (order=3)" begin
        spop, _ = create_bilocal_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(spop, config)

        @test result.objective ≈ EXPECTED_BILOCAL.TS_d3.opt atol = 1e-6
        @test_broken flatten_sizes(result.moment_matrix_sizes) == EXPECTED_BILOCAL.TS_d3.sides
        @test result.n_unique_moment_matrix_elements == EXPECTED_BILOCAL.TS_d3.nuniq
    end
end
