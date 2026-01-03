# =============================================================================
# Bilocal Network Tests
# =============================================================================
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
# NOTE: This test requires --local (Mosek license) and order=3.
# =============================================================================

using NCTSSoS, Test
using JuMP

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "..", "standalone_setup.jl"))

@testset "Bilocal Networks" begin
    @testset "Example 8.1.1 - Bilocal Quantum Bound" begin
        # Create unipotent variables for three parties
        # Index mapping:
        #   1,2 → x[1], x[2] (Alice)
        #   3,4 → y[1], y[2] (Bob)
        #   5,6 → z[1], z[2] (Charlie)
        reg, (x, y, z) = create_unipotent_variables([("x", 1:2), ("y", 1:2), ("z", 1:2)])

        # Helper to create state expectation from indices
        # Maps: 1,2 → x[1],x[2]; 3,4 → y[1],y[2]; 5,6 → z[1],z[2]
        vars = [x[1], x[2], y[1], y[2], z[1], z[2]]
        function make_expectation(indices)
            m = one(typeof(x[1]))
            for idx in indices
                m = m * vars[idx]
            end
            return ς(m)
        end

        # Support structure from Example 8.1.1
        # Each element is either:
        #   - [[a;b;c]] → single expectation ς(xyz)
        #   - [[a;b;c], [d;e;f]] → product of expectations ς(xyz) * ς(xyz')
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
        # Quadratic terms: 1/8 * [1, 1, ..., -2] (36 terms)
        coe_quad = (1 / 8) .* [
            1, 1, 1, 1, 1, 1, 1, 1,           # Squared terms (8)
            2, 2, 2, -2, 2, 2, -2,             # Cross terms row 1 (7)
            2, 2, -2, 2, 2, -2,                # Cross terms row 2 (6)
            2, -2, 2, 2, -2,                   # Cross terms row 3 (5)
            -2, 2, 2, -2,                      # Cross terms row 4 (4)
            -2, -2, 2, 2, -2, -2               # Cross terms row 5 (6)
        ]
        # Linear terms: -[1, 1, 1, 1, 1, -1, -1, 1] (8 terms)
        coe_linear = Float64[-1, -1, -1, -1, -1, 1, 1, -1]

        # Build the state polynomial by accumulation
        # Start with the first quadratic term
        sp = coe_quad[1] * make_expectation(supp[1][1]) * make_expectation(supp[1][2])

        # Add remaining quadratic terms (products of state expectations)
        for i in 2:36
            term = supp[i]
            sp = sp + coe_quad[i] * make_expectation(term[1]) * make_expectation(term[2])
        end

        # Add linear terms (single state expectations)
        for i in 1:8
            term = supp[36+i]
            sp = sp + coe_linear[i] * make_expectation(term[1])
        end

        # Convert to NCStatePolynomial for optimization
        spop = polyopt(sp * one(typeof(x[1])), reg)

        # Solve with order=3 relaxation
        solver_config = SolverConfig(optimizer=SOLVER, order=3)
        result = cs_nctssos(spop, solver_config)

        # Expected quantum bound is -4
        @test result.objective ≈ -4.0 atol = 1e-4
    end
end
