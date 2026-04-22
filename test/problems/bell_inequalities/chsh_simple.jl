# test/problems/bell_inequalities/chsh_simple.jl
# Tests: CHSH Bell inequality - basic sparsity configurations (order=1)
#
# Coverage: Dense, Correlative Sparsity (MF), Term Sparsity (MMD),
# symmetry-reduced dense relaxation
# Expected optimal value: -2sqrt(2) ≈ -2.8284 (quantum bound)

using Test, NCTSSoS, JuMP

# Expectations in test/data/expectations/chsh_simple.toml

function create_chsh_problem()
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    pop = polyopt(-f, reg)
    return pop, reg, x, y
end

function create_chsh_symmetry(x, y)
    alice_swap = SignedPermutation(
        x[1].word[1] => x[2].word[1],
        x[2].word[1] => x[1].word[1],
        y[2].word[1] => (-1, y[2].word[1]),
    )
    bob_swap = SignedPermutation(
        x[2].word[1] => (-1, x[2].word[1]),
        y[1].word[1] => y[2].word[1],
        y[2].word[1] => y[1].word[1],
    )
    party_swap = SignedPermutation(
        x[1].word[1] => y[1].word[1],
        x[2].word[1] => y[2].word[1],
        y[1].word[1] => x[1].word[1],
        y[2].word[1] => x[2].word[1],
    )
    return SymmetrySpec(alice_swap, bob_swap, party_swap)
end

@testset "CHSH Simple (order=1)" begin

    @testset "Dense" begin
        oracle = expectations_oracle("expectations/chsh_simple.toml", "Dense_d1")
        pop, _, _, _ = create_chsh_problem()
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

    @testset "Correlative Sparsity (MF)" begin
        oracle = expectations_oracle("expectations/chsh_simple.toml", "CS_d1")
        pop, _, _, _ = create_chsh_problem()
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

    @testset "Term Sparsity (MMD)" begin
        oracle = expectations_oracle("expectations/chsh_simple.toml", "TS_d1")
        pop, _, _, _ = create_chsh_problem()
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

    @testset "Dense + Symmetry" begin
        oracle = TestExpectations.expectations_case(
            TestExpectations.expectations_load("expectations/chsh_simple.toml"),
            "Symmetry_d1",
        )["expected"]
        pop, _, x, y = create_chsh_problem()
        Q = [one(x[1]), x[1], x[2], y[1], y[2]]
        config = SolverConfig(
            optimizer=SOLVER,
            moment_basis=Q,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=create_chsh_symmetry(x, y),
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ oracle["objective"] atol = 1e-6
        @test !isnothing(result.symmetry)
        @test result.symmetry.group_order == oracle["group_order"]
        @test result.symmetry.invariant_moment_count == oracle["invariant_moment_count"]
        @test result.symmetry.psd_block_sizes == oracle["psd_block_sizes"]

        shown = sprint(show, result)
        @test occursin("Moment Matrix Block Sizes (pre-symmetry): [5]", shown)
        @test occursin("Symmetry-Reduced PSD Block Sizes: [1, 1, 1]", shown)
    end

end
