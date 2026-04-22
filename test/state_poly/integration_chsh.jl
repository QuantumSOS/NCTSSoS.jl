# test/state_poly/integration_chsh.jl
# Tests: CHSH Bell inequality - State Polynomial formulation
#
# Uses state polynomial formulation with ς() state function.

using Test, NCTSSoS, JuMP

# Expectations in test/data/expectations/chsh_state.toml

@testset "CHSH State Polynomial" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp = -ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
    spop = polyopt(sp * one(typeof(x[1])), reg)

    @testset "Dense" begin
        oracle = expectations_oracle("expectations/chsh_state.toml", "Dense")
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(spop, config)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Term Sparsity (MMD)" begin
        oracle = expectations_oracle("expectations/chsh_state.toml", "TS")
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(spop, config)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Symmetry is rejected on the state-polynomial solver path" begin
        symmetry = SymmetrySpec(SignedPermutation(
            x[1].word[1] => x[2].word[1],
            x[2].word[1] => x[1].word[1],
        ))

        err = try
            cs_nctssos(
                spop,
                SolverConfig(
                    optimizer=SOLVER,
                    order=1,
                    cs_algo=NoElimination(),
                    ts_algo=NoElimination(),
                    symmetry=symmetry,
                ),
            )
            nothing
        catch caught
            caught
        end

        @test err isa ArgumentError
        @test occursin(
            "state/trace polynomial optimization",
            sprint(showerror, err),
        )
    end
end
