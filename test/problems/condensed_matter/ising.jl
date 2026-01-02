# =============================================================================
# 1D Transverse Field Ising Model Tests
# =============================================================================
# Tests polynomial optimization on the transverse field Ising model.
# Requires --local flag (Mosek) for sufficient precision.
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "..", "setup.jl"))

if USE_LOCAL
    @testset "1D Transverse Field Ising Model" begin
        N = 3
        registry, (sx, sy, sz) = create_pauli_variables(1:N)

        J = 1.0
        h = 2.0
        for (periodic, true_ans) in zip((true, false), (-1.0175918, -1.0104160))
            ham = sum(-complex(J / 4) * sz[i] * sz[mod1(i + 1, N)] for i in 1:(periodic ? N : N - 1)) + sum(-h / 2 * sx[i] for i in 1:N)

            pop = polyopt(ham, registry)

            solver_config = SolverConfig(optimizer=SOLVER, order=2)

            res = cs_nctssos(pop, solver_config)
            @test res.objective / N ≈ true_ans atol = 1e-6
        end
    end
end
