# =============================================================================
# 1D Heisenberg Chain Tests
# =============================================================================
# Tests polynomial optimization on the 1D Heisenberg chain model.
# Requires --local flag (Mosek) for sufficient precision.
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "..", "standalone_setup.jl"))

if USE_LOCAL
    @testset "1D Heisenberg Chain" begin
        N = 6
        registry, (sx, sy, sz) = create_pauli_variables(1:N)

        ham = sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [sx, sy, sz] for i in 1:N)

        pop = polyopt(ham, registry)

        solver_config = SolverConfig(optimizer=SOLVER, order=2)

        res = cs_nctssos(pop, solver_config)

        @test res.objective / N ≈ -0.467129 atol = 1e-6
    end
end
