# 1D Transverse Field Ising Model Tests
# Tests polynomial optimization on the transverse field Ising model.

using Test, NCTSSoS, JuMP

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

# USE_LOCAL fallback for standalone
if !@isdefined(USE_LOCAL)
    const USE_LOCAL = true  # Standalone assumes Mosek
end

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
