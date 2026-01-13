# 1D Heisenberg Chain Tests
# Tests polynomial optimization on the 1D Heisenberg chain model.

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
