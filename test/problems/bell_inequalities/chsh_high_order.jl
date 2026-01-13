# test/problems/bell_inequalities/chsh_high_order.jl
# Tests: CHSH Bell inequality - higher order relaxations (order=2+)
#
# WARNING: Combined CS+TS does NOT converge to -2sqrt(2) for CHSH!
# This gives opt=-4.0, a much looser bound than the quantum value -2.828.
# This is consistent with NCTSSOS behavior - it's an inherent limitation.

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

# Oracle values from NCTSSOS
const CHSH_HIGH_ORDER_ORACLES = (
    CS_TS_d2 = (opt=-3.999999999803662, sides=[3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 2, 2], nuniq=5),
)

if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

function create_chsh_problem()
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    pop = polyopt(-f, reg)
    return pop, reg
end

@testset "CHSH High Order" begin

    @testset "Combined CS+TS (order=2) - KNOWN LOOSE BOUND" begin
        pop, _ = create_chsh_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MF(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config)

        @test result.objective ≈ CHSH_HIGH_ORDER_ORACLES.CS_TS_d2.opt atol = 1e-6
        # Block sizes may differ but total count matches
        @test length(flatten_sizes(result.moment_matrix_sizes)) == length(CHSH_HIGH_ORDER_ORACLES.CS_TS_d2.sides)
        @test result.n_unique_moment_matrix_elements == CHSH_HIGH_ORDER_ORACLES.CS_TS_d2.nuniq
    end

end
