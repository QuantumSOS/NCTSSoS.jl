# test/problems/nc_polynomial/nc_correlative.jl
# Tests: NC Correlative Sparsity - Constrained example (3 variables)
#
# Tests correlative sparsity exploitation with constraints.

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
const NC_CORRELATIVE_ORACLES = (
    CS_d3 = (opt=0.9975308091952613, sides=[15, 15], nuniq=149),
    TS_d3 = (opt=0.9975305666745705, nuniq=123),  # sides vary
)

if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

@testset "NC Correlative Sparsity" begin
    n = 3
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])

    f = 1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
        2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
        x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
        9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
        54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0 / 3 for i = 1:n])
    pop = polyopt(f, reg; ineq_constraints=cons)

    @testset "CS MF (order=3)" begin
        oracle = NC_CORRELATIVE_ORACLES.CS_d3
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "CS MF (SOS)" begin
        oracle = NC_CORRELATIVE_ORACLES.CS_d3
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=MF()
        )
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ oracle.opt atol = 1e-5
    end

    @testset "TS MMD (order=3)" begin
        oracle = NC_CORRELATIVE_ORACLES.TS_d3
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config; dualize=false)
        result = cs_nctssos_higher(pop, result, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "TS MMD (SOS)" begin
        oracle = NC_CORRELATIVE_ORACLES.TS_d3
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ oracle.opt atol = 1e-4
    end
end
