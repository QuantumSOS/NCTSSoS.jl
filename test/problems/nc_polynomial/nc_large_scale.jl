# =============================================================================
# test/problems/nc_polynomial/nc_large_scale.jl
# =============================================================================
# Tests: NC Large-scale sparsity example (n=10 variables)
# Dependencies: SOLVER
# Requires --local: yes (large SDP, needs Mosek for performance)
#
# Validates CS+TS on large-scale problems.
# =============================================================================

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

# USE_LOCAL flag for large-scale tests
@isdefined(USE_LOCAL) || (USE_LOCAL = true)

# Oracle values from NCTSSOS
const NC_LARGE_SCALE_ORACLES = (
    CS_TS_d3 = (opt=3.011288353315061, nblocks=1544, nuniq=1165),
)

flatten_sizes(sizes) = reduce(vcat, sizes)

if USE_LOCAL
    @testset "NC Large Scale (n=10)" begin
        n = 10
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])

        # Build polynomial
        f = Polynomial{NonCommutativeAlgebra,UInt8,Float64}(
            Term{NormalMonomial{NonCommutativeAlgebra,UInt8},Float64}[]
        )
        for i = 1:n
            jset = max(1, i - 5):min(n, i + 1)
            jset = setdiff(jset, i)
            f += (2.0 * x[i] + 5.0 * x[i]^3 + 1)^2
            f -= sum([
                4.0 * x[i] * x[j] +
                10.0 * x[i]^3 * x[j] +
                2.0 * x[j] +
                4.0 * x[i] * x[j]^2 +
                10.0 * x[i]^3 * x[j]^2 +
                2.0 * x[j]^2 for j in jset
            ])
            f += sum([
                1.0 * x[j] * x[k] + 2.0 * x[j]^2 * x[k] + 1.0 * x[j]^2 * x[k]^2
                for j in jset for k in jset
            ])
        end

        cons = vcat([(1.0 - x[i]^2) for i = 1:n], [(1.0 * x[i] - 1.0 / 3) for i = 1:n])
        pop = polyopt(f, reg; ineq_constraints=cons)

        oracle = NC_LARGE_SCALE_ORACLES.CS_TS_d3

        @testset "Moment Method (order=3)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config; dualize=false)
            @test result.objective ≈ oracle.opt atol = 1e-4
            @test length(flatten_sizes(result.moment_matrix_sizes)) == oracle.nblocks
            @test result.n_unique_moment_matrix_elements == oracle.nuniq
        end

        @testset "SOS Dualization (order=3)" begin
            config = SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config; dualize=true)
            @test result.objective ≈ oracle.opt atol = 1e-6
        end
    end
end
