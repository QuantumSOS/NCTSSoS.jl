# =============================================================================
# test/small_tests.jl - Fast Problem Tests
# =============================================================================
# Run: julia --project -e 'include("test/TestUtils.jl"); include("test/small_tests.jl")'
# With Mosek: julia --project -e 'push!(ARGS, "--local"); include("test/TestUtils.jl"); include("test/small_tests.jl")'
#
# Fast problem tests (~5-7 min with COSMO, faster with Mosek).
# Solver: COSMO (default), Mosek if --local flag passed.
#
# Requires: SOLVER must be defined (include TestUtils.jl first)
# =============================================================================

using Test, NCTSSoS

# Setup solver if not already defined (for standalone execution)
if !@isdefined(SOLVER)
    include("TestUtils.jl")
end

@testset "Small Tests" begin
    # Bell inequalities - fast subset
    @testset "Bell Inequalities" begin
        include("problems/bell_inequalities/chsh_simple.jl")
        include("problems/bell_inequalities/i3322.jl")
    end

    # NC polynomial - fast examples
    @testset "NC Polynomial" begin
        include("problems/nc_polynomial/nc_example1.jl")
    end

    # Benchmarks - sparsity validation
    @testset "Benchmarks" begin
        include("problems/benchmarks/ncpop_benchmarks.jl")
    end

    # State polynomial
    @testset "State Polynomial" begin
        include("problems/state_polynomial/state_polynomial.jl")
    end

    # Trace polynomial
    @testset "Trace Polynomial" begin
        include("problems/trace_polynomial/trace_polynomial.jl")
    end
end
