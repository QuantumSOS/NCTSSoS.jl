# =============================================================================
# test/large_tests.jl - Large/Slow Tests (Mosek Required)
# =============================================================================
# Run: julia --project -e 'include("test/large_tests.jl")'
#
# Slow tests and tests requiring Mosek solver.
# Will fail loudly if Mosek is not available.
# =============================================================================

using Test, NCTSSoS, JuMP

# Require Mosek - fail loudly if unavailable
try
    using MosekTools
catch e
    error("large_tests.jl requires Mosek. Install MosekTools and obtain a license.")
end

if !@isdefined(SOLVER)
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

@testset "Large Tests (Mosek)" begin
    # Extended Bell inequality tests
    @testset "Bell Inequalities (Extended)" begin
        include("problems/bell_inequalities/chsh_high_order.jl")
        include("problems/bell_inequalities/chsh_state.jl")
        include("problems/bell_inequalities/chsh_trace.jl")
        include("problems/bell_inequalities/bell_inequalities.jl")
    end

    # Extended NC polynomial tests
    @testset "NC Polynomial (Extended)" begin
        include("problems/nc_polynomial/nc_example2.jl")
        include("problems/nc_polynomial/nc_correlative.jl")
        include("problems/nc_polynomial/nc_large_scale.jl")
        include("problems/nc_polynomial/nc_readme.jl")
    end

    # Trace polynomial extended
    @testset "Trace Polynomial (Extended)" begin
        include("problems/trace_polynomial/sparsity_variants.jl")
    end

    # Condensed matter physics
    @testset "Condensed Matter" begin
        include("problems/condensed_matter/ising.jl")
        include("problems/condensed_matter/heisenberg_chain.jl")
        include("problems/condensed_matter/heisenberg_star.jl")
        include("problems/condensed_matter/heisenberg.jl")
        include("problems/condensed_matter/xy_model.jl")
        include("problems/condensed_matter/bose_hubbard.jl")
        include("problems/condensed_matter/pxp.jl")
    end

    # Fermionic systems
    @testset "Fermionic" begin
        include("problems/fermionic/fermionic.jl")
        include("problems/fermionic/fermionic_chain.jl")
    end

    # Quantum networks
    @testset "Quantum Networks" begin
        include("problems/quantum_networks/bilocal_networks.jl")
    end
end
