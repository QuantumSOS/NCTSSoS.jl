# =============================================================================
# Problem-Based Tests Runner
# =============================================================================
# Runs all problem-based tests organized by domain:
#
# 1. Bell Inequalities: CHSH, I_3322, general Bell inequalities
# 2. Condensed Matter: Heisenberg, Ising, XY, Bose-Hubbard, PXP
# 3. NC Polynomial: Example 1, 2, correlative sparsity
# 4. State Polynomial: Section 7.2.x examples
# 5. Trace Polynomial: Section 6.x examples
# 6. Quantum Networks: Bilocal networks
# 7. Fermionic: Fermionic systems and chains
# 8. Benchmarks: Classical optimization benchmarks (Rosenbrock, etc.)
#
# Prerequisites:
#   - SOLVER must be defined (from test/TestUtils.jl or parent runtests.jl)
#   - USE_LOCAL determines whether Mosek or COSMO is used
#
# NOTE: Many tests require --local flag (Mosek license) for sufficient
#       precision on larger SDP problems.
# =============================================================================

using Test

@testset "Problems" begin
    # Bell inequalities (CHSH, I_3322)
    @testset "Bell Inequalities" begin
        include("bell_inequalities/runtests.jl")
    end

    # NC polynomial examples
    @testset "NC Polynomial" begin
        include("nc_polynomial/runtests.jl")
    end

    # State polynomial examples
    @testset "State Polynomial" begin
        include("state_polynomial/runtests.jl")
    end

    # Trace polynomial examples
    @testset "Trace Polynomial" begin
        include("trace_polynomial/runtests.jl")
    end

    # Benchmarks (Rosenbrock, etc.)
    @testset "Benchmarks" begin
        include("benchmarks/runtests.jl")
    end

    # Tests requiring --local (Mosek)
    if USE_LOCAL
        # Condensed matter physics
        @testset "Condensed Matter" begin
            include("condensed_matter/runtests.jl")
        end

        # Quantum networks
        @testset "Quantum Networks" begin
            include("quantum_networks/runtests.jl")
        end

        # Fermionic systems
        @testset "Fermionic" begin
            include("fermionic/runtests.jl")
        end
    end
end
