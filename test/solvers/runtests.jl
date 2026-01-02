# =============================================================================
# Solver Tests Runner
# =============================================================================
# Runs all solver-related tests organized by:
#
# 1. Component Tests (utilities, algorithms):
#    - interface.jl: PolyOpt constructor, basic solver tests
#    - sos.jl: SOS dualization components (Cαj)
#    - sparsity.jl: Correlative/term sparsity components
#    - gns.jl: GNS reconstruction
#
# 2. Problem-Based Tests (in problems/ subdirectory):
#    - chsh.jl: CHSH Bell inequality (moment, SOS, state, trace)
#    - i3322.jl: I_3322 Bell inequality
#    - heisenberg_star.jl: Heisenberg model on star graph
#    - nc_examples.jl: NC polynomial examples (Example 1, 2, correlative sparsity)
#    - state_polynomial.jl: State polynomial examples (7.2.x)
#    - trace_polynomial.jl: Trace polynomial examples (6.x)
#
# 3. Benchmark Tests:
#    - ncpop_benchmarks.jl: Classical optimization benchmarks
#
# Prerequisites:
#   - SOLVER must be defined (from test/setup.jl or parent runtests.jl)
#   - USE_LOCAL determines whether Mosek or COSMO is used
# =============================================================================

using Test

@testset "Solvers" begin
    # =========================================================================
    # Component Tests
    # =========================================================================
    @testset "Components" begin
        include("interface.jl")
        include("sos.jl")
        include("sparsity.jl")
        include("gns.jl")
    end

    # =========================================================================
    # Problem-Based Tests
    # =========================================================================
    @testset "Problems" begin
        include("problems/chsh.jl")
        include("problems/i3322.jl")
        include("problems/heisenberg_star.jl")
        include("problems/nc_examples.jl")
        include("problems/state_polynomial.jl")
        include("problems/trace_polynomial.jl")
    end

    # =========================================================================
    # Benchmark Tests
    # =========================================================================
    include("ncpop_benchmarks.jl")
end
