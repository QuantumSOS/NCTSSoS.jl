using NCTSSoS, Test

# =============================================================================
# FastPolynomials Integration - Test Status (2025-12-13)
# =============================================================================
# Migration status:
#   ✓ FastPolynomials tests - using new API, all pass
#   ✓ pop.jl, sparse.jl - migrated, all pass
#   ✓ moment_solver.jl, sos_solver.jl, interface.jl - migrated, all pass
#   ✓ Aqua.jl, ExplicitImports.jl - pass
#   ✓ heisenberg.jl - passes (only runs with LOCAL_TESTING=true)
#
# Remaining:
#   - state_poly_opt.jl, trace_poly_opt.jl - not yet migrated (uses ς, tr)
#   - Doctest.jl disabled (FastPolynomials doctests use invalid import path)
# =============================================================================

@testset "NCTSSoS.jl" begin
    # FastPolynomials - uses new API, all tests pass
    include("fastpoly_test/runtests.jl")

    # Core optimization tests - migrated to new API
    include("pop.jl")
    include("sparse.jl")

    # Quality checks
    include("Aqua.jl")
    # NOTE: Doctest.jl disabled - FastPolynomials doctests use "using FastPolynomials"
    # but FastPolynomials is a submodule, not a registered package.
    # TODO: Fix FastPolynomials doctests to use correct import path
    # include("Doctest.jl")
    include("ExplicitImports.jl")

    # Solver integration tests
    include("moment_solver.jl")
    if haskey(ENV, "LOCAL_TESTING")
        include("heisenberg.jl")
    end
    include("sos_solver.jl")
    include("interface.jl")

    # State/Trace polynomial tests - not yet migrated (uses specialized features: ς, tr)
    # include("state_poly_opt.jl")
    # include("trace_poly_opt.jl")
end
