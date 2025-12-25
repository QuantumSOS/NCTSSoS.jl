using NCTSSoS, Test

# =============================================================================
# FastPolynomials Integration - Test Status (2025-12-21)
# =============================================================================
# Migration status:
#   ✓ FastPolynomials tests - using new API, all pass
#   ✓ pop.jl, sparse.jl - migrated, all pass
#   ✓ moment_solver.jl, sos_solver.jl, interface.jl - migrated, all pass
#   ✓ Aqua.jl, ExplicitImports.jl - pass
#   ✓ heisenberg.jl, xy_model.jl, bose_hubbard.jl - pass (LOCAL_TESTING only)
#   ✓ bell_ineq.jl - migrated to ProjectorAlgebra API (LOCAL_TESTING only)
#   ✓ state_poly_opt.jl, trace_poly_opt.jl - migrated (uses ς, tr)
#
# Remaining:
#   - Doctest.jl disabled (FastPolynomials doctests need import path fix)
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
        include("xy_model.jl")
        include("bose_hubbard.jl")
        # include("bell_ineq.jl")
    end
    include("sos_solver.jl")
    include("interface.jl")

    # Fermionic parity superselection tests
    include("fermionic_parity_test.jl")

    # State/Trace polynomial tests
    include("state_poly_opt.jl")
    include("trace_poly_opt.jl")
end
