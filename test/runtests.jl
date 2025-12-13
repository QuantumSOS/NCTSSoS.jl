using NCTSSoS, Test

# =============================================================================
# FastPolynomials Integration - Test Re-enablement Plan
# =============================================================================
# Tests are commented out during the fastpoly-integration refactor.
# Re-enable incrementally as each phase completes:
#
#   Phase 1.1 complete → pop.jl
#   Phase 1.2 complete → sparse.jl
#   Phase 1.3 complete → moment_solver.jl
#   Phase 2 complete   → algebra_constructors.jl
#   Phase 3 complete   → interface.jl, sos_solver.jl
#   Phase 4 complete   → All remaining tests
#
# FastPolynomials tests should pass throughout (uses new API).
# =============================================================================

@testset "NCTSSoS.jl" begin
    # FastPolynomials - uses new API, should pass throughout refactor
    include("fastpoly_test/runtests.jl")

    # === Commented out during fastpoly-integration ===
    # Core optimization (re-enable after Phase 1)
    include("pop.jl")
    include("sparse.jl")
    # include("solver_utils.jl")

    # Moment solver (re-enable after Phase 1, LOCAL_TESTING only)
    # if haskey(ENV, "LOCAL_TESTING")
    #     include("moment_solver.jl")
    #     include("heisenberg.jl")
    # end

    # Algebra constructors (re-enable after Phase 2)
    # include("algebra_constructors.jl")

    # Interface & SOS (re-enable after Phase 3)
    # include("sos_solver.jl")
    # include("interface.jl")

    # Advanced features (re-enable after Phase 3)
    # include("state_poly_opt.jl")
    # include("trace_poly_opt.jl")

    # Quality checks (re-enable after Phase 4)
    # include("Aqua.jl")
    # include("Doctest.jl")
    # include("ExplicitImports.jl")
end
