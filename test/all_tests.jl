# =============================================================================
# test/all_tests.jl - Full Test Suite
# =============================================================================
# Run: julia --project -e 'include("test/all_tests.jl")'
#
# Runs everything: polynomials, quality, relaxations, small tests, large tests.
# Requires Mosek for large_tests.jl.
# =============================================================================

using Test, NCTSSoS

@testset "NCTSSoS.jl Full Suite" begin
    # Polynomial algebra (no solver needed)
    include("polynomials.jl")

    # Quality checks
    include("quality.jl")

    # Relaxation components (includes TestUtils.jl if needed)
    include("relaxations.jl")

    # Fast problem tests (uses SOLVER from relaxations)
    include("small_tests.jl")

    # Slow/Mosek problem tests
    include("large_tests.jl")
end
