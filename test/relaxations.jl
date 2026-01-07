# =============================================================================
# test/relaxations.jl - Relaxation Algorithm Tests
# =============================================================================
# Run: julia --project -e 'include("test/TestUtils.jl"); include("test/relaxations.jl")'
#
# Tests relaxation algorithm components.
# Requires: SOLVER (include TestUtils.jl first)
# =============================================================================

using Test, NCTSSoS

# Setup solver if not already defined (for standalone execution)
if !@isdefined(SOLVER)
    include("TestUtils.jl")
end

@testset "Relaxations" begin
    include("relaxations/interface.jl")
    include("relaxations/sos.jl")
    include("relaxations/sparsity.jl")
    include("relaxations/gns.jl")
    include("relaxations/dualization.jl")
end
