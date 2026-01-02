# =============================================================================
# Relaxation Component Tests Runner
# =============================================================================
# Tests for relaxation algorithm components:
#   - interface.jl: PolyOpt/StatePolyOpt constructors, dualization
#   - sos.jl: SOS dualization components (Cαj)
#   - sparsity.jl: Correlative/term sparsity components
#   - gns.jl: GNS reconstruction
#
# Prerequisites:
#   - SOLVER must be defined (from test/setup.jl or parent runtests.jl)
# =============================================================================

using Test

@testset "Relaxations" begin
    include("interface.jl")
    include("sos.jl")
    include("sparsity.jl")
    include("gns.jl")
end
