# =============================================================================
# test/problems/condensed_matter/runtests.jl
# =============================================================================
# Tests: Condensed matter physics models (Heisenberg, Ising, XY, Bose-Hubbard, PXP)
# Dependencies: SOLVER
# Requires --local: yes (parent runtests.jl guards with `if USE_LOCAL`)
# =============================================================================

using Test

@testset "Condensed Matter" begin
    include("ising.jl")
    include("heisenberg_chain.jl")
    include("heisenberg_star.jl")
    include("heisenberg.jl")
    include("xy_model.jl")
    include("bose_hubbard.jl")
    include("pxp.jl")
end
