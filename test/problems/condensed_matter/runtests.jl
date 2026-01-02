# =============================================================================
# Condensed Matter Physics Tests Runner
# =============================================================================
# Tests for condensed matter physics models:
#   - Heisenberg models (chain, star graph)
#   - Ising model (transverse field)
#   - XY model
#   - Bose-Hubbard model
#   - PXP model (Rydberg atoms)
#
# NOTE: These tests are commented out as they take too long or require
#       --local flag (Mosek license) and are not yet validated.
# =============================================================================

using Test

#=
@testset "Condensed Matter" begin
    include("ising.jl")
    include("heisenberg_chain.jl")
    include("heisenberg_star.jl")
    include("heisenberg.jl")
    include("xy_model.jl")
    include("bose_hubbard.jl")
    include("pxp.jl")
end
=#
