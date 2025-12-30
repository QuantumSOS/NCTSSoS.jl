# =============================================================================
# Physics Model Tests Runner
# =============================================================================
# Tests for real-world physics applications:
#   - Bell inequalities (CHSH, I_3322)
#   - Bilocal networks (quantum network nonlocality)
#   - Heisenberg model
#   - XY model
#   - Bose-Hubbard model
#   - Fermionic systems
#   - PXP model (Rydberg atoms)
#
# NOTE: These tests require --local flag (Mosek license)
#       as they involve larger SDP problems.
# =============================================================================

using Test

@testset "Physics Models" begin
    include("heisenberg.jl")
    include("xy_model.jl")
    include("bose_hubbard.jl")
    include("bell_inequalities.jl")
    include("bilocal_networks.jl")
    include("fermionic.jl")
    include("fermionic_chain.jl")
    include("pxp.jl")
end
