# =============================================================================
# Quantum Networks Tests Runner
# =============================================================================
# Tests for quantum network nonlocality problems:
#   - Bilocal networks
# =============================================================================

using Test

@testset "Quantum Networks" begin
    include("bilocal_networks.jl")
end
