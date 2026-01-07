# =============================================================================
# test/problems/fermionic/runtests.jl
# =============================================================================
# Tests: Fermionic systems (parity superselection, fermionic chains)
# Dependencies: SOLVER
# Requires --local: yes (parent runtests.jl guards with `if USE_LOCAL`)
# =============================================================================

using Test

@testset "Fermionic Systems" begin
    include("fermionic.jl")
    include("fermionic_chain.jl")
end
