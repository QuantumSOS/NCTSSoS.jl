# test/problems/fermionic/runtests.jl
# Tests: Fermionic systems (parity superselection, fermionic chains)

using Test

@testset "Fermionic Systems" begin
    include("fermionic.jl")
    include("fermionic_chain.jl")
end
