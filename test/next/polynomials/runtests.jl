module PublicPolynomialChecks
using Test, LinearAlgebra
using NCTSSoS.Next

include("matrix_oracles.jl")
include("simplify.jl")
include("monomials.jl")

function check()
    @testset "Public polynomial oracles" begin
        matrix_oracles()
        simplify_oracles()
        monomial_powers()
    end
end
end

PublicPolynomialChecks.check()
