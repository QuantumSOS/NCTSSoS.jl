module PublicPolynomialChecks
using Test, LinearAlgebra, Random
using NCTSSoS

include("matrix_oracles.jl")
include("simplify.jl")
include("monomials.jl")
include("basis.jl")
include("composed_monomial.jl")
include("arithmetic.jl")

function check()
    @testset "Public polynomial oracles" begin
        matrix_oracles()
        site_matrix_oracles()
        simplify_oracles()
        monomial_powers()
        basis_contracts()
        mixed_contracts()
        arithmetic_contracts()
    end
end
end

PublicPolynomialChecks.check()
