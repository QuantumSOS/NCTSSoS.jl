using Test, NCTSSoS

@testset "Polynomials" begin
    include("algebra_types.jl")
    include("variables.jl")
    include("monomials.jl")
    include("term.jl")
    include("composed_monomial.jl")
    include("polynomial.jl")
    include("arithmetic.jl")
    include("compare.jl")
    include("canonicalization.jl")
    include("simplify.jl")
    include("matrix_oracles.jl")
    include("basis.jl")
    include("state_word.jl")
    include("statepolynomial.jl")
    include("state_basis.jl")
    include("utils.jl")
    include("allocations.jl")
end
