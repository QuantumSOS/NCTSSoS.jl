# FastPolynomials test suite
# This loads FastPolynomials directly (NCTSSoS migration is Phase 3)

include("setup.jl")

@testset "FastPolynomials" begin
    include("variables.jl")
    include("monomials.jl")
    include("polynomial.jl")
    include("arithmetic.jl")
    include("compare.jl")
    include("simplify.jl")
    include("state_word.jl")
    include("statepolynomial.jl")
    include("utils.jl")
    include("allocations.jl")
end
