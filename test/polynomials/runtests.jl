using Test
# NCTSSoS is already loaded by parent runtests.jl via `using NCTSSoS`
# All polynomial types are exported directly from NCTSSoS

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
    include("basis.jl")
    include("state_word.jl")
    include("statepolynomial.jl")
    include("state_basis.jl")
    include("utils.jl")
    include("allocations.jl")
end
