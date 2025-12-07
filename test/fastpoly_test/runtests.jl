using Test
# Use FastPolynomials from NCTSSoS (loaded by runtests.jl via `using NCTSSoS`)
using NCTSSoS.FastPolynomials

# Import @ncpolyvar explicitly to avoid ambiguity
import NCTSSoS.FastPolynomials: @ncpolyvar

@testset "FastPolynomials" begin
    include("variables.jl")
    include("monomials.jl")
    include("polynomial.jl")
    include("arithmetic.jl")
    include("compare.jl")
    include("simplify.jl")
    include("basis.jl")
    include("state_word.jl")
    include("statepolynomial.jl")
    include("utils.jl")
    include("allocations.jl")
end
