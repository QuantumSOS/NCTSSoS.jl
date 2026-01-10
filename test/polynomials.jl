# =============================================================================
# test/polynomials.jl - Polynomial Algebra Tests
# =============================================================================
# Run: julia --project -e 'include("test/polynomials.jl")'
#
# Tests core polynomial algebra (types, arithmetic, simplification).
# No solver needed.
# =============================================================================

using Test, NCTSSoS

@testset "Polynomials" begin
    include("polynomials/algebra_types.jl")
    include("polynomials/variables.jl")
    include("polynomials/monomials.jl")
    include("polynomials/composed_monomial.jl")
    include("polynomials/polynomial.jl")
    include("polynomials/arithmetic.jl")
    include("polynomials/compare.jl")
    include("polynomials/canonicalization.jl")
    include("polynomials/simplify.jl")
    include("polynomials/matrix_oracles.jl")
    include("polynomials/basis.jl")
    include("polynomials/state_word.jl")
    include("polynomials/statepolynomial.jl")
    include("polynomials/state_basis.jl")
    include("polynomials/utils.jl")
    include("polynomials/allocations.jl")
end
