# =============================================================================
# test/polynomials/runtests.jl
# =============================================================================
# Tests: Core polynomial algebra (types, arithmetic, simplification)
# Dependencies: none (no solver needed)
# Requires --local: no
#
# Coverage:
# ---------
# - algebra_types.jl: Algebra type system (Pauli, Fermionic, Bosonic, etc.)
# - variables.jl: Variable creation and registry
# - monomials.jl: Monomial operations
# - term.jl: Term (coefficient * monomial) operations
# - composed_monomial.jl: Trace/state monomials
# - polynomial.jl: Polynomial construction and access
# - arithmetic.jl: Addition, multiplication, scalar operations
# - compare.jl: Equality and ordering
# - canonicalization.jl: Canonical form computation
# - simplify.jl: Algebra-specific simplification rules
# - matrix_oracles.jl: Matrix representation utilities
# - basis.jl: Basis generation
# - state_word.jl, statepolynomial.jl, state_basis.jl: State polynomials
# - utils.jl, allocations.jl: Utilities and performance tests
# =============================================================================

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
