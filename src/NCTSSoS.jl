module NCTSSoS

using SparseArrays, LinearAlgebra, JuMP
using CliqueTrees, ChordalGraph, Graphs
using CliqueTrees: EliminationAlgorithm, SupernodeType
import CliqueTrees.cliquetree

include("FastPolynomials/src/FastPolynomials.jl")
using .FastPolynomials

# Import types (new and legacy compatibility)
using .FastPolynomials: Monomial, Polynomial, Term
using .FastPolynomials: AlgebraType, NonCommutativeAlgebra, PauliAlgebra
using .FastPolynomials: UnipotentAlgebra, ProjectorAlgebra, FermionicAlgebra, BosonicAlgebra
using .FastPolynomials: StateWord, NCStateWord, StatePolynomial, NCStatePolynomial
using .FastPolynomials: Arbitrary, MaxEntangled
using .FastPolynomials: VariableRegistry, symbols, indices, index_type, algebra_type, subregistry
using .FastPolynomials: variable_indices

# Import legacy compatibility types
using .FastPolynomials: AbstractPolynomial, Variable

# Import functions
using .FastPolynomials: sorted_union, sorted_unique
using .FastPolynomials: monomials, coefficients, terms, degree, maxdegree, variables
using .FastPolynomials: get_basis, get_ncbasis, get_ncbasis_deg
using .FastPolynomials: simplify, simplify!, canonicalize
using .FastPolynomials: neat_dot, _neat_dot3, expval
using .FastPolynomials: star, star!
using .FastPolynomials: monomial

# Export the macro (needs special handling)
using .FastPolynomials: @ncpolyvar

# Re-export key types and functions for downstream users
export @ncpolyvar, ς

# Note: Core types (Monomial, Polynomial, Term, Variable)
# are NOT exported from NCTSSoS to avoid ambiguity when FastPolynomials
# is also loaded. Access them via NCTSSoS.FastPolynomials.Monomial etc.
# or using NCTSSoS.FastPolynomials: Monomial, Polynomial, ...

export polyopt, cpolyopt
export SolverConfig
export NoElimination, MF, MMD, AsIsElimination, MaximalElimination
export cs_nctssos, cs_nctssos_higher
export reconstruct
export pauli_algebra

include("pop.jl")

include("elimination.jl")

include("solver_utils.jl")

include("sparse.jl")

include("moment_solver.jl")

include("complex_moment_solver.jl")

include("sos_solver.jl")

include("gns.jl")

include("interface.jl")

include("algebra_constructors.jl")
end
