module NCTSSoS

using SparseArrays, LinearAlgebra, JuMP
using CliqueTrees, ChordalGraph, Graphs
using CliqueTrees: EliminationAlgorithm, SupernodeType
import CliqueTrees.cliquetree

include("FastPolynomials/src/FastPolynomials.jl")
using .FastPolynomials
using .FastPolynomials: AbstractPolynomial, Variable, Monomial

using .FastPolynomials: sorted_union, monomials, sorted_unique, maxdegree, get_basis, neat_dot, _neat_dot3, monomials, coefficients, terms, expval

export @ncpolyvar, ς
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
