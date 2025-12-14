module NCTSSoS

using SparseArrays, LinearAlgebra, JuMP
using CliqueTrees, ChordalGraph, Graphs
using CliqueTrees: EliminationAlgorithm, SupernodeType
import CliqueTrees.cliquetree

include("FastPolynomials/src/FastPolynomials.jl")
using .FastPolynomials

# Import types (core types for internal use)
using .FastPolynomials: Monomial, Polynomial
using .FastPolynomials: AlgebraType, NonCommutativeAlgebra, PauliAlgebra
using .FastPolynomials: UnipotentAlgebra, ProjectorAlgebra, FermionicAlgebra, BosonicAlgebra
using .FastPolynomials: VariableRegistry, symbols, indices, subregistry
using .FastPolynomials: variable_indices
using .FastPolynomials: create_pauli_variables, create_fermionic_variables, create_bosonic_variables
using .FastPolynomials: create_projector_variables, create_unipotent_variables, create_noncommutative_variables

# Import type alias for backward compatibility
using .FastPolynomials: AbstractPolynomial

# Import functions (only those actually used in NCTSSoS)
using .FastPolynomials: sorted_union, sorted_unique
using .FastPolynomials: monomials, coefficients, terms, maxdegree
using .FastPolynomials: get_ncbasis
using .FastPolynomials: symmetric_canon
using .FastPolynomials: neat_dot, _neat_dot3, expval

# Re-export key types and functions for downstream users
export ς

# Note: Core types (Monomial, Polynomial, Term)
# are NOT exported from NCTSSoS to avoid ambiguity when FastPolynomials
# is also loaded. Access them via NCTSSoS.FastPolynomials.Monomial etc.
# or using NCTSSoS.FastPolynomials: Monomial, Polynomial, ...

export polyopt, cpolyopt
export SolverConfig
export NoElimination, MF, MMD, AsIsElimination, MaximalElimination
export cs_nctssos, cs_nctssos_higher
export reconstruct

# Algebra constructors (registry-based)
export pauli_algebra, fermionic_algebra, bosonic_algebra
export projector_algebra, unipotent_algebra, noncommutative_algebra

# Variable creation functions (public API)
export create_pauli_variables, create_fermionic_variables, create_bosonic_variables
export create_projector_variables, create_unipotent_variables, create_noncommutative_variables

include("pop.jl")

include("elimination.jl")

include("sparse.jl")

include("moment_solver.jl")

# NOTE: complex_moment_solver.jl has been merged into moment_solver.jl
# The unified MomentProblem{A,T,M,P} handles both real and complex algebras.
# complex_moment_relax() is now a backward-compatibility alias in moment_solver.jl.

include("sos_solver.jl")

include("gns.jl")

include("interface.jl")

include("algebra_constructors.jl")
end
