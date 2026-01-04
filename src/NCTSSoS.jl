module NCTSSoS

using SparseArrays, LinearAlgebra, JuMP
using CliqueTrees, ChordalGraph, Graphs
using CliqueTrees: EliminationAlgorithm, SupernodeType
import CliqueTrees.cliquetree

# ============================================================================
# Core Types (formerly FastPolynomials)
# ============================================================================

# Algebra types must come first (defines AlgebraType hierarchy)
include("types/algebra.jl")

# Variable registry (depends on algebra types)
include("types/registry.jl")

# Monomial type (depends on algebra types and registry)
include("types/monomial.jl")

# Term type (depends on monomial)
include("types/term.jl")

# Polynomial type (depends on term)
include("types/polynomial.jl")

# ============================================================================
# Simplification Algorithms (algebra-specific)
# ============================================================================

include("simplification/projector.jl")
include("simplification/unipotent.jl")
include("simplification/noncommutative.jl")
include("simplification/pauli.jl")
include("simplification/fermionic.jl")
include("simplification/bosonic.jl")

# ============================================================================
# Composed Types (depends on Polynomial for type checking in simplify)
# ============================================================================

include("types/composed.jl")

# ============================================================================
# Algorithms (canonicalization, basis generation)
# ============================================================================

include("algorithms/canonicalization.jl")
include("algorithms/basis.jl")

# ============================================================================
# State Polynomial Types
# ============================================================================

include("states/types.jl")
include("states/word.jl")
include("states/polynomial.jl")

# ============================================================================
# Utility Functions
# ============================================================================

include("util/helpers.jl")

# ============================================================================
# Optimization Framework
# ============================================================================

include("optimization/problem.jl")
include("optimization/elimination.jl")
include("optimization/sparsity.jl")
include("optimization/moment.jl")
include("optimization/sos.jl")
include("optimization/gns.jl")
include("optimization/interface.jl")

# ============================================================================
# Exports - User-Facing API Only
# ============================================================================

# Problem Definition
export PolyOpt, polyopt, PolyOptResult, SolverConfig

# Solver Interface
export cs_nctssos, cs_nctssos_higher, reconstruct

# Elimination Strategies
export NoElimination, MF, MMD, AsIsElimination, MaximalElimination

# Variable Creation (primary user entry point)
export create_pauli_variables, create_fermionic_variables, create_bosonic_variables
export create_projector_variables, create_unipotent_variables, create_noncommutative_variables

# Core Types (users need these for type annotations and construction)
export Monomial, Polynomial, Term, VariableRegistry
export ComposedMonomial

# Algebra Types (users need for dispatch)
export AlgebraType
export NonCommutativeAlgebra, PauliAlgebra, FermionicAlgebra
export BosonicAlgebra, ProjectorAlgebra, UnipotentAlgebra

# State Polynomial Operations
export ς, tr
export StateSymbol, StateWord, StatePolynomial
export NCStateWord, MaxEntangled, Arbitrary

# Polynomial Operations (commonly used)
export degree, monomials, coefficients, terms, variables
export simplify, simplify!
export variable_indices

# Canonicalization (user-facing)
export symmetric_canon, cyclic_canon, canonicalize

# Basis Generation
export get_ncbasis, get_state_basis

# Registry Helpers
export symbols, indices, subregistry

# Fermionic Helper
export has_even_parity

# Coefficient Type
export coeff_type

end
