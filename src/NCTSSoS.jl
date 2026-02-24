"""
    NCTSSoS

Non-Commutative Term-Sparsity Sum-of-Squares optimization in Julia.

This package provides tools for polynomial optimization over non-commutative
algebras using semidefinite programming (SDP) relaxations with term sparsity
exploitation.

# Supported Algebras

NCTSSoS supports multiple algebra types organized by their normal form structure:

- **MonoidAlgebra**: Normal form is a single monomial
  - `NonCommutativeAlgebra`: Free algebra with no simplification rules
  - `ProjectorAlgebra`: Projectors satisfying P² = P
  - `UnipotentAlgebra`: Involutory operators satisfying U² = 𝟙

- **TwistedGroupAlgebra**: Normal form is scalar × monomial
  - `PauliAlgebra`: Pauli spin operators with σ² = I and cyclic products

- **PBWAlgebra**: Normal form is a sum of monomials (Poincaré-Birkhoff-Witt)
  - `FermionicAlgebra`: Fermionic creation/annihilation operators
  - `BosonicAlgebra`: Bosonic creation/annihilation operators

# Quick Start

```julia
using NCTSSoS

# Create Pauli variables on 2 sites
reg, (σx, σy, σz) = create_pauli_variables(1:2)

# Define an objective (auto-simplified by PauliAlgebra)
obj = 0.25 * (σx[1] * σx[2] + σy[1] * σy[2] + σz[1] * σz[2])

# Build the optimization problem (add constraints as needed)
pop = polyopt(obj, reg)
```

See also: [`polyopt`](@ref), [`SolverConfig`](@ref), [`cs_nctssos`](@ref)
"""
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

# Polynomial type (depends on monomial)
include("types/polynomial.jl")

# ============================================================================
# Simplification Algorithms (algebra-specific)
# ============================================================================

include("simplification/site_helpers.jl")
include("simplification/projector.jl")
include("simplification/unipotent.jl")
include("simplification/noncommutative.jl")
include("simplification/pauli.jl")
include("simplification/fermionic.jl")
include("simplification/bosonic.jl")

# ============================================================================
# Arithmetic Operations (depends on simplification for NormalMonomial × NormalMonomial)
# ============================================================================

include("types/arithmetic.jl")

# ============================================================================
# Composed Types (depends on Polynomial for type checking in simplify)
# ============================================================================

include("types/composed.jl")

# ============================================================================
# Algorithms (canonicalization, basis generation)
# ============================================================================

# State type tags are needed by canonicalization (defines Arbitrary/MaxEntangled)
include("states/types.jl")

include("algorithms/canonicalization.jl")
include("algorithms/basis.jl")

# ============================================================================
# State Polynomial Types
# ============================================================================

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
export SparsityResult, compute_sparsity

# Solver Interface
export cs_nctssos, cs_nctssos_higher, reconstruct

# Elimination Strategies
export NoElimination, MF, MMD, AsIsElimination, MaximalElimination

# Variable Creation (primary user entry point)
export create_pauli_variables, create_fermionic_variables, create_bosonic_variables
export create_projector_variables, create_unipotent_variables, create_noncommutative_variables

# Core Types (users need these for type annotations and construction)
export NormalMonomial, Polynomial, VariableRegistry
export ComposedMonomial
export AbstractTensorMonomial, AbstractMonomial

# Algebra Types (users need for dispatch)
export AlgebraType
export MonoidAlgebra, TwistedGroupAlgebra, PBWAlgebra
export NonCommutativeAlgebra, PauliAlgebra, FermionicAlgebra
export BosonicAlgebra, ProjectorAlgebra, UnipotentAlgebra

# State Polynomial Operations
export ς, varsigma, tr
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
