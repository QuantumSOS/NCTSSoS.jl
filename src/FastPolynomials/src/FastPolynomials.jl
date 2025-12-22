module FastPolynomials

# Core types
include("algebra_types.jl")
include("variable_registry.jl")
include("monomial.jl")
include("term.jl")

# Polynomial type (needed by simplification algorithms and composed_monomial)
include("polynomial.jl")

# Simplification algorithms
include("simplification/projector.jl")
include("simplification/unipotent.jl")
include("simplification/noncommutative.jl")
include("simplification/pauli.jl")
include("simplification/fermionic.jl")
include("simplification/bosonic.jl")

# Composed types (depends on Polynomial for type checking in simplify)
include("composed_monomial.jl")

# Canonicalization algorithms
include("canonicalization.jl")

# Basis generation
include("basis.jl")

# State polynomial types
include("state_types.jl")
include("state_word.jl")
include("state_polynomial.jl")

# Utility functions (including legacy compatibility)
include("utils.jl")

# Exports
export VariableRegistry
export AbstractMonomial, Monomial, ComposedMonomial
export Term
export AbstractPolynomial, Polynomial
export degree
export coefficients, monomials, terms, variables
export simplify, simplify!
export symmetric_canon, cyclic_canon, cyclic_symmetric_canon, canonicalize

# Basis generation
export get_ncbasis, get_ncbasis_deg, get_state_basis

# Algebra types (singleton types for dispatch)
export AlgebraType
export NonCommutativeAlgebra, PauliAlgebra, FermionicAlgebra, BosonicAlgebra
export ProjectorAlgebra, UnipotentAlgebra
export default_coeff_type

# Fermionic algebra helpers
export has_even_parity

# State types
export StateType, Arbitrary, MaxEntangled
export StateWord, NCStateWord
export StatePolynomial, NCStatePolynomial
export expval, neat_dot
export ς, tr



# Variable creation functions
export create_pauli_variables, create_fermionic_variables, create_bosonic_variables
export create_projector_variables, create_unipotent_variables, create_noncommutative_variables

# Variable registry helpers
export symbols, indices, subregistry

# Polynomial helpers
export variable_indices


# Utility functions (kept from legacy API)
export sorted_union, sorted_unique
export _neat_dot3

end
