module FastPolynomials

# Core types
include("algebra_types.jl")
include("variable_registry.jl")
include("monomial.jl")
include("term.jl")

# Simplification algorithms
include("simplification/projector.jl")
include("simplification/unipotent.jl")
include("simplification/noncommutative.jl")
include("simplification/pauli.jl")
include("simplification/fermionic.jl")
include("simplification/bosonic.jl")

# Composed types
include("composed_monomial.jl")

# Polynomial type
include("polynomial.jl")

# Canonicalization algorithms
include("canonicalization.jl")

# Basis generation
include("basis.jl")

# State polynomial types
include("state_types.jl")
include("state_word.jl")
include("state_polynomial.jl")

# Exports
export VariableRegistry
export AbstractMonomial, Monomial, ComposedMonomial
export Term
export Polynomial
export degree, maxdegree
export coefficients, monomials, terms, support, variables
export simplify, simplify!
export star, star!, adjoint!
export symmetric_canon, cyclic_canon, cyclic_symmetric_canon, canonicalize

# Basis generation
export get_ncbasis, get_ncbasis_deg, has_consecutive_repeats

# Algebra types (singleton types for dispatch)
export AlgebraType
export NonCommutativeAlgebra, PauliAlgebra, FermionicAlgebra, BosonicAlgebra
export ProjectorAlgebra, UnipotentAlgebra

# State types
export StateType, Arbitrary, MaxEntangled
export StateWord, NCStateWord
export StatePolynomial, NCStatePolynomial
export expval, neat_dot
export ς, tr

# Polynomial utilities
export is_symmetric

# Variable creation functions
export create_pauli_variables, create_fermionic_variables, create_bosonic_variables
export create_projector_variables, create_unipotent_variables, create_noncommutative_variables

# Variable registry helpers
export symbols, indices

end
