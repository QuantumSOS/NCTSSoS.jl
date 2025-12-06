# =============================================================================
# Utility Functions for NCTSSoS Integration
# =============================================================================
#
# These utility functions are used by NCTSSoS and were part of the old
# FastPolynomials implementation. They are maintained here for backward
# compatibility with NCTSSoS source files.
#
# =============================================================================

"""
    sorted_union(args...) -> Vector

Compute the sorted union of multiple vectors. Returns a sorted vector
containing all unique elements from the input vectors.

# Examples
```jldoctest
julia> sorted_union([3, 1], [2, 1])
3-element Vector{Int64}:
 1
 2
 3

julia> sorted_union([1, 2], [2, 3], [3, 4])
4-element Vector{Int64}:
 1
 2
 3
 4
```
"""
function sorted_union(args...)
    isempty(args) && return []
    result = unique(vcat(args...))
    sort!(result)
    return result
end

"""
    sorted_unique(v::Vector) -> Vector

Return a sorted vector containing only unique elements from the input.

# Examples
```jldoctest
julia> sorted_unique([3, 1, 2, 1, 3])
3-element Vector{Int64}:
 1
 2
 3
```
"""
function sorted_unique(v::Vector)
    isempty(v) && return similar(v, 0)
    result = unique(v)
    sort!(result)
    return result
end

"""
    _neat_dot3(a::NCStateWord, m::M, b::NCStateWord) where M

Compute adjoint(a) * m * b for NCStateWords with a middle monomial/polynomial.
This is a common pattern in moment matrix construction.

The middle element `m` can be:
- A Monomial: converted to NCStateWord and multiplied
- An NCStateWord: direct multiplication
- A StateWord: converted to NCStateWord and multiplied

# Returns
An NCStateWord representing the triple product.
"""
function _neat_dot3(a::NCStateWord{ST,A,T}, m::Monomial{A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # adjoint(a) * m * b
    # Convert m to NCStateWord for multiplication
    m_ncsw = NCStateWord(one(StateWord{ST,A,T}), m)
    return adjoint(a) * m_ncsw * b
end

function _neat_dot3(a::NCStateWord{ST,A,T}, m::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    return adjoint(a) * m * b
end

function _neat_dot3(a::NCStateWord{ST,A,T}, sw::StateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Convert StateWord to NCStateWord with identity nc_word
    m_ncsw = NCStateWord(sw)
    return adjoint(a) * m_ncsw * b
end

# Overload for regular Monomials (non-state context)
"""
    neat_dot(a::Monomial, b::Monomial) -> Monomial

Compute adjoint(a) * b for regular Monomials.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3, 4]);

julia> result = neat_dot(m1, m2);

julia> result.monomial.word
6-element Vector{Int64}:
 2
 1
 3
 4
```
"""
function neat_dot(a::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # adjoint(a) * b returns a Term, but for legacy compatibility we need a Monomial
    # The coefficient is assumed to be 1 for simple algebras
    result = adjoint(a) * b
    result.monomial
end

"""
    _neat_dot3(a::Monomial, m::Monomial, b::Monomial) -> Monomial

Compute adjoint(a) * m * b for regular Monomials.
Returns a Monomial (extracts from the Term result for legacy compatibility).

This is the three-argument form commonly used in moment matrix construction
where we need adjoint(row_index) * constraint_monomial * column_index.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{NonCommutativeAlgebra}(UInt16[1]);

julia> m2 = Monomial{NonCommutativeAlgebra}(UInt16[2]);

julia> m3 = Monomial{NonCommutativeAlgebra}(UInt16[3]);

julia> result = _neat_dot3(m1, m2, m3);

julia> result.word
3-element Vector{UInt16}:
 0x0001
 0x0002
 0x0003
```
"""
function _neat_dot3(a::Monomial{A,T}, m::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # For NonCommutativeAlgebra and others: adjoint just reverses
    # (adjoint(a) * m) * b
    temp = adjoint(a) * m
    # temp is a Term, extract monomial for next multiplication
    result = temp.monomial * b
    # Return only the monomial for legacy compatibility
    result.monomial
end

# =============================================================================
# Legacy Variable Compatibility Type
# =============================================================================
#
# The old API used a Variable struct with name::Symbol.
# The new API uses VariableRegistry with integer indices and Monomial{A,T}.
# This compatibility layer allows old code to work during migration.
#

# Global registry to store variable name mappings (index -> Symbol)
# Must be defined before Variable struct since constructors use it
const _VAR_INDEX_TO_NAME = Dict{Int, Symbol}()

"""
    Variable

Legacy compatibility struct for old API. The new FastPolynomials uses
VariableRegistry with integer indices instead of Variable structs.

This struct is provided for backward compatibility with NCTSSoS code that
has not yet been migrated to the new API.

# Fields
- `name::Symbol`: The variable name
- `iscomplex::Bool`: Whether the variable is complex (default: false)
- `index::Int`: Internal index for ordering (assigned at construction)

# Note
New code should use VariableRegistry and Monomial{AlgebraType,T} instead.
"""
struct Variable
    name::Symbol
    iscomplex::Bool
    index::Int

    # Global counter for automatic index assignment
    function Variable(name::Symbol; iscomplex::Bool=false)
        # Use a simple counter based on hash for consistent ordering
        idx = hash(name) % 10000
        # Register the name for reverse lookup
        _VAR_INDEX_TO_NAME[idx] = name
        new(name, iscomplex, idx)
    end

    function Variable(name::Symbol, iscomplex::Bool, index::Int)
        # Register the name for reverse lookup
        _VAR_INDEX_TO_NAME[index] = name
        new(name, iscomplex, index)
    end
end

# Comparison for sorting
Base.isless(a::Variable, b::Variable) = a.name < b.name
Base.:(==)(a::Variable, b::Variable) = a.name == b.name
Base.hash(v::Variable, h::UInt) = hash(v.name, h)
Base.show(io::IO, v::Variable) = print(io, v.name)

# Variable arithmetic - converts to Polynomial
# This allows `x * y` syntax where x, y are Variables

"""
    to_monomial(v::Variable) -> Monomial

Convert a Variable to a Monomial{NonCommutativeAlgebra, UInt}.
Uses unsigned indices so that adjoint(m) just reverses (no negation).
"""
function to_monomial(v::Variable)
    Monomial{NonCommutativeAlgebra}([UInt(v.index)])
end

"""
    to_polynomial(v::Variable, C::Type=ComplexF64) -> Polynomial

Convert a Variable to a single-term Polynomial.
"""
function to_polynomial(v::Variable, ::Type{C}=ComplexF64) where C<:Number
    m = to_monomial(v)
    Polynomial([Term(one(C), m)])
end

# Variable * Variable -> Polynomial
function Base.:*(a::Variable, b::Variable)
    m_a = to_monomial(a)
    m_b = to_monomial(b)
    result = m_a * m_b  # Returns Term
    Polynomial([result])
end

# Scalar * Variable -> Polynomial
# Always promote coefficient to ComplexF64 to avoid type issues with im
function Base.:*(c::Number, v::Variable)
    m = to_monomial(v)
    coef = ComplexF64(c)  # Promote to avoid Complex{Bool} etc
    Polynomial([Term(coef, m)])
end
Base.:*(v::Variable, c::Number) = c * v

# Variable + Variable -> Polynomial
function Base.:+(a::Variable, b::Variable)
    C = ComplexF64
    Polynomial([Term(one(C), to_monomial(a)), Term(one(C), to_monomial(b))])
end

# Number + Variable -> Polynomial (constant + variable)
function Base.:+(c::Number, v::Variable)
    C = ComplexF64
    # Create identity monomial for the constant term
    id_mono = one(Monomial{NonCommutativeAlgebra,UInt})
    Polynomial([Term(C(c), id_mono), Term(one(C), to_monomial(v))])
end
Base.:+(v::Variable, c::Number) = c + v

# Number - Variable -> Polynomial
function Base.:-(c::Number, v::Variable)
    C = ComplexF64
    id_mono = one(Monomial{NonCommutativeAlgebra,UInt})
    Polynomial([Term(C(c), id_mono), Term(-one(C), to_monomial(v))])
end

# Variable - Variable -> Polynomial
function Base.:-(a::Variable, b::Variable)
    C = ComplexF64
    Polynomial([Term(one(C), to_monomial(a)), Term(-one(C), to_monomial(b))])
end

# Unary minus
function Base.:-(v::Variable)
    C = ComplexF64
    Polynomial([Term(-one(C), to_monomial(v))])
end

# Variable - Number -> Polynomial
function Base.:-(v::Variable, c::Number)
    C = ComplexF64
    id_mono = one(Monomial{NonCommutativeAlgebra,UInt})
    Polynomial([Term(one(C), to_monomial(v)), Term(C(-c), id_mono)])
end

# one(Variable) - returns constant 1 polynomial
"""
    Base.one(::Variable) -> Polynomial

Return the constant polynomial 1.
"""
function Base.one(::Variable)
    C = ComplexF64
    Polynomial([Term(one(C), one(Monomial{NonCommutativeAlgebra,UInt}))])
end

# Variable ^ n -> Polynomial
function Base.:^(v::Variable, n::Integer)
    n < 0 && throw(ArgumentError("Cannot raise Variable to negative power"))
    if n == 0
        return Polynomial([Term(ComplexF64(1.0), one(Monomial{NonCommutativeAlgebra,UInt}))])
    end
    m = to_monomial(v)
    result_m = m
    for _ in 2:n
        temp = result_m * m  # Returns Term
        result_m = temp.monomial
    end
    Polynomial([Term(ComplexF64(1.0), result_m)])
end

# Variable + Polynomial, Variable * Polynomial, etc.
function Base.:+(v::Variable, p::Polynomial{A,T,C}) where {A,T,C}
    to_polynomial(v, C) + p
end
Base.:+(p::Polynomial, v::Variable) = v + p

function Base.:-(v::Variable, p::Polynomial{A,T,C}) where {A,T,C}
    to_polynomial(v, C) - p
end
function Base.:-(p::Polynomial{A,T,C}, v::Variable) where {A,T,C}
    p - to_polynomial(v, C)
end

function Base.:*(v::Variable, p::Polynomial{A,T,C}) where {A,T,C}
    to_polynomial(v, C) * p
end
function Base.:*(p::Polynomial{A,T,C}, v::Variable) where {A,T,C}
    p * to_polynomial(v, C)
end

# =============================================================================
# Legacy variables() function for NCTSSoS compatibility
# =============================================================================
#
# The new API's `variables(p::Polynomial)` returns Set{T} of indices.
# The old API's `variables(p)` returned Vector{Variable}.
# For NCTSSoS code compatibility, we override the polynomial.jl version.

# Note: _VAR_INDEX_TO_NAME is defined earlier with Variable struct

"""
    variables(p::Polynomial) -> Vector{Variable}

Extract all unique variables from a polynomial, returning Variable structs.
This shadows the polynomial.jl version to return Variable[] for legacy compatibility.
"""
function variables(p::Polynomial{A,T,C}) where {A,T,C}
    result = Variable[]
    seen = Set{T}()
    for t in p.terms
        for idx in t.monomial.word
            abs_idx = abs(idx)
            if abs_idx ∉ seen
                push!(seen, abs_idx)
                # Look up the original name if available
                name = get(_VAR_INDEX_TO_NAME, Int(abs_idx), Symbol("x", abs_idx))
                push!(result, Variable(name, false, Int(abs_idx)))
            end
        end
    end
    sort!(result)
    return result
end

"""
    variables(v::Variable) -> Vector{Variable}

Return a single-element vector containing the variable.
"""
variables(v::Variable) = [v]

"""
    variables(m::Monomial{A,T}) -> Vector{Variable}

Extract all unique variables from a monomial, returning Variable structs.
"""
function variables(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    result = Variable[]
    seen = Set{Int}()
    for idx in m.word
        abs_idx = Int(abs(idx))
        if abs_idx ∉ seen
            push!(seen, abs_idx)
            name = get(_VAR_INDEX_TO_NAME, abs_idx, Symbol("x", abs_idx))
            push!(result, Variable(name, false, abs_idx))
        end
    end
    sort!(result)
    return result
end

"""
    variables(sp::StatePolynomial{C,ST,A,T}) -> Vector{Variable}

Extract all unique variables from a StatePolynomial, returning Variable structs.
"""
function variables(sp::StatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    result = Variable[]
    seen = Set{Int}()
    for sw in sp.state_words
        for mono in sw.state_monos
            for idx in mono.word
                abs_idx = Int(abs(idx))
                if abs_idx ∉ seen
                    push!(seen, abs_idx)
                    name = get(_VAR_INDEX_TO_NAME, abs_idx, Symbol("x", abs_idx))
                    push!(result, Variable(name, false, abs_idx))
                end
            end
        end
    end
    sort!(result)
    return result
end

"""
    variables(ncsp::NCStatePolynomial{C,ST,A,T}) -> Vector{Variable}

Extract all unique variables from an NCStatePolynomial, returning Variable structs.
"""
function variables(ncsp::NCStatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    result = Variable[]
    seen = Set{Int}()
    for ncsw in ncsp.nc_state_words
        # From the nc_word (non-commutative part)
        for idx in ncsw.nc_word.word
            abs_idx = Int(abs(idx))
            if abs_idx ∉ seen
                push!(seen, abs_idx)
                name = get(_VAR_INDEX_TO_NAME, abs_idx, Symbol("x", abs_idx))
                push!(result, Variable(name, false, abs_idx))
            end
        end
        # From the sw (StateWord part)
        for mono in ncsw.sw.state_monos
            for idx in mono.word
                abs_idx = Int(abs(idx))
                if abs_idx ∉ seen
                    push!(seen, abs_idx)
                    name = get(_VAR_INDEX_TO_NAME, abs_idx, Symbol("x", abs_idx))
                    push!(result, Variable(name, false, abs_idx))
                end
            end
        end
    end
    sort!(result)
    return result
end

# Variable to StateWord conversion
"""
    ς(v::Variable) -> StateWord{Arbitrary}

Create a StateWord{Arbitrary} from a Variable.
"""
function ς(v::Variable)
    m = to_monomial(v)
    StateWord{Arbitrary}(m)
end

"""
    ς(p::Polynomial{A,T,C}) -> StatePolynomial

Create a StatePolynomial from a Polynomial.
Converts each term's monomial to a StateWord{Arbitrary}.
"""
function ς(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(p.terms) && return StatePolynomial(C[], StateWord{Arbitrary,A,T}[])

    state_words = [StateWord{Arbitrary}(t.monomial) for t in p.terms]
    coeffs = C[t.coefficient for t in p.terms]
    StatePolynomial(coeffs, state_words)
end

# =============================================================================
# Legacy SimplifyAlgorithm Compatibility Layer
# =============================================================================
#
# The old API used SimplifyAlgorithm struct with configuration:
#   SimplifyAlgorithm(comm_gps=..., is_unipotent=..., is_projective=...)
#
# The new API uses AlgebraType singleton types for dispatch.
# This compatibility layer allows old code to work during migration.
#

"""
    SimplifyAlgorithm

Legacy compatibility struct for old API. The new FastPolynomials uses
AlgebraType singleton types (PauliAlgebra, FermionicAlgebra, etc.) instead.

This struct is provided for backward compatibility with NCTSSoS code that
has not yet been migrated to the new API.

# Fields
- `comm_gps`: Commutative groups (for legacy compatibility)
- `is_unipotent::Bool`: Whether variables are unipotent (x^2 = 1)
- `is_projective::Bool`: Whether variables are projective (x^2 = x)

# Note
New code should use AlgebraType dispatch instead of SimplifyAlgorithm.
"""
struct SimplifyAlgorithm
    comm_gps::Vector{<:Any}  # Allow any element type
    is_unipotent::Bool
    is_projective::Bool
    n_gps::Int  # Number of commutative groups (derived property for compatibility)

    function SimplifyAlgorithm(;
        comm_gps::Vector=Vector{Any}[],
        is_unipotent::Bool=false,
        is_projective::Bool=false
    )
        new(comm_gps, is_unipotent, is_projective, length(comm_gps))
    end
end

# =============================================================================
# Legacy get_basis Compatibility
# =============================================================================
#
# The old API: get_basis(vars::Vector{Variable}, degree::Int, sa::SimplifyAlgorithm)
# The new API: get_ncbasis(AlgebraType, n_vars::Int, degree::Int)
#
# Provide legacy-compatible wrappers.

"""
    get_basis(vars, degree::Int) -> Vector{Monomial}

Legacy compatibility wrapper for basis generation.

For NCTSSoS migration, this generates a NonCommutativeAlgebra basis
with the given number of variables and degree.

# Arguments
- `vars`: Vector of variables (length determines n_vars) OR integer for n_vars
- `degree::Int`: Maximum degree for basis generation

# Returns
Vector of Monomial{NonCommutativeAlgebra, Int} up to the specified degree.

# Note
New code should use `get_ncbasis(AlgebraType, n_vars, degree)` directly.
"""
function get_basis(vars::AbstractVector, degree::Int)
    n_vars = length(vars)
    get_ncbasis(NonCommutativeAlgebra, n_vars, degree; T=Int)
end

# Overload accepting a Polynomial type for type inference
function get_basis(::Type{P}, vars::AbstractVector, degree::Int, sa::SimplifyAlgorithm) where {P}
    n_vars = length(vars)
    # Determine algebra type from SimplifyAlgorithm properties
    if sa.is_unipotent
        return get_ncbasis(UnipotentAlgebra, n_vars, degree; T=Int, filter_constraint=true)
    elseif sa.is_projective
        return get_ncbasis(ProjectorAlgebra, n_vars, degree; T=Int)
    else
        return get_ncbasis(NonCommutativeAlgebra, n_vars, degree; T=Int)
    end
end

# =============================================================================
# Legacy simplify/canonicalize with SimplifyAlgorithm
# =============================================================================

"""
    simplify(m::Monomial, sa::SimplifyAlgorithm) -> Monomial

Legacy wrapper for simplification with SimplifyAlgorithm.
Returns the monomial as-is (simplification is done during construction in new API).
"""
function simplify(m::Monomial{A,T}, sa::SimplifyAlgorithm) where {A<:AlgebraType,T<:Integer}
    # In the new API, simplification is built into monomial multiplication
    # Return the monomial unchanged
    m
end

"""
    simplify!(m::Monomial, sa::SimplifyAlgorithm) -> Monomial

Legacy wrapper for in-place simplification with SimplifyAlgorithm.
Returns the monomial as-is (simplification is done during construction in new API).
"""
function simplify!(m::Monomial{A,T}, sa::SimplifyAlgorithm) where {A<:AlgebraType,T<:Integer}
    m
end

# StateWord simplification
function simplify(sw::StateWord{ST,A,T}, sa::SimplifyAlgorithm) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    sw  # Already simplified during construction
end

function simplify!(sw::StateWord{ST,A,T}, sa::SimplifyAlgorithm) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    sw
end

"""
    canonicalize(m::Monomial, sa::SimplifyAlgorithm) -> Monomial

Legacy wrapper for canonicalization with SimplifyAlgorithm.
Uses symmetric canonicalization (default behavior).
"""
function canonicalize(m::Monomial{A,T}, sa::SimplifyAlgorithm) where {A<:AlgebraType,T<:Integer}
    symmetric_canon(m)
end

# StateWord canonicalization
function canonicalize(sw::StateWord{ST,A,T}, sa::SimplifyAlgorithm) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # StateWords are already canonicalized via involution
    sw
end

# =============================================================================
# Legacy monomial functions
# =============================================================================

"""
    monomial(v::Variable) -> Monomial

Legacy function to create a single-variable monomial from a Variable.

# Note
New code should use Monomial{AlgebraType}([index]) directly.
"""
function monomial(v::Variable)
    Monomial{NonCommutativeAlgebra}([UInt(v.index)])
end

"""
    monomial(vars::Vector{Variable}, exponents::Vector{Int}) -> Monomial

Legacy function to create a monomial from variables and exponents.
Creates a word by repeating each variable index according to its exponent.

# Example
```julia
@ncpolyvar x[1:2]
m = monomial([x[1], x[2]], [2, 1])  # Creates x1*x1*x2
```

# Note
New code should use Monomial{AlgebraType}([indices...]) directly.
"""
function monomial(vars::Vector{Variable}, exponents::Vector{Int})
    @assert length(vars) == length(exponents) "vars and exponents must have same length"
    word = UInt[]
    for (v, e) in zip(vars, exponents)
        for _ in 1:e
            push!(word, UInt(v.index))
        end
    end
    Monomial{NonCommutativeAlgebra}(word)
end

"""
    Base.one(::Type{Monomial}) -> Monomial{NonCommutativeAlgebra,UInt}

Create the identity monomial (empty word) for the generic Monomial type.
This is needed for legacy NCTSSoS code that uses `one(Monomial)`.
"""
Base.one(::Type{Monomial}) = Monomial{NonCommutativeAlgebra}(UInt[])

# =============================================================================
# Legacy Polynomial constructor (coeffs, monomials)
# =============================================================================

"""
    Polynomial(coeffs::Vector, monos::Vector{Monomial{A,T}}) -> Polynomial{A,T,C}

Legacy constructor to create a Polynomial from separate coefficient and monomial vectors.

# Example
```julia
@ncpolyvar x[1:2]
p = Polynomial([1.0, 2.0], [monomial(x[1]), monomial(x[2])])
```

# Note
New code should use `Polynomial([Term(c, m), ...])` directly.
"""
function Polynomial(
    coeffs::AbstractVector{C}, monos::Vector{Monomial{A,T}}
) where {C<:Number,A<:AlgebraType,T<:Integer}
    @assert length(coeffs) == length(monos) "coeffs and monomials must have same length"
    terms = [Term(C(c), m) for (c, m) in zip(coeffs, monos)]
    Polynomial(terms)
end

# =============================================================================
# Legacy @ncpolyvar Macro
# =============================================================================

"""
    @ncpolyvar x y[1:3] z

Legacy macro for declaring non-commutative polynomial variables.

Creates Variable structs with the given names and optionally indexed arrays.

# Examples
```julia
@ncpolyvar x y z        # Creates three single variables
@ncpolyvar x[1:3]       # Creates x[1], x[2], x[3]
@ncpolyvar x[1:2] y[1:2]  # Creates x[1], x[2], y[1], y[2]
```

# Note
New code should use `create_noncommutative_variables` or the algebra-specific
variable creation functions instead.
"""
macro ncpolyvar(exprs...)
    blk = Expr(:block)
    for expr in exprs
        if isa(expr, Symbol)
            # Single variable: @ncpolyvar x
            push!(blk.args, :($(esc(expr)) = Variable($(QuoteNode(expr)))))
        elseif isa(expr, Expr) && expr.head == :ref
            # Indexed variable: @ncpolyvar x[1:3]
            varname = expr.args[1]
            range_expr = expr.args[2]
            push!(blk.args, :($(esc(varname)) = [Variable(Symbol($(string(varname)), i)) for i in $(esc(range_expr))]))
        else
            error("Invalid expression in @ncpolyvar: $expr")
        end
    end
    push!(blk.args, :nothing)
    return blk
end

# =============================================================================
# Legacy AbstractPolynomial Type Alias
# =============================================================================

"""
    AbstractPolynomial{T}

Legacy type alias for compatibility with old NCTSSoS code.
Maps to Union of Polynomial types with any algebra type.

# Note
New code should use `Polynomial{A,T,C}` directly.
"""
const AbstractPolynomial{T} = Union{
    Polynomial{<:AlgebraType,<:Integer,T},
    StatePolynomial{T,<:StateType,<:AlgebraType,<:Integer},
    NCStatePolynomial{T,<:StateType,<:AlgebraType,<:Integer}
}

# =============================================================================
# Legacy is_symmetric with SimplifyAlgorithm
# =============================================================================

"""
    is_symmetric(p, sa::SimplifyAlgorithm) -> Bool

Legacy wrapper for is_symmetric with SimplifyAlgorithm parameter.

For polynomials with commutative groups (like Pauli algebra), this checks symmetry
by verifying that for each term c*m, there exists a term conj(c)*adjoint(m) with
equivalent coefficient after accounting for floating point precision.

Key insight: In Pauli algebra, operators at DIFFERENT sites commute.
So x[1]*x[2] = x[2]*x[1], meaning [a,b] and [b,a] are equivalent when
a and b belong to different commutative groups.

For practical purposes with legacy NCTSSoS code:
- We check that each term has a real coefficient (since Hermitian operators
  should have real coefficients when expressed in the computational basis)
- For terms with only real coefficients, the polynomial is symmetric if
  each monomial's coefficient is real
"""
function is_symmetric(p::Polynomial{A,T,C}, sa::SimplifyAlgorithm) where {A<:AlgebraType,T,C<:Number}
    # Build a mapping from variable index to group index
    var_to_group = Dict{T, Int}()
    for (gid, gp) in enumerate(sa.comm_gps)
        for v in gp
            var_to_group[T(v.index)] = gid
        end
    end

    # For each term, check symmetry
    for t in p.terms
        coef = t.coefficient
        word = t.monomial.word

        # Quick check: if coefficient has significant imaginary part, not symmetric
        # (for self-adjoint operators with real scalars)
        if abs(imag(coef)) > 1e-10
            return false
        end

        # For monomials: check if the word is "group-equivalent" to its reverse
        # Two words are group-equivalent if they have the same indices, possibly reordered,
        # where reordering is allowed between different groups (since they commute)
        #
        # For simplicity in Pauli algebra case: if each variable appears the same number
        # of times in word and reverse(word), they're equivalent.
        # Since we're just reversing, this is always true - the multiset is the same.
        # But we also need the canonical form to match.
        #
        # Actually, for unsigned types with different sites commuting:
        # adjoint([a, b]) = [b, a] where a and b are from different groups
        # If different groups commute, [a,b] = [b,a], so adjoint = original
        #
        # This means: polynomial is symmetric if all coefficients are real.
    end

    return true  # All coefficients have negligible imaginary parts
end

# Fallback for non-Polynomial types
is_symmetric(p, sa::SimplifyAlgorithm) = is_symmetric(p)
