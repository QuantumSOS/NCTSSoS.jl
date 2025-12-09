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

Compute adjoint(a) * b for regular Monomials by concatenating words.

Returns a Monomial with the adjoint of a's word followed by b's word.
Does NOT apply simplification - callers should simplify! explicitly if needed.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3, 4]);

julia> result = neat_dot(m1, m2);

julia> result.word
4-element Vector{Int64}:
 -2
 -1
  3
  4
```
"""
function neat_dot(a::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # adjoint(a) * b - concatenate adjoint(a).word with b.word
    adjoint(a) * b
end

"""
    _neat_dot3(a::Monomial, m::Monomial, b::Monomial) -> Monomial

Compute adjoint(a) * m * b for regular Monomials by concatenating words.

Returns a Monomial with adjoint(a).word concatenated with m.word and b.word.
Does NOT apply simplification - callers should simplify! explicitly if needed.

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
    # adjoint(a) * m * b - just concatenate the words
    adjoint(a) * m * b
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
    result = m_a * m_b  # Returns Monomial (word concatenation)
    Polynomial([Term(one(ComplexF64), result)])
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
        result_m = result_m * m  # Returns Monomial (word concatenation)
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

# NOTE: We intentionally do NOT override variables() for StatePolynomial here.
# StatePolynomial tests use integer indices directly, so they need the original
# implementation from state_polynomial.jl which returns Set{T}.

# However, NCStatePolynomial is used by NCTSSoS via polyopt() which expects
# Vector{Variable} for compatibility with comm_gps. So we override it here.

"""
    variables(ncsp::NCStatePolynomial{C,ST,A,T}) -> Vector{Variable}

Extract all unique variables from an NCStatePolynomial, returning Variable structs.
This is needed for NCTSSoS polyopt() compatibility with comm_gps parameter.
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
# Legacy get_basis Compatibility
# =============================================================================
#
# The old API: get_basis(vars::Vector{Variable}, degree::Int)
# The new API: get_ncbasis(VariableRegistry, degree::Int)
#
# Provide legacy-compatible wrappers that generate bases directly.

"""
    get_basis(vars, degree::Int) -> Vector{Monomial}

Legacy compatibility wrapper for basis generation.

For NCTSSoS migration, this generates a NonCommutativeAlgebra basis
with the given number of variables and degree.

# Arguments
- `vars`: Vector of variables
- `degree::Int`: Maximum degree for basis generation

# Returns
Vector of Monomial{NonCommutativeAlgebra, Int} up to the specified degree.

# Note
New code should use `get_ncbasis(VariableRegistry, degree)` directly.
"""
function get_basis(vars::AbstractVector{Variable}, degree::Int)
    # Extract indices from Variable structs
    indices = [v.index for v in vars]

    # Generate all words (monomial representations) up to the given degree
    result = Monomial{NonCommutativeAlgebra,Int}[]

    # Helper function to generate all words of a specific degree
    function generate_words(current_word::Vector{Int}, remaining_deg::Int)
        if remaining_deg == 0
            push!(result, Monomial{NonCommutativeAlgebra}(current_word))
            return
        end
        for idx in indices
            generate_words(vcat(current_word, idx), remaining_deg - 1)
        end
    end

    # Generate words for each degree from 0 to d
    for d in 0:degree
        if d == 0
            push!(result, Monomial{NonCommutativeAlgebra}(Int[]))  # Identity monomial
        else
            generate_words(Int[], d)
        end
    end

    return result
end

# Fallback for non-Variable vectors (e.g., Monomial vectors)
function get_basis(vars::AbstractVector, degree::Int)
    # If vars are Monomials, extract unique variable indices from them
    indices = Int[]
    for v in vars
        if v isa Variable
            push!(indices, v.index)
        elseif v isa Monomial
            for idx in v.word
                push!(indices, Int(abs(idx)))
            end
        end
    end
    indices = sort(unique(indices))

    # Generate all words up to the given degree
    result = Monomial{NonCommutativeAlgebra,Int}[]

    function generate_words(current_word::Vector{Int}, remaining_deg::Int)
        if remaining_deg == 0
            push!(result, Monomial{NonCommutativeAlgebra}(current_word))
            return
        end
        for idx in indices
            generate_words(vcat(current_word, idx), remaining_deg - 1)
        end
    end

    for d in 0:degree
        if d == 0
            push!(result, Monomial{NonCommutativeAlgebra}(Int[]))
        else
            generate_words(Int[], d)
        end
    end

    return result
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
