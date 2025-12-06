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
        new(name, iscomplex, idx)
    end

    function Variable(name::Symbol, iscomplex::Bool, index::Int)
        new(name, iscomplex, index)
    end
end

# Comparison for sorting
Base.isless(a::Variable, b::Variable) = a.name < b.name
Base.:(==)(a::Variable, b::Variable) = a.name == b.name
Base.hash(v::Variable, h::UInt) = hash(v.name, h)
Base.show(io::IO, v::Variable) = print(io, v.name)

"""
    variables(p::Polynomial) -> Vector{Variable}

Extract the variables used in a polynomial as legacy Variable structs.
This is for backward compatibility - new code should work with integer indices.
"""
function variables(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Get unique variable indices from the polynomial
    var_indices = Set{T}()
    for t in p.terms
        for idx in t.monomial.word
            push!(var_indices, abs(idx))  # abs for fermionic/bosonic
        end
    end

    # Convert to sorted Vector of Variables with synthetic names
    sorted_indices = sort(collect(var_indices))
    return [Variable(Symbol("x_", i), false, i) for i in sorted_indices]
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

    function SimplifyAlgorithm(;
        comm_gps::Vector=Vector{Any}[],
        is_unipotent::Bool=false,
        is_projective::Bool=false
    )
        new(comm_gps, is_unipotent, is_projective)
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
# Legacy monomial(var) function
# =============================================================================

"""
    monomial(v::Variable) -> Monomial

Legacy function to create a single-variable monomial from a Variable.

# Note
New code should use Monomial{AlgebraType}([index]) directly.
"""
function monomial(v::Variable)
    Monomial{NonCommutativeAlgebra}([v.index])
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
