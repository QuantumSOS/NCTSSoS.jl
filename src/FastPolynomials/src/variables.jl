"""
    Variable

A `Variable` represents a symbolic variable in a polynomial expression.

# Fields
- `name::Symbol`: The name of the variable
- `iscomplex::Bool`: Whether the variable is complex
"""
struct Variable
    name::Symbol
    iscomplex::Bool

    function Variable(name::Symbol; iscomplex::Bool=false)
        return new(name, iscomplex)
    end
end

"""
    polyarrayvar(prefix, indices...; iscomplex=false)

Creates an array of variables with indexed names.

# Arguments
- `prefix::String`: The base name for the variables
- `indices...`: Variable number of index ranges like `i1, i2, ...`

# Keyword Arguments
- `iscomplex::Bool`: Whether the variables are complex

# Returns
- Array of `Variable` objects with names formatted as `prefix[i1,i2,...]`
"""
function polyarrayvar(prefix, indices...; iscomplex=false)
    subscripts = ("₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉")
    return map(
        idxs -> Variable(
            Symbol(
                prefix,
                join(map(i -> join(reverse(subscripts[digits(i) .+ 1])), idxs), ","),
            );
            iscomplex=iscomplex,
        ),
        Iterators.product(indices...),
    )
end

"""
    @ncpolyvar(args...)

Macro to create non-commutative polynomial variables (hermitian operators).

# Arguments
- `args...`: Variable arguments specifying variable names (symbols or indexed expressions)

# Returns
- Tuple of created Variable objects

# Example
```jldoctest; setup=:(using NCTSSoS.FastPolynomials)
julia> @ncpolyvar x y z     # Creates three real variables
(x, y, z)

julia> @ncpolyvar u[1:3]     # Creates array of variables u[1], u[2], u[3]
(Variable[u₁, u₂, u₃],)
```

"""
macro ncpolyvar(args...)
    vars, exprs = _buildpolyvars(args, false)
    # calls the exprs that initializes the variables
    return :($(foldl((x, y) -> :($x; $y), exprs; init=:()));
    $(Expr(:tuple, esc.(vars)...))) # returns the variables
end

# Builds multiple polynomial variable declarations from a collection of arguments.
# Arguments
# - `args`: Collection of symbols or expressions for variable names
# - `iscomplex::Bool`: Whether the variables are complex
# Returns
# - Tuple of (variable_names_array, expressions_array) for variable creation
function _buildpolyvars(args, iscomplex)
    vars = Symbol[]
    exprs = []
    for arg in args
        var, expr = buildpolyvar(arg, iscomplex)
        push!(vars, var)
        push!(exprs, expr)
    end
    return vars, exprs
end

"""
    buildpolyvar(var, iscomplex)

Builds a polynomial variable declaration from a symbol or expression.

# Arguments
- `var`: Either a Symbol for single variable or Expr for array variables
- `iscomplex::Bool`: Whether the variables are complex

# Returns
- Tuple of (variable_name, expression) for variable creation
"""
function buildpolyvar(var, iscomplex)
    @match var begin
        ::Symbol =>
            (var, :($(esc(var)) = $Variable($(QuoteNode(var)); iscomplex=$iscomplex)))
        :($varname[$(indices...)]) => begin
            (
                varname,
                :(
                    $(esc(varname)) = $polyarrayvar(
                        $(QuoteNode(varname)), $(esc.(indices)...); iscomplex=$iscomplex
                    )
                ),
            )
        end
        _ => error("Expected $var to be a variable name")
    end
end

# from https://scientificcoder.com/user-defined-show-method-in-julia
# 2-argument show, used by Array show, print(obj) and repr(obj), keep it short
function Base.show(io::IO, obj::Variable)
    return print_object(io, obj; multiline=true)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::Variable)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

function print_object(io::IO, obj::Variable; multiline::Bool)
    return multiline ? print(io, "$(obj.name)") : Base.show_default(io, obj)
end

function Base.hash(v::Variable, u::UInt)
    return hash(v.name, u)
end

Base.one(::Variable) = Monomial(Variable[], Int[])

"""
    monomials(vars::Vector{Variable}, cur_d::Int)

Generates all monomials of a specific degree from given variables.

# Arguments
- `vars::Vector{Variable}`: Variables to use in monomials
- `cur_d::Int`: Degree of monomials to generate

# Returns
- `Vector{Monomial}`: All monomials of degree `cur_d` in the given variables
"""
function monomials(vars::Vector{Variable}, ::Val{D}) where {D}
    return map(cartesian_product(vars, Val{D}())) do cur_vars
        monomial(collect(Variable, cur_vars), ones(Int, length(cur_vars)))
    end
end

function cartesian_product(vars::Vector{T}, ::Val{D}) where {T,D}
    iterable = Iterators.product(ntuple(_ -> vars, Val{D}())...)
    a = collect(iterable)
    return reshape(a, :)
end

"""
    get_basis(vars::Vector{Variable}, d::Int)

Generates a sorted basis of all monomials up to a given degree.

# Arguments
- `vars::Vector{Variable}`: Variables to use in the basis
- `d::Int`: Maximum degree of monomials

# Returns
- `Vector{Monomial}`: Sorted basis containing all monomials of degrees `0` through `d`
"""
function get_basis(vars::Vector{Variable}, d::Int)
    vars = sort(vars)
    d == 0 && return [Monomial(Variable[], Int[])]

    num_vars = length(vars)
    basis_length = sum(num_vars^i for i in 0:d)
    vars_vec = [Variable[]]
    sizehint!(vars_vec, basis_length)
    expos_vec = [Int[]]
    sizehint!(expos_vec, basis_length)

    for v in vars 
        push!(vars_vec, [v])
        push!(expos_vec, [1])
    end

    d == 1 && return Monomial.(vars_vec, expos_vec) 

    last_end_idx = 1

    @inbounds for i in 2:d
        # index of last monomial
        last_end_idx += num_vars^(i - 1)
        for j in 1:num_vars
            for k in 1:(num_vars^(i - 1))
                if vars_vec[last_end_idx - k + 1][end] == vars[j]
                    push!(vars_vec, copy(vars_vec[last_end_idx - k + 1]))
                    push!(expos_vec, copy(expos_vec[last_end_idx - k + 1]))
                    expos_vec[end][end] += 1
                else
                    push!(vars_vec, [vars_vec[last_end_idx - k + 1]; vars[j]])
                    push!(expos_vec, [expos_vec[last_end_idx - k + 1]; 1])
                end
            end
        end
    end
    # technically speaking, if I constructed it correctly, it should already be sorted, could be related to order of vars
    return Monomial.(vars_vec, expos_vec)
end

function Base.:(^)(a::Variable, expo::Int)
    @assert expo >= 0 "Exponent must be non-negative."
    return iszero(expo) ? one(a) : Monomial([a], [expo])
end
