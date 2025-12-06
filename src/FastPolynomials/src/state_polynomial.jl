"""
    StatePolynomial{C<:Number, ST<:StateType, A<:AlgebraType, T<:Integer}

A polynomial in state words with coefficients.

Represents a sum of state words with coefficients: sum_i c_i * sw_i

# Fields
- `coeffs::Vector{C}`: Coefficients for each state word
- `state_words::Vector{StateWord{ST,A,T}}`: Sorted unique state words

# Invariants
- `state_words` is sorted by the StateWord ordering
- All state words are unique (combined during construction)
- All coefficients are non-zero (zeros removed during construction)
- `length(coeffs) == length(state_words)`

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3]);

julia> sw1 = StateWord{Arbitrary}([m1]);

julia> sw2 = StateWord{Arbitrary}([m2]);

julia> sp = StatePolynomial([1.0, 2.0], [sw1, sw2]);

julia> length(sp.state_words)
2
```

See also: [`StateWord`](@ref), [`NCStatePolynomial`](@ref)
"""
struct StatePolynomial{C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    coeffs::Vector{C}
    state_words::Vector{StateWord{ST,A,T}}

    function StatePolynomial(
        coeffs::Vector{C}, state_words::Vector{StateWord{ST,A,T}}
    ) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
        # Handle empty case
        isempty(state_words) && return new{C,ST,A,T}(C[], StateWord{ST,A,T}[])

        # Sort, combine duplicates, filter zeros
        perm = sortperm(state_words)
        sorted_sws = state_words[perm]
        sorted_coeffs = coeffs[perm]

        result_coeffs = C[]
        result_sws = StateWord{ST,A,T}[]

        current_sw = sorted_sws[1]
        current_coef = sorted_coeffs[1]

        for i in 2:length(sorted_sws)
            if sorted_sws[i] == current_sw
                current_coef += sorted_coeffs[i]
            else
                if !iszero(current_coef)
                    push!(result_coeffs, current_coef)
                    push!(result_sws, current_sw)
                end
                current_sw = sorted_sws[i]
                current_coef = sorted_coeffs[i]
            end
        end
        if !iszero(current_coef)
            push!(result_coeffs, current_coef)
            push!(result_sws, current_sw)
        end

        new{C,ST,A,T}(result_coeffs, result_sws)
    end
end

# =============================================================================
# Accessors
# =============================================================================

"""
    coefficients(sp::StatePolynomial) -> Vector{C}

Get the coefficients of the StatePolynomial.
"""
coefficients(sp::StatePolynomial) = sp.coeffs

"""
    monomials(sp::StatePolynomial) -> Vector{StateWord}

Get the state words of the StatePolynomial (these are the "monomials" in state space).
"""
monomials(sp::StatePolynomial) = sp.state_words

"""
    terms(sp::StatePolynomial) -> Zip

Iterate over (coefficient, state_word) pairs.
"""
terms(sp::StatePolynomial) = zip(sp.coeffs, sp.state_words)

"""
    degree(sp::StatePolynomial) -> Int

Compute the maximum degree of all state words in the polynomial.
Returns 0 for the zero polynomial.
"""
function degree(sp::StatePolynomial)
    isempty(sp.state_words) && return 0
    maximum(degree, sp.state_words)
end

"""
    variables(sp::StatePolynomial) -> Set

Get the set of all variable indices used in the StatePolynomial.
"""
function variables(sp::StatePolynomial{C,ST,A,T}) where {C,ST,A,T}
    result = Set{T}()
    for sw in sp.state_words
        union!(result, variables(sw))
    end
    result
end

# =============================================================================
# zero and one
# =============================================================================

"""
    Base.zero(::Type{StatePolynomial{C,ST,A,T}}) where {C,ST,A,T}

Create the zero StatePolynomial (no terms).
"""
function Base.zero(::Type{StatePolynomial{C,ST,A,T}}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    StatePolynomial(C[], StateWord{ST,A,T}[])
end

"""
    Base.zero(sp::StatePolynomial{C,ST,A,T}) where {C,ST,A,T}

Create the zero StatePolynomial for the same type as `sp`.
"""
function Base.zero(::StatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    zero(StatePolynomial{C,ST,A,T})
end

"""
    Base.iszero(sp::StatePolynomial) -> Bool

Check if a StatePolynomial is zero (has no terms).
"""
Base.iszero(sp::StatePolynomial) = isempty(sp.state_words)

"""
    Base.one(::Type{StatePolynomial{C,ST,A,T}}) where {C,ST,A,T}

Create the identity StatePolynomial (1 * identity StateWord).
"""
function Base.one(::Type{StatePolynomial{C,ST,A,T}}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    StatePolynomial([one(C)], [one(StateWord{ST,A,T})])
end

"""
    Base.one(sp::StatePolynomial{C,ST,A,T}) where {C,ST,A,T}

Create the identity StatePolynomial for the same type as `sp`.
"""
function Base.one(::StatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    one(StatePolynomial{C,ST,A,T})
end

"""
    Base.isone(sp::StatePolynomial) -> Bool

Check if a StatePolynomial is the identity.
"""
function Base.isone(sp::StatePolynomial)
    length(sp.state_words) == 1 || return false
    isone(sp.coeffs[1]) && isone(sp.state_words[1])
end

# =============================================================================
# Equality and Hashing
# =============================================================================

"""
    Base.:(==)(a::StatePolynomial, b::StatePolynomial) -> Bool

Check if two StatePolynomials are equal.
"""
function Base.:(==)(
    a::StatePolynomial{C1,ST1,A1,T1}, b::StatePolynomial{C2,ST2,A2,T2}
) where {C1,C2,ST1,ST2,A1,A2,T1,T2}
    ST1 !== ST2 && return false
    A1 !== A2 && return false
    T1 !== T2 && return false
    a.state_words != b.state_words && return false
    a.coeffs != b.coeffs && return false
    return true
end

"""
    Base.hash(sp::StatePolynomial, h::UInt) -> UInt

Hash function for StatePolynomial.
"""
function Base.hash(sp::StatePolynomial, h::UInt)
    h = hash(:StatePolynomial, h)
    for (c, sw) in zip(sp.coeffs, sp.state_words)
        h = hash(c, hash(sw, h))
    end
    return h
end

# =============================================================================
# Arithmetic - Addition
# =============================================================================

"""
    Base.:(+)(a::StatePolynomial{C1,ST,A,T}, b::StatePolynomial{C2,ST,A,T}) where {C1,C2,ST,A,T}

Add two StatePolynomials.
"""
function Base.:(+)(
    a::StatePolynomial{C1,ST,A,T}, b::StatePolynomial{C2,ST,A,T}
) where {C1<:Number,C2<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    C = promote_type(C1, C2)
    StatePolynomial(
        C[C.(a.coeffs); C.(b.coeffs)],
        [a.state_words; b.state_words]
    )
end

# =============================================================================
# Arithmetic - Subtraction
# =============================================================================

"""
    Base.:(-)(sp::StatePolynomial)

Negate a StatePolynomial.
"""
function Base.:(-)(sp::StatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    StatePolynomial(-sp.coeffs, copy(sp.state_words))
end

"""
    Base.:(-)(a::StatePolynomial, b::StatePolynomial)

Subtract two StatePolynomials.
"""
function Base.:(-)(
    a::StatePolynomial{C1,ST,A,T}, b::StatePolynomial{C2,ST,A,T}
) where {C1<:Number,C2<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    a + (-b)
end

# =============================================================================
# Arithmetic - Scalar Multiplication
# =============================================================================

"""
    Base.:(*)(c::Number, sp::StatePolynomial)

Scalar multiplication (scalar on left).
"""
function Base.:(*)(c::Number, sp::StatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    iszero(c) && return zero(StatePolynomial{promote_type(typeof(c), C),ST,A,T})
    NC = promote_type(typeof(c), C)
    StatePolynomial(NC.(c .* sp.coeffs), copy(sp.state_words))
end

"""
    Base.:(*)(sp::StatePolynomial, c::Number)

Scalar multiplication (scalar on right).
"""
Base.:(*)(sp::StatePolynomial, c::Number) = c * sp

# =============================================================================
# Arithmetic - Polynomial Multiplication
# =============================================================================

"""
    Base.:(*)(a::StatePolynomial{C1,ST,A,T}, b::StatePolynomial{C2,ST,A,T}) where {C1,C2,ST,A,T}

Multiply two StatePolynomials.
"""
function Base.:(*)(
    a::StatePolynomial{C1,ST,A,T}, b::StatePolynomial{C2,ST,A,T}
) where {C1<:Number,C2<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    # Handle zero cases
    isempty(a.state_words) && return zero(StatePolynomial{promote_type(C1, C2),ST,A,T})
    isempty(b.state_words) && return zero(StatePolynomial{promote_type(C1, C2),ST,A,T})

    C = promote_type(C1, C2)
    result_coeffs = C[]
    result_sws = StateWord{ST,A,T}[]

    for (ca, swa) in zip(a.coeffs, a.state_words)
        for (cb, swb) in zip(b.coeffs, b.state_words)
            # StateWord multiplication returns Term{StateWord, Float64}
            prod_term = swa * swb
            push!(result_coeffs, C(ca * cb * prod_term.coefficient))
            push!(result_sws, prod_term.monomial)
        end
    end

    StatePolynomial(result_coeffs, result_sws)
end

# =============================================================================
# Display
# =============================================================================

"""
    Base.show(io::IO, sp::StatePolynomial{C,ST,A,T}) where {C,ST,A,T}

Display a StatePolynomial.
"""
function Base.show(io::IO, sp::StatePolynomial{C,ST,A,T}) where {C,ST,A,T}
    if isempty(sp.state_words)
        print(io, "0")
        return
    end

    first_term = true
    for (c, sw) in zip(sp.coeffs, sp.state_words)
        if !first_term
            print(io, " + ")
        end
        print(io, c, " * ", sw)
        first_term = false
    end
end

# =============================================================================
# NCStatePolynomial
# =============================================================================

"""
    NCStatePolynomial{C<:Number, ST<:StateType, A<:AlgebraType, T<:Integer}

A polynomial in non-commutative state words with coefficients.

Represents a sum of NCStateWords with coefficients: sum_i c_i * ncsw_i

# Fields
- `coeffs::Vector{C}`: Coefficients for each NC state word
- `nc_state_words::Vector{NCStateWord{ST,A,T}}`: Sorted unique NC state words

# Invariants
Same as StatePolynomial: sorted, unique, non-zero coefficients.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1]);

julia> sw = StateWord{Arbitrary}([m1]);

julia> ncsw = NCStateWord(sw, m1);

julia> ncsp = NCStatePolynomial([1.0], [ncsw]);

julia> length(ncsp.nc_state_words)
1
```

See also: [`NCStateWord`](@ref), [`StatePolynomial`](@ref)
"""
struct NCStatePolynomial{C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    coeffs::Vector{C}
    nc_state_words::Vector{NCStateWord{ST,A,T}}

    function NCStatePolynomial(
        coeffs::Vector{C}, nc_state_words::Vector{NCStateWord{ST,A,T}}
    ) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
        # Handle empty case
        isempty(nc_state_words) && return new{C,ST,A,T}(C[], NCStateWord{ST,A,T}[])

        # Sort, combine duplicates, filter zeros
        perm = sortperm(nc_state_words)
        sorted_ncsws = nc_state_words[perm]
        sorted_coeffs = coeffs[perm]

        result_coeffs = C[]
        result_ncsws = NCStateWord{ST,A,T}[]

        current_ncsw = sorted_ncsws[1]
        current_coef = sorted_coeffs[1]

        for i in 2:length(sorted_ncsws)
            if sorted_ncsws[i] == current_ncsw
                current_coef += sorted_coeffs[i]
            else
                if !iszero(current_coef)
                    push!(result_coeffs, current_coef)
                    push!(result_ncsws, current_ncsw)
                end
                current_ncsw = sorted_ncsws[i]
                current_coef = sorted_coeffs[i]
            end
        end
        if !iszero(current_coef)
            push!(result_coeffs, current_coef)
            push!(result_ncsws, current_ncsw)
        end

        new{C,ST,A,T}(result_coeffs, result_ncsws)
    end
end

# =============================================================================
# NCStatePolynomial - Accessors
# =============================================================================

"""
    coefficients(ncsp::NCStatePolynomial) -> Vector{C}

Get the coefficients of the NCStatePolynomial.
"""
coefficients(ncsp::NCStatePolynomial) = ncsp.coeffs

"""
    monomials(ncsp::NCStatePolynomial) -> Vector{NCStateWord}

Get the NC state words of the NCStatePolynomial.
"""
monomials(ncsp::NCStatePolynomial) = ncsp.nc_state_words

"""
    terms(ncsp::NCStatePolynomial) -> Zip

Iterate over (coefficient, nc_state_word) pairs.
"""
terms(ncsp::NCStatePolynomial) = zip(ncsp.coeffs, ncsp.nc_state_words)

"""
    degree(ncsp::NCStatePolynomial) -> Int

Compute the maximum degree of all NC state words in the polynomial.
"""
function degree(ncsp::NCStatePolynomial)
    isempty(ncsp.nc_state_words) && return 0
    maximum(degree, ncsp.nc_state_words)
end

"""
    variables(ncsp::NCStatePolynomial) -> Set

Get the set of all variable indices used in the NCStatePolynomial.
"""
function variables(ncsp::NCStatePolynomial{C,ST,A,T}) where {C,ST,A,T}
    result = Set{T}()
    for ncsw in ncsp.nc_state_words
        union!(result, variables(ncsw))
    end
    result
end

# =============================================================================
# NCStatePolynomial - zero and one
# =============================================================================

"""
    Base.zero(::Type{NCStatePolynomial{C,ST,A,T}}) where {C,ST,A,T}

Create the zero NCStatePolynomial.
"""
function Base.zero(::Type{NCStatePolynomial{C,ST,A,T}}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStatePolynomial(C[], NCStateWord{ST,A,T}[])
end

"""
    Base.zero(ncsp::NCStatePolynomial{C,ST,A,T}) where {C,ST,A,T}

Create the zero NCStatePolynomial for the same type.
"""
function Base.zero(::NCStatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    zero(NCStatePolynomial{C,ST,A,T})
end

"""
    Base.iszero(ncsp::NCStatePolynomial) -> Bool

Check if an NCStatePolynomial is zero.
"""
Base.iszero(ncsp::NCStatePolynomial) = isempty(ncsp.nc_state_words)

"""
    Base.one(::Type{NCStatePolynomial{C,ST,A,T}}) where {C,ST,A,T}

Create the identity NCStatePolynomial.
"""
function Base.one(::Type{NCStatePolynomial{C,ST,A,T}}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStatePolynomial([one(C)], [one(NCStateWord{ST,A,T})])
end

"""
    Base.one(ncsp::NCStatePolynomial{C,ST,A,T}) where {C,ST,A,T}

Create the identity NCStatePolynomial for the same type.
"""
function Base.one(::NCStatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    one(NCStatePolynomial{C,ST,A,T})
end

"""
    Base.isone(ncsp::NCStatePolynomial) -> Bool

Check if an NCStatePolynomial is the identity.
"""
function Base.isone(ncsp::NCStatePolynomial)
    length(ncsp.nc_state_words) == 1 || return false
    isone(ncsp.coeffs[1]) && isone(ncsp.nc_state_words[1])
end

# =============================================================================
# NCStatePolynomial - Equality and Hashing
# =============================================================================

"""
    Base.:(==)(a::NCStatePolynomial, b::NCStatePolynomial) -> Bool

Check if two NCStatePolynomials are equal.
"""
function Base.:(==)(
    a::NCStatePolynomial{C1,ST1,A1,T1}, b::NCStatePolynomial{C2,ST2,A2,T2}
) where {C1,C2,ST1,ST2,A1,A2,T1,T2}
    ST1 !== ST2 && return false
    A1 !== A2 && return false
    T1 !== T2 && return false
    a.nc_state_words != b.nc_state_words && return false
    a.coeffs != b.coeffs && return false
    return true
end

"""
    Base.hash(ncsp::NCStatePolynomial, h::UInt) -> UInt

Hash function for NCStatePolynomial.
"""
function Base.hash(ncsp::NCStatePolynomial, h::UInt)
    h = hash(:NCStatePolynomial, h)
    for (c, ncsw) in zip(ncsp.coeffs, ncsp.nc_state_words)
        h = hash(c, hash(ncsw, h))
    end
    return h
end

# =============================================================================
# NCStatePolynomial - Arithmetic
# =============================================================================

"""
    Base.:(+)(a::NCStatePolynomial, b::NCStatePolynomial)

Add two NCStatePolynomials.
"""
function Base.:(+)(
    a::NCStatePolynomial{C1,ST,A,T}, b::NCStatePolynomial{C2,ST,A,T}
) where {C1<:Number,C2<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    C = promote_type(C1, C2)
    NCStatePolynomial(
        C[C.(a.coeffs); C.(b.coeffs)],
        [a.nc_state_words; b.nc_state_words]
    )
end

"""
    Base.:(-)(ncsp::NCStatePolynomial)

Negate an NCStatePolynomial.
"""
function Base.:(-)(ncsp::NCStatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    NCStatePolynomial(-ncsp.coeffs, copy(ncsp.nc_state_words))
end

"""
    Base.:(-)(a::NCStatePolynomial, b::NCStatePolynomial)

Subtract two NCStatePolynomials.
"""
function Base.:(-)(
    a::NCStatePolynomial{C1,ST,A,T}, b::NCStatePolynomial{C2,ST,A,T}
) where {C1<:Number,C2<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    a + (-b)
end

"""
    Base.:(*)(c::Number, ncsp::NCStatePolynomial)

Scalar multiplication (scalar on left).
"""
function Base.:(*)(c::Number, ncsp::NCStatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    iszero(c) && return zero(NCStatePolynomial{promote_type(typeof(c), C),ST,A,T})
    NC = promote_type(typeof(c), C)
    NCStatePolynomial(NC.(c .* ncsp.coeffs), copy(ncsp.nc_state_words))
end

"""
    Base.:(*)(ncsp::NCStatePolynomial, c::Number)

Scalar multiplication (scalar on right).
"""
Base.:(*)(ncsp::NCStatePolynomial, c::Number) = c * ncsp

# =============================================================================
# NCStatePolynomial - expval
# =============================================================================

"""
    expval(ncsp::NCStatePolynomial) -> StatePolynomial

Compute the expectation value of an NCStatePolynomial.

Converts each NCStateWord to a StateWord using `expval` and combines with coefficients.
"""
function expval(ncsp::NCStatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    isempty(ncsp.nc_state_words) && return zero(StatePolynomial{C,ST,A,T})

    result_sws = [expval(ncsw) for ncsw in ncsp.nc_state_words]
    StatePolynomial(copy(ncsp.coeffs), result_sws)
end

# =============================================================================
# NCStatePolynomial - Display
# =============================================================================

"""
    Base.show(io::IO, ncsp::NCStatePolynomial{C,ST,A,T}) where {C,ST,A,T}

Display an NCStatePolynomial.
"""
function Base.show(io::IO, ncsp::NCStatePolynomial{C,ST,A,T}) where {C,ST,A,T}
    if isempty(ncsp.nc_state_words)
        print(io, "0")
        return
    end

    first_term = true
    for (c, ncsw) in zip(ncsp.coeffs, ncsp.nc_state_words)
        if !first_term
            print(io, " + ")
        end
        print(io, c, " * ", ncsw)
        first_term = false
    end
end
