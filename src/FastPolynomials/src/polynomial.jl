abstract type AbstractPolynomial{T} end
"""
    Polynomial{T}

A polynomial represented as a sum of monomials with coefficients.
Maintains sorted, unique monomials with non-zero coefficients.

# Fields
- `coeffs::Vector{T}`: Non-zero coefficients
- `monos::Vector{Monomial}`: Sorted, unique monomials

# Invariants
- `monos` are always sorted and non-repeating
- `coeffs` are always non-zero
- Length of `coeffs` equals length of `monos`
"""
struct Polynomial{T} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    # perhaps make it into an ordered set?
    monos::Vector{Monomial}

    function Polynomial(a::Vector{T}, x::Vector{Monomial}) where {T}
        length(a) == length(x) ||
            throw(ArgumentError("There should be as many coefficient than monomials"))
        sorted_x_idx = sortperm(x)
        sorted_x = x[sorted_x_idx]
        unique_monos = sort!(unique(x))
        unique_coeffs = map(unique_monos) do mono
            mapreduce(
                idx -> a[sorted_x_idx[idx]],
                +,
                searchsortedfirst(sorted_x, mono):searchsortedlast(sorted_x, mono);
                init=zero(T),
            )
        end
        nz_idx = findall(!iszero, unique_coeffs)
        return new{T}(unique_coeffs[nz_idx], unique_monos[nz_idx])
    end
end

Polynomial(a) = Polynomial(Float64, a)
Polynomial(::Type{T}, a::Variable) where {T} = Polynomial([one(T)], [monomial(a)])
Polynomial(::Type{T}, a::Monomial) where {T} = Polynomial([one(T)], [a])
function Polynomial(::Type{T1}, a::T2) where {T1<:Number,T2<:Number}
    return Polynomial([promote_type(T1, T2)(a)], [one(Monomial)])
end
Polynomial(a::Polynomial) = a

function Base.convert(::Type{Polynomial{T}}, a::Variable) where {T}
    return Polynomial([one(T)], [Monomial([a], [1])])
end
Base.convert(::Type{Polynomial{T}}, a::Monomial) where {T} = Polynomial([one(T)], [a])

"""
    variables(p::Polynomial)

Extracts all variables appearing in a polynomial.

# Arguments
- `p::Polynomial`: The polynomial

# Returns
- Set of all variables in the polynomial
"""
variables(p::Polynomial) = union(variables.(p.monos)...)

function Base.show(io::IO, obj::Polynomial)
    return print_object(io, obj; multiline=true)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::Polynomial)
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

function print_object(io::IO, obj::Polynomial{T}; multiline::Bool) where {T}
    if multiline
        return print(io, join(
            map(zip(obj.coeffs, obj.monos)) do (coef, mono)
                T <: Real ? ("$(coef) * $(string(mono))") : ("($(coef)) * $(string(mono))")
            end,
            " + ",
        ))
    else
        return Base.show_default(io, obj)
    end
end

function Base.hash(p::Polynomial, u::UInt)
    return hash(p.coeffs, hash(p.monos, u))
end

maxdegree(p::Polynomial) = maximum(degree.(p.monos))

Base.:(^)(p::Polynomial, n::Int) = Base.power_by_squaring(p, n)

Base.one(::Type{Polynomial{T}}) where {T} = Polynomial([one(T)], [one(Monomial)])
Base.one(::Polynomial{T}) where {T} = Polynomial([one(T)], [one(Monomial)])
Base.zero(::Polynomial{T}) where {T} = Polynomial([zero(T)], [one(Monomial)])

Base.real(p::Polynomial) = Polynomial(real.(p.coeffs), p.monos)

coefficients(p::Polynomial) = p.coeffs
monomials(p::Polynomial) = p.monos
terms(p::Polynomial) = zip(p.coeffs, p.monos)
expval(p::Polynomial) = p

"""
    support(poly::Polynomial{T}) where {T}

Computes the support of a polynomial after canonicalization.

# Arguments
- `poly::Polynomial{T}`: The polynomial

# Returns
- `Vector{Monomial}`: Unique monomials from the support of a polynomial
"""
function support(poly::Polynomial{T}) where {T}
    return unique!(poly.monos)
end
