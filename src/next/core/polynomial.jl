"Caller-owned storage. Reserve ahead for allocation-free kernels; otherwise storage grows geometrically."
mutable struct PolynomialBuffer{W<:AbstractWord,C<:Number}
    words::Vector{W}
    coefficients::Vector{C}
    n::Int
    canonical::Bool
    function PolynomialBuffer{W,C}(capacity::Integer=0) where {W<:AbstractWord,C<:Number}
        capacity >= 0 || throw(ArgumentError("Polynomial capacity must be nonnegative."))
        new{W,C}(Vector{W}(undef, capacity), Vector{C}(undef, capacity), 0, true)
    end
end
Base.length(p::PolynomialBuffer) = p.n
Base.isempty(p::PolynomialBuffer) = p.n == 0
Base.empty!(p::PolynomialBuffer) = (p.n = 0; p.canonical = true; p)

"Explicitly allocate capacity outside arithmetic kernels; never shrinks live data."
function reserve!(p::PolynomialBuffer, capacity::Integer)
    capacity >= p.n || throw(ArgumentError("Capacity cannot be smaller than the live polynomial."))
    resize!(p.words, capacity)
    resize!(p.coefficients, capacity)
    return p
end

function _grow!(p::PolynomialBuffer, needed::Int)
    needed <= length(p.words) || reserve!(p, max(needed, max(4, 2length(p.words))))
    return p
end

"Append a term without sorting or coalescing; call canonicalize! before arithmetic."
function push_term!(p::PolynomialBuffer{W,C}, c::Number, w::W) where {W,C}
    cc = convert(C, c)
    p.canonical = false
    iszero(cc) && return p
    _grow!(p, p.n + 1)
    p.n += 1
    @inbounds p.words[p.n], p.coefficients[p.n] = w, cc
    return p
end

@inline function _swap_terms!(p::PolynomialBuffer, i, j)
    @inbounds begin
        p.words[i], p.words[j] = p.words[j], p.words[i]
        p.coefficients[i], p.coefficients[j] = p.coefficients[j], p.coefficients[i]
    end
end
function _sift_terms!(p::PolynomialBuffer, root::Int, n::Int)
    while 2root <= n
        child = 2root
        @inbounds child < n && isless(p.words[child], p.words[child + 1]) && (child += 1)
        @inbounds isless(p.words[root], p.words[child]) || break
        _swap_terms!(p, root, child)
        root = child
    end
    return p
end

"In-place paired heapsort, exact duplicate reduction, and zero removal; no scratch allocation."
function canonicalize!(p::PolynomialBuffer)
    p.canonical && return p
    for root in (p.n ÷ 2):-1:1
        _sift_terms!(p, root, p.n)
    end
    for stop in p.n:-1:2
        _swap_terms!(p, 1, stop)
        _sift_terms!(p, 1, stop - 1)
    end
    n, i = 0, 1
    @inbounds while i <= p.n
        w, c = p.words[i], p.coefficients[i]
        i += 1
        while i <= p.n && p.words[i] == w
            c += p.coefficients[i]
            i += 1
        end
        if !iszero(c)
            n += 1
            p.words[n], p.coefficients[n] = w, c
        end
    end
    p.n = n
    p.canonical = true
    return p
end

"Copy a polynomial into caller storage, growing it when needed."
function Base.copyto!(out::PolynomialBuffer{W,C}, p::PolynomialBuffer{W,C}) where {W,C}
    _grow!(out, p.n)
    copyto!(out.words, 1, p.words, 1, p.n)
    copyto!(out.coefficients, 1, p.coefficients, 1, p.n)
    out.n = p.n
    out.canonical = p.canonical
    return out
end

"out = α*a + β*b for canonical inputs. Destination must not alias either input."
function add!(out::PolynomialBuffer{W,C}, a::PolynomialBuffer{W,C}, b::PolynomialBuffer{W,C}, α::Number=one(C), β::Number=one(C)) where {W,C}
    out !== a && out !== b || throw(ArgumentError("add! needs a destination distinct from both inputs."))
    a.canonical && b.canonical || throw(ArgumentError("add! requires canonical inputs; call canonicalize! on both inputs first."))
    empty!(out)
    i = j = 1
    @inbounds while i <= a.n || j <= b.n
        if j > b.n || (i <= a.n && isless(a.words[i], b.words[j]))
            push_term!(out, α * a.coefficients[i], a.words[i]); i += 1
        elseif i > a.n || isless(b.words[j], a.words[i])
            push_term!(out, β * b.coefficients[j], b.words[j]); j += 1
        else
            push_term!(out, α * a.coefficients[i] + β * b.coefficients[j], a.words[i]); i += 1; j += 1
        end
    end
    out.canonical = true
    return out
end

"out = a*b. Scratch holds unmerged products; out may alias an input, scratch may not."
function mul!(out::PolynomialBuffer{W,C}, a::PolynomialBuffer{W,C}, b::PolynomialBuffer{W,C}, scratch::PolynomialBuffer{W,C}) where {W,C}
    scratch !== a && scratch !== b && scratch !== out || throw(ArgumentError("Product scratch must be distinct from the inputs and destination."))
    empty!(scratch)
    @inbounds for j in 1:b.n, i in 1:a.n
        mul!(scratch, a.words[i], b.words[j], a.coefficients[i] * b.coefficients[j])
    end
    canonicalize!(scratch)
    return copyto!(out, scratch)
end

"out = adjoint(p); destination may be the input."
function star!(out::PolynomialBuffer{W,C}, p::PolynomialBuffer{W,C}) where {W,C}
    _grow!(out, p.n)
    out.canonical = false
    @inbounds for i in 1:p.n
        phase, w = star(p.words[i])
        out.words[i], out.coefficients[i] = w, scale(phase, conj(p.coefficients[i]))
    end
    out.n = p.n
    return canonicalize!(out)
end
