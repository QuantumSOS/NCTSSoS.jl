"User polynomial with an owning registry. Arithmetic allocates; PolynomialBuffer kernels remain caller-owned."
struct Polynomial{W,C,R}
    registry::R
    data::PolynomialBuffer{W,C}
end

_locations(w::PauliWord) = ((div(l-1,3)+1,mod(l-1,3)) for l in legacy_letters(w))
_locations(w::FermiWord) = ((abs(l),sign(l)) for l in legacy_letters(w))
_locations(w::SiteWord{A,U,B,S}) where {A,U,B,S} =
    ((_site(typeof(w),_letter(w,i)),Int(_letter(w,i)>>S)) for i in 1:degree(w))
"Display a word in the names of its registry, without changing kernel word storage."
function Base.show(io::IO,w::AbstractWord,r::Registry)
    isone(w) && return print(io,"1")
    for (i,location) in enumerate(_locations(w))
        i>1 && print(io," ")
        name = findfirst(==(location),r.names)
        if isnothing(name)
            print(io,"operator(site=",location[1],", label=",location[2],")")
        else
            print(io,name)
        end
    end
end
_word_string(r::Registry,w::AbstractWord) = sprint(show,w,r)
function Base.show(io::IO,p::Polynomial)
    isempty(p.data) && return print(io,"0")
    for (i,(c,w)) in enumerate(terms(p))
        i>1 && print(io," + ")
        if isone(w)
            print(io,isreal(c) ? real(c) : c)
        else
            c==1 || print(io,"(",isreal(c) ? real(c) : c,") ")
            show(io,w,p.registry)
        end
    end
end
Base.show(io::IO,::MIME"text/plain",p::Polynomial) = show(io,p)

terms(p::Polynomial) = ((p.data.coefficients[i],p.data.words[i]) for i in 1:length(p.data))
degree(p::Polynomial) = maximum(degree,p.data.words[1:length(p.data)];init=0)
Base.iszero(p::Polynomial) = isempty(p.data)
function _polynomial(r,data::PolynomialBuffer)
    all(isfinite,data.coefficients[1:length(data)]) || throw(ArgumentError("Polynomial coefficients must be finite; replace NaN or Inf with a finite coefficient."))
    return Polynomial(r,canonicalize!(data))
end
function _converted(p::Polynomial{W},::Type{C}) where {W,C}
    b = PolynomialBuffer{W,C}(length(p.data))
    for (c,w) in terms(p); push_term!(b,c,w); end
    return canonicalize!(b)
end
function _constant(p::Polynomial{W,C},c::Number) where {W,C}
    b = PolynomialBuffer{W,promote_type(C,typeof(c))}(1)
    push_term!(b,c,one(W))
    return _polynomial(p.registry,b)
end
Base.zero(p::Polynomial) = _constant(p,0)
Base.one(p::Polynomial) = _constant(p,1)
function _same_registry(a::Polynomial,b::Polynomial)
    a.registry === b.registry || throw(ArgumentError("Polynomials use different registries; construct all variables of this problem from the same registry."))
end
function _add(a::Polynomial{W,C},b::Polynomial{W,D},sign) where {W,C,D}
    _same_registry(a,b)
    T = promote_type(C,D)
    out = PolynomialBuffer{W,T}(length(a.data)+length(b.data))
    add!(out,_converted(a,T),_converted(b,T),one(T),sign)
    return _polynomial(a.registry,out)
end
Base.:+(a::Polynomial,b::Polynomial) = _add(a,b,1)
Base.:-(a::Polynomial,b::Polynomial) = _add(a,b,-1)
function Base.:*(a::Polynomial{W,C},b::Polynomial{W,D}) where {W,C,D}
    _same_registry(a,b)
    T = promote_type(C,D)
    out,scratch = PolynomialBuffer{W,T}(),PolynomialBuffer{W,T}()
    mul!(out,_converted(a,T),_converted(b,T),scratch)
    return _polynomial(a.registry,out)
end
# Different word capacities/algebras must produce a registry error, not a dispatch error.
_add(a::Polynomial,b::Polynomial,sign) = (_same_registry(a,b); throw(ArgumentError("Polynomial word types differ; use variables constructed together.")))
Base.:*(a::Polynomial,b::Polynomial) = (_same_registry(a,b); throw(ArgumentError("Polynomial word types differ; use variables constructed together.")))
Base.:+(a::Polynomial,b::Number) = a+_constant(a,b)
Base.:+(a::Number,b::Polynomial) = _constant(b,a)+b
Base.:-(a::Polynomial,b::Number) = a-_constant(a,b)
Base.:-(a::Number,b::Polynomial) = _constant(b,a)-b
Base.:*(a::Polynomial,b::Number) = a*_constant(a,b)
Base.:*(a::Number,b::Polynomial) = _constant(b,a)*b
Base.:-(p::Polynomial) = -1*p
function Base.:/(p::Polynomial,c::Number)
    isfinite(c) && !iszero(c) || throw(ArgumentError("Polynomial divisor must be finite and nonzero; replace the divisor."))
    return p*inv(c)
end
function Base.adjoint(p::Polynomial{W,C}) where {W,C}
    out = PolynomialBuffer{W,C}(length(p.data))
    star!(out,p.data)
    return _polynomial(p.registry,out)
end
function Base.:^(p::Polynomial,n::Integer)
    n >= 0 || throw(ArgumentError("Polynomial exponent $n is negative; use a nonnegative integer exponent."))
    result = one(p)
    while n > 0
        isodd(n) && (result=result*p)
        n >>= 1
        n>0 && (p=p*p)
    end
    return result
end
Base.inv(::Polynomial) = throw(ArgumentError("Polynomial inverse is not supported; use nonnegative integer powers instead of a negative exponent."))
function Base.:(==)(a::Polynomial,b::Polynomial)
    a.registry === b.registry || return false
    return length(a.data)==length(b.data) && all(a.data.words[i]==b.data.words[i] && a.data.coefficients[i]==b.data.coefficients[i] for i in 1:length(a.data))
end
LinearAlgebra.ishermitian(p::Polynomial) = p==adjoint(p)

_model_word_type(r::Registry) = word_type(r,0)
function _model_word_type(r::Registry{<:SiteAlgebra})
    bits = max(1,ndigits(r.nsites;base=2))+ndigits(r.algebra.nletters;base=2)
    return word_type(r,128÷bits)
end
"Named variables on an existing registry; site words use the full supported modeling capacity."
function variables(r::Registry)
    W = _model_word_type(r)
    return Dict(name=>_variable(r,name,W) for name in keys(r.names))
end
function _variable(r,name,W)
    b = PolynomialBuffer{W,ComplexF64}(1)
    push_term!(b,1,word(r,name,W))
    return _polynomial(r,b)
end
function _site_count(sites)
    ids = sites isa Integer ? collect(1:sites) : collect(sites)
    !isempty(ids) && ids==collect(1:length(ids)) || throw(ArgumentError("Variable sites must be 1:n with n > 0; use Registry and variables for custom site names or locations."))
    return length(ids)
end
function pauli_variables(sites)
    n = _site_count(sites)
    r = Registry(Pauli(),n,[(Symbol(axis,i),(i,a-1)) for (a,axis) in enumerate((:X,:Y,:Z)) for i in 1:n])
    v = variables(r)
    return r,Tuple([v[Symbol(axis,i)] for i in 1:n] for axis in (:X,:Y,:Z))
end
function fermionic_variables(sites)
    n = _site_count(sites)
    r = Registry(Fermionic(),n,[(Symbol(name,i),(i,label)) for (name,label) in ((:a,1),(:ad,-1)) for i in 1:n])
    v = variables(r)
    return r,Tuple([v[Symbol(name,i)] for i in 1:n] for name in (:a,:ad))
end
"Self-adjoint operators at each commuting site; output[site][letter]."
function site_variables(rule::SiteRule,sites,nletters::Integer)
    n = _site_count(sites)
    r = Registry(SiteAlgebra(rule,nletters),n,[(Symbol(:s,i,:_,j),(i,j)) for i in 1:n for j in 1:nletters])
    v = variables(r)
    return r,[[v[Symbol(:s,i,:_,j)] for j in 1:nletters] for i in 1:n]
end
