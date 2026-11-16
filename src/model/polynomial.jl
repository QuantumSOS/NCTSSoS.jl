"User polynomial with an owning registry. Arithmetic allocates; PolynomialBuffer kernels remain caller-owned."
struct Polynomial{W,C,R}
    registry::R
    data::PolynomialBuffer{W,C}
end

_locations(w::PauliWord) = ((only(sites(l)),all(iszero,l.x) ? 2 : all(iszero,l.z) ? 0 : 1) for l in word_letters(w))
_locations(w::FermiWord) = ((only(sites(l)),all(iszero,l.c) ? 1 : -1) for l in word_letters(w))
_locations(w::BosonWord) = ((only(sites(l)),all(iszero,l.c) ? 1 : -1) for l in word_letters(w))
_locations(w::SiteWord{A,U,B,S}) where {A,U,B,S} =
    ((_site(typeof(w),_letter(w,i)),Int(_letter(w,i)>>S)) for i in 1:degree(w))
"Display a word in the names of its registry, without changing kernel word storage."
function Base.show(io::IO,w::AbstractWord,r::Union{Registry,Nothing})
    isone(w) && return print(io,"1")
    for (i,location) in enumerate(_locations(w))
        i>1 && print(io," ")
        name = isnothing(r) ? nothing : findfirst(==(location),r.names)
        if isnothing(name)
            print(io,"operator(site=",location[1],", label=",location[2],")")
        else
            print(io,name)
        end
    end
end
function Base.show(io::IO,w::Union{PauliWord,FermiWord,BosonWord,SiteWord,FunctionalWord,MixedWord})
    isone(w) && return print(io,"1")
    registry=get(io,:registry,nothing)
    isnothing(registry) && print(io,nameof(typeof(w)),"(")
    show(io,w,registry)
    isnothing(registry) && print(io,")")
end
_word_string(r::Registry,w::AbstractWord) = sprint(show,w,r)
function Base.show(io::IO,p::Polynomial)
    isempty(p.data) && return print(io,"0")
    for (i,(c,w)) in enumerate(terms(p))
        i>8 && get(io,:limit,true) && return print(io," + … (",length(p.data)," terms)")
        i>1 && print(io," + ")
        shown=isreal(c) ? real(c) : complex(iszero(real(c)) ? zero(real(c)) : real(c),imag(c))
        if isone(w)
            print(io,shown)
        else
            c==1 || print(io,"(",shown,") ")
            show(io,w,p.registry)
        end
    end
end
Base.show(io::IO,::MIME"text/plain",p::Polynomial) = show(io,p)

"""
Borrowed, sized iterator of (coefficient, word) pairs. Read it without mutating the polynomial.
A display shows at most eight terms, including print and string; iterate NCTSSoS.terms(p)
for all terms, or show with IOContext(io, :limit => false).
To print a returned word with its operator names, use
`sprint(show, w; context=:registry => registry)` with the polynomial's registry.
"""
terms(p::Polynomial) = ((p.data.coefficients[i],p.data.words[i]) for i in 1:length(p.data))
degree(p::Polynomial) = maximum((degree(w) for (_,w) in terms(p));init=0)
Base.iszero(p::Polynomial) = isempty(p.data)
function _polynomial(r,data::PolynomialBuffer)
    all(isfinite(data.coefficients[i]) for i in 1:length(data)) || throw(ArgumentError("Polynomial coefficients must be finite; replace NaN or Inf with a finite coefficient."))
    return Polynomial(r,canonicalize!(data))
end
# Arithmetic kernels borrow their inputs without mutating them.
_converted(p::Polynomial{W,C},::Type{C}) where {W,C} = p.data
function _converted(p::Polynomial{W},::Type{C}) where {W,C}
    b = PolynomialBuffer{W,C}(length(p.data))
    for (c,w) in terms(p); push_term!(b,c,w); end
    b.canonical = true
    return b
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
# Match value equality across coefficient types and signed zeros. As with any
# mutable key, do not change the backing terms while a polynomial is in a Set/Dict.
function Base.hash(p::Polynomial,h::UInt)
    h=hash(objectid(p.registry),hash(:Polynomial,h))
    for (c,w) in terms(p)
        re,im=real(c),imag(c)
        h=hash(iszero(im) ? 0 : im,hash(iszero(re) ? 0 : re,hash(w,h)))
    end
    return h
end
LinearAlgebra.ishermitian(p::Polynomial) = p==adjoint(p)

"""
    NCTSSoS.polynomial(registry, summands)

Construct a canonical polynomial from an iterator of polynomials belonging to
registry, without mutating any summand. An empty iterator gives a registry-owned
zero. Coefficient types are promoted as summands arrive; duplicate words are
added in encounter order, with no tolerance or Hermitian projection.
Equality is by value (`==`); signed zeros are not preserved.

For fixed-width words and coefficients, accumulation uses expected O(N) work
for N input terms and O(K) storage for K distinct words. Fermionic words use
O(L*K) radix sorting for L key bytes; other domains use O(K log K) heap sorting.
Repeated addition copies the growing prefix and can allocate quadratically:
for a sum of many terms use `NCTSSoS.polynomial(registry, summands)`.
"""
function polynomial(r::Registry,summands)
    item=iterate(summands)
    W=_model_word_type(r)
    isnothing(item) && return Polynomial(r,PolynomialBuffer{W,ComplexF64}())
    p=first(item)
    p isa Polynomial || throw(ArgumentError("Summands must be polynomials; construct each term from variables in registry."))
    return _sum_polynomials(r,summands,item,Dict{W,eltype(p.data.coefficients)}())
end
function _sum_polynomials(r,summands,item,coefficients::Dict{W,C}) where {W,C}
    while !isnothing(item)
        p,state=item
        p isa Polynomial || throw(ArgumentError("Summands must be polynomials; construct each term from variables in registry."))
        p.registry === r || throw(ArgumentError("Summands use different registries; construct all variables from the supplied registry."))
        eltype(p.data.words) === W || throw(ArgumentError("Polynomial word types differ; use variables constructed together."))
        T=promote_type(C,eltype(p.data.coefficients))
        T===C || return _sum_polynomials(r,summands,item,Dict{W,T}(coefficients))
        for (c,w) in terms(p)
            coefficients[w]=get(coefficients,w,zero(C))+c
        end
        item=iterate(summands,state)
    end
    out=PolynomialBuffer{W,C}(length(coefficients))
    for (w,c) in coefficients; push_term!(out,c,w); end
    return _polynomial(r,_sort_unique!(out))
end

"Sum a polynomial array without repeated prefix copies. Empty arrays require a registry-owned polynomial init."
function Base.sum(ps::AbstractArray{<:Polynomial};init=nothing,dims=:)
    if !(dims isa Colon)
        return isnothing(init) ? mapreduce(identity,+,ps;dims) : mapreduce(identity,+,ps;init,dims)
    end
    if isempty(ps)
        init isa Polynomial || throw(ArgumentError("An empty polynomial sum needs a registry; use NCTSSoS.polynomial(registry, ()) or sum(array; init=zero(variable))."))
        return polynomial(init.registry,(init,))
    end
    r=first(ps).registry
    isnothing(init) && return polynomial(r,ps)
    initial=init isa Number ? _constant(first(ps),init) : init
    return polynomial(r,Iterators.flatten(((initial,),ps)))
end

_model_word_type(r::Registry) = word_type(r,0)
_model_word_type(r::Registry{Bosonic}) = word_type(r,255)
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
"""
    registry, (x, y, z) = pauli_variables(sites)

Return a registry and three vectors of Pauli operator polynomials.
Pass a positive site count `n` or the range `1:n` as `sites`.
Read the operators at site `i` as `x[i]`, `y[i]` and `z[i]`.
They satisfy `x[i]*y[i] = im*z[i]` and cyclic permutations, with each square equal to 1.
Operators at different sites commute.
The registry names are `X1`, `Y1`, `Z1`, and so on.

Example: `using NCTSSoS; registry, (x, y, z) = pauli_variables(2)`.
"""
function pauli_variables(sites)
    n = _site_count(sites)
    r = Registry(Pauli(),n,[(Symbol(axis,i),(i,a-1)) for (a,axis) in enumerate((:X,:Y,:Z)) for i in 1:n])
    v = variables(r)
    return r,Tuple([v[Symbol(axis,i)] for i in 1:n] for axis in (:X,:Y,:Z))
end
"""
    registry, (a, ad) = fermionic_variables(sites)

Return a registry and vectors of annihilation operators `a` and creation operators `ad`.
Pass a positive mode count `n` or the range `1:n` as `sites`.
They satisfy `a[i]*a[j] + a[j]*a[i] = 0` and
`a[i]*ad[j] + ad[j]*a[i] = δ[i,j]`, where `ad[i] = a[i]'`.
The symbol `δ[i,j]` is 1 for equal indices and 0 otherwise.
The registry names are `a1`, `ad1`, and so on.

Example: `using NCTSSoS; registry, (a, ad) = fermionic_variables(2)`.
"""
function fermionic_variables(sites)
    n = _site_count(sites)
    r = Registry(Fermionic(),n,[(Symbol(name,i),(i,label)) for (name,label) in ((:a,1),(:ad,-1)) for i in 1:n])
    v = variables(r)
    return r,Tuple([v[Symbol(name,i)] for i in 1:n] for name in (:a,:ad))
end
"Bosonic annihilation and creation variables; modeling exponents have capacity 255 per mode."
function bosonic_variables(sites)
    n=_site_count(sites)
    r=Registry(Bosonic(),n,[(Symbol(name,i),(i,label)) for (name,label) in ((:b,1),(:bd,-1)) for i in 1:n])
    v=variables(r)
    return r,Tuple([v[Symbol(name,i)] for i in 1:n] for name in (:b,:bd))
end
"Self-adjoint operators at each commuting site; output[site][letter]. Select a real functional once with functional=NCTSSoS.RealState() or NCTSSoS.RealTrace()."
function site_variables(rule::SiteRule,sites,nletters::Integer;functional=nothing)
    n = _site_count(sites)
    algebra=SiteAlgebra(rule,nletters)
    isnothing(functional) || (algebra=FunctionalAlgebra(functional,algebra))
    r = Registry(algebra,n,[(Symbol(:s,i,:_,j),(i,j)) for i in 1:n for j in 1:nletters])
    v = variables(r)
    return r,[[v[Symbol(:s,i,:_,j)] for j in 1:nletters] for i in 1:n]
end
