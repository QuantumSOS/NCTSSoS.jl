"""
    NCTSSoS.RealState()

Normalized positive functionals on the real algebra of self-adjoint generators.
Symbols of a word and its adjoint are equal; symbols are not cyclic.

For real-coefficient problems whose scalar symbols all have self-adjoint
arguments, this domain has the same optimum as complex representations:
replace each complex matrix Z by [real(Z) -imag(Z); imag(Z) real(Z)]. This
preserves products, adjoints, positivity and operator constraints. Realifying
a density matrix and dividing it by two preserves the real expectation of
every self-adjoint argument, hence also products of these expectations.
Covariance Bell objectives with commuting parties satisfy this condition.
No equivalence is claimed for general non-self-adjoint symbol arguments or
constraints with complex coefficients. This statement concerns true optima,
not equality of finite-order relaxation bounds.
"""
struct RealState end

"""
    NCTSSoS.RealTrace()

Normalized real positive traces, with cyclicity and equality under adjoint.
Matrix dimension is not fixed. The doubling argument in `RealState` also
preserves normalized traces of self-adjoint arguments, so real-coefficient
problems using only such symbols have the same true optimum as complex
representations. More generally it preserves any trace known to be real,
including the trace of a product of two self-adjoint matrices. Thus the
single-group trace covariance Bell objective also has this guarantee.
No such claim is made for arbitrary non-self-adjoint trace arguments.
"""
struct RealTrace end

struct FunctionalAlgebra{F,A<:SiteRule}
    base::SiteAlgebra{A}
end
FunctionalAlgebra(f::Union{RealState,RealTrace},a::SiteAlgebra{A}) where {A} =
    FunctionalAlgebra{typeof(f),A}(a)
FunctionalAlgebra(f,a) = throw(ArgumentError("RealState and RealTrace require Free, Projector or Unipotent site variables; Pauli, fermionic and bosonic registries are unsupported."))

# Binary header width, including one bit for a zero-capacity layout. Unlike
# generic ndigits, this folds with the type-fixed capacity in word kernels.
@inline _fheader(d) = max(1,8sizeof(d)-leading_zeros(d))
@inline function _functional_capacity(b,k)
    d=64k÷(b+1)
    while d*(b+1)+2*_fheader(d)>64k
        d-=1
    end
    d
end

"Tokens, scalar endpoints and two lengths; no registry or heap-owned fields."
struct FunctionalWord{F,A<:SiteRule,B,S,K} <: AbstractWord
    bits::NTuple{K,UInt64}
    function FunctionalWord{F,A,B,S,K}(bits::NTuple{K,UInt64}) where {F,A,B,S,K}
        W=FunctionalWord{F,A,B,S,K}
        1<=K<=8 && 0<S<B<=64 && F in (RealState,RealTrace) || throw(ArgumentError("Invalid functional word layout."))
        w=new{F,A,B,S,K}(bits)
        _check_functional_word(w)
        w
    end
    FunctionalWord{F,A,B,S,K}(bits::NTuple{K,UInt64},::Val{:normal}) where {F,A,B,S,K} = new{F,A,B,S,K}(bits)
end

@inline _fcap(::Type{<:FunctionalWord{F,A,B,S,K}}) where {F,A,B,S,K} = _functional_capacity(B,K)
@inline _fmask(b) = (UInt64(1)<<b)-UInt64(1)
@inline function _fget(bits,offset,b)
    i,s=divrem(offset,64)
    v=bits[i+1]>>s
    s+b>64 && (v|=bits[i+2]<<(64-s))
    v&_fmask(b)
end
@inline function _fput(bits,offset,b,v::UInt64)
    i,s=divrem(offset,64)
    bits=_replace(bits,i+1,(bits[i+1]&~(_fmask(b)<<s))|(v<<s))
    if s+b>64
        m=_fmask(s+b-64)
        bits=_replace(bits,i+2,(bits[i+2]&~m)|((v>>(64-s))&m))
    end
    bits
end
@inline _ftoken(bits,i,b) = _fget(bits,(i-1)*b,b)
@inline _fsettoken(bits,i,b,v) = _fput(bits,(i-1)*b,b,v)
@inline _ftoken(w::FunctionalWord{F,A,B},i) where {F,A,B} = _ftoken(w.bits,i,B)
@inline function degree(w::W) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    d=_fcap(W)
    Int(_fget(w.bits,d*(B+1),_fheader(d)))
end
@inline function _ftail(w::W) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    d=_fcap(W); h=_fheader(d)
    Int(_fget(w.bits,d*(B+1)+h,h))
end
@inline function _flengths(w::W,n,t) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    d=_fcap(W); h=_fheader(d)
    bits=_fput(w.bits,d*(B+1),h,UInt64(n))
    W(_fput(bits,d*(B+1)+h,h,UInt64(t)),Val(:normal))
end
@inline _fend(w::W,i) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}} = !iszero(_fget(w.bits,_fcap(W)*B+i-1,1))
@inline _fmark(w::W,i) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}} = W(_fput(w.bits,_fcap(W)*B+i-1,1,UInt64(1)),Val(:normal))
@inline _fsite(v,s) = Int(v&_fmask(s))
@inline _fzero(::Type{<:FunctionalWord{F,A,B,S,K}}) where {F,A,B,S,K} = ntuple(_->UInt64(0),Val(K))
Base.one(::Type{W}) where {W<:FunctionalWord} = W(_fzero(W),Val(:normal))
Base.isone(w::FunctionalWord) = iszero(degree(w))
Base.:(==)(a::W,b::W) where {W<:FunctionalWord} = a.bits==b.bits
Base.hash(w::FunctionalWord,h::UInt) = hash(w.bits,h)
function Base.isless(a::W,b::W) where {W<:FunctionalWord}
    degree(a)!=degree(b) && return degree(a)<degree(b)
    _ftail(a)!=_ftail(b) && return _ftail(a)<_ftail(b)
    a.bits<b.bits
end
ProductStyle(::Type{<:FunctionalWord}) = SingleTerm()

function word_type(a::FunctionalAlgebra{F,A},nsites::Integer,max_letters::Integer) where {F,A}
    nsites>=0 && max_letters>=0 || throw(ArgumentError("Site count and required degree must be nonnegative."))
    s=max(1,ndigits(nsites;base=2)); b=s+ndigits(a.base.nletters;base=2)
    cap=_functional_capacity(b,8)
    b<=64 && max_letters<=cap || throw(ArgumentError("Functional word requires degree $max_letters with $b-bit tokens; capacity is $cap letters in eight limbs (64 bytes); no larger word type is supported."))
    for k in 1:8
        max_letters<=_functional_capacity(b,k) && return FunctionalWord{F,A,b,s,k}
    end
    error("unreachable functional width")
end
@noinline function _foverflow(::Type{W},required) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    throw(ArgumentError("Functional word requires degree $required with $B-bit tokens; capacity is $(_fcap(W)) letters in $K limbs. Select a wider word_type before multiplication."))
end
@inline function _fappend(w::W,v::UInt64) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    n=degree(w)+1
    n<=_fcap(W) || _foverflow(W,n)
    _flengths(W(_fsettoken(w.bits,n,B,v),Val(:normal)),n,_ftail(w))
end
function site_letter(::Type{W},site::Integer,label::Integer) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    1<=site<=_fmask(S) && 1<=label<=_fmask(B-S) || throw(ArgumentError("Site or label exceeds the functional token width."))
    _flengths(_fappend(one(W),(UInt64(label)<<S)|UInt64(site)),1,1)
end

@inline function _fstop(w,i)
    j=i
    while !_fend(w,j); j+=1; end
    j
end
@inline function _fslice(w::W,i,j) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    bits=_fzero(W)
    for k in i:j; bits=_fsettoken(bits,k-i+1,B,_ftoken(w,k)); end
    bits
end
@inline function _frawless(a,na,b,nb,B)
    na!=nb && return na<nb
    for i in 1:na
        x,y=_ftoken(a,i,B),_ftoken(b,i,B)
        x!=y && return x<y
    end
    false
end
function _freverse(bits,n,B,S)
    out=map(zero,bits); start=1
    while start<=n
        stop=start
        while stop<n && _fsite(_ftoken(bits,stop+1,B),S)==_fsite(_ftoken(bits,start,B),S); stop+=1; end
        for i in start:stop; out=_fsettoken(out,i,B,_ftoken(bits,stop-i+start,B)); end
        start=stop+1
    end
    out
end
function _fcyclic(bits,n,::Type{A},B,S) where {A}
    out=map(zero,bits); used=0; start=1
    while start<=n
        stop=start
        while stop<n && _fsite(_ftoken(bits,stop+1,B),S)==_fsite(_ftoken(bits,start,B),S); stop+=1; end
        lo,hi=start,stop
        if A!==Free
            while lo<hi && _ftoken(bits,lo,B)==_ftoken(bits,hi,B)
                hi-=1
                A===Unipotent && (lo+=1)
            end
        end
        len=hi-lo+1; best=0
        for shift in 1:len-1
            for j in 0:len-1
                x=_ftoken(bits,lo+mod(shift+j,len),B)
                y=_ftoken(bits,lo+mod(best+j,len),B)
                x==y && continue
                x<y && (best=shift)
                break
            end
        end
        for j in 0:len-1
            used+=1; out=_fsettoken(out,used,B,_ftoken(bits,lo+mod(best+j,len),B))
        end
        start=stop+1
    end
    out,used
end
function _ffactor(bits,n,::Type{F},::Type{A},B,S) where {F,A}
    back=_freverse(bits,n,B,S)
    if F===RealTrace
        bits,m=_fcyclic(bits,n,A,B,S)
        back,k=_fcyclic(back,n,A,B,S)
        return _frawless(back,k,bits,m,B) ? (back,k) : (bits,m)
    end
    _frawless(back,n,bits,n,B) ? (back,n) : (bits,n)
end
function _fcopyfactor(out,w,i,j)
    for k in i:j; out=_fappend(out,_ftoken(w,k)); end
    _fmark(out,degree(out))
end
function _finsert(w::W,bits,n) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    n==0 && return w
    degree(w)+n<=_fcap(W) || _foverflow(W,degree(w)+n)
    out=one(W)
    for i in 1:_ftail(w); out=_fappend(out,_ftoken(w,i)); end
    out=_flengths(out,degree(out),degree(out))
    i=_ftail(w)+1; inserted=false
    while i<=degree(w)
        j=_fstop(w,i)
        if !inserted && _frawless(bits,n,_fslice(w,i,j),j-i+1,B)
            for k in 1:n; out=_fappend(out,_ftoken(bits,k,B)); end
            out=_fmark(out,degree(out)); inserted=true
        end
        out=_fcopyfactor(out,w,i,j); i=j+1
    end
    if !inserted
        for k in 1:n; out=_fappend(out,_ftoken(bits,k,B)); end
        out=_fmark(out,degree(out))
    end
    out
end

function mul(a::W,b::W) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    # Twice the payload capacity is needed transiently: a long prefix can cancel.
    raw=ntuple(_->UInt64(0),Val(2K)); n=0; i=j=1
    na,nb=_ftail(a),_ftail(b)
    while i<=na || j<=nb
        if j>nb || (i<=na && _fsite(_ftoken(a,i),S)<=_fsite(_ftoken(b,j),S))
            v=_ftoken(a,i); i+=1
        else
            v=_ftoken(b,j); j+=1
        end
        if A!==Free && n>0 && _ftoken(raw,n,B)==v
            A===Unipotent && (raw=_fsettoken(raw,n,B,UInt64(0)); n-=1)
        else
            n+=1; raw=_fsettoken(raw,n,B,v)
        end
    end
    required=n+degree(a)-na+degree(b)-nb
    required<=_fcap(W) || _foverflow(W,required)
    out=one(W)
    for k in 1:n; out=_fappend(out,_ftoken(raw,k,B)); end
    out=_flengths(out,n,n); i=na+1; j=nb+1
    while i<=degree(a) || j<=degree(b)
        ai=i<=degree(a) ? _fstop(a,i) : i-1
        bj=j<=degree(b) ? _fstop(b,j) : j-1
        if j>degree(b) || (i<=degree(a) && !_frawless(_fslice(b,j,bj),bj-j+1,_fslice(a,i,ai),ai-i+1,B))
            out=_fcopyfactor(out,a,i,ai); i=ai+1
        else
            out=_fcopyfactor(out,b,j,bj); j=bj+1
        end
    end
    one(Phase),out
end
function star(w::W) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    raw=_freverse(w.bits,_ftail(w),B,S); bits=w.bits
    for i in 1:_ftail(w); bits=_fsettoken(bits,i,B,_ftoken(raw,i,B)); end
    one(Phase),W(bits,Val(:normal))
end
function moment_key(w::W) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    _ftail(w)==0 && return w,one(Phase),false
    raw,n=_ffactor(_fslice(w,1,_ftail(w)),_ftail(w),F,A,B,S)
    out=one(W); i=_ftail(w)+1
    while i<=degree(w)
        j=_fstop(w,i); out=_fcopyfactor(out,w,i,j); i=j+1
    end
    _finsert(out,raw,n),one(Phase),false
end

struct FunctionalSupport{W}
    word::W
end
Base.eltype(::Type{FunctionalSupport{W}}) where {W} = W
Base.length(x::FunctionalSupport) = degree(x.word)
function Base.iterate(x::FunctionalSupport{W},i=1) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    i>degree(x.word) && return nothing
    _flengths(_fappend(one(W),_ftoken(x.word,i)),1,1),i+1
end
support_letters(w::FunctionalWord) = FunctionalSupport(w)
function word_letters(w::FunctionalWord)
    _ftail(w)==degree(w) || throw(ArgumentError("A word with scalar factors is not a product of operator letters; use support_letters for its dependencies."))
    support_letters(w)
end
struct FunctionalSites{W}
    word::W
end
Base.eltype(::Type{<:FunctionalSites}) = Int
Base.length(x::FunctionalSites) = count(_->true,x)
function Base.iterate(x::FunctionalSites{W},previous=0) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    next=typemax(Int)
    for i in 1:degree(x.word)
        site=_fsite(_ftoken(x.word,i),S)
        previous<site<next && (next=site)
    end
    next==typemax(Int) ? nothing : (next,next)
end
sites(w::FunctionalWord) = FunctionalSites(w)

function _check_functional_word(w::W) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    n,t=degree(w),_ftail(w)
    0<=t<=n<=_fcap(W) || throw(ArgumentError("Invalid functional word lengths."))
    rebuilt=one(W); start=1
    for i in 1:n
        v=_ftoken(w,i)
        _fsite(v,S)>0 && v>>S>0 || throw(ArgumentError("Functional tokens need nonzero sites and labels."))
        if i>start
            prev=_ftoken(w,i-1)
            _fsite(prev,S)<=_fsite(v,S) && (A===Free || prev!=v) || throw(ArgumentError("Functional segment is not in base normal form."))
        end
        if i<=t
            !_fend(w,i) || throw(ArgumentError("Operator tail cannot have scalar endpoints."))
            rebuilt=_fappend(rebuilt,v)
            i==t && (rebuilt=_flengths(rebuilt,t,t); start=i+1)
        elseif _fend(w,i)
            raw,m=_ffactor(_fslice(w,start,i),i-start+1,F,A,B,S)
            rebuilt=_finsert(rebuilt,raw,m); start=i+1
        end
    end
    rebuilt.bits==w.bits || throw(ArgumentError("Functional word is not canonical, has missing endpoints or nonzero padding."))
    nothing
end
