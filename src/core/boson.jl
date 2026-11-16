"""
    algebra = Bosonic()

Return the bosonic algebra description.
Its operators satisfy `b[i]*b[j] = b[j]*b[i]` and
`b[i]*bd[j] - bd[j]*b[i] = δ[i,j]`, where `bd[i] = b[i]'`.
The symbol `δ[i,j]` is 1 for equal indices and 0 otherwise.
Use `bosonic_variables` to get a registry and operator vectors.

Example: `using NCTSSoS; registry, (b, bd) = bosonic_variables(2)`.
"""
struct Bosonic end

"Normal-ordered bosons: creation and annihilation exponents in B-bit lanes."
struct BosonWord{B,K} <: AbstractWord
    c::NTuple{K,UInt64}
    a::NTuple{K,UInt64}
    function BosonWord{B,K}(c::NTuple{K,UInt64},a::NTuple{K,UInt64}) where {B,K}
        1<=B<=8 && K>0 || throw(ArgumentError("Bosonic storage needs 1:8 bits per exponent and at least one limb."))
        mask=typemax(UInt64) >> (64-B*(64÷B))
        all(x->iszero(x & ~mask),c) && all(x->iszero(x & ~mask),a) ||
            throw(ArgumentError("Bosonic word has nonzero unused lane bits."))
        new{B,K}(c,a)
    end
end
function word_type(::Bosonic,n::Integer,d::Integer)
    n>=0 && 0<=d<=255 || throw(ArgumentError("Bosonic word capacity requires nonnegative modes and degree at most 255."))
    b=max(1,ndigits(d;base=2))
    return BosonWord{b,max(1,cld(n,64÷b))}
end
Base.one(::Type{BosonWord{B,K}}) where {B,K} = BosonWord{B,K}(ntuple(_->UInt64(0),Val(K)),ntuple(_->UInt64(0),Val(K)))
Base.isone(w::BosonWord) = all(iszero,w.c) && all(iszero,w.a)
Base.:(==)(u::W,v::W) where {W<:BosonWord} = u.c==v.c && u.a==v.a
Base.hash(w::BosonWord,h::UInt) = hash(w.a,hash(w.c,h))
Base.isless(u::W,v::W) where {W<:BosonWord} = (degree(u),u.c,u.a)<(degree(v),v.c,v.a)
@inline _boson_modes(::Type{BosonWord{B,K}}) where {B,K} = K*(64÷B)
@inline function _boson_count(::Type{BosonWord{B,K}},pack,mode) where {B,K}
    limb,lane=divrem(mode-1,64÷B)
    return Int((pack[limb+1]>>(B*lane)) & ((UInt64(1)<<B)-1))
end
@inline function _boson_set(::Type{BosonWord{B,K}},pack,mode,count) where {B,K}
    0<=count<(1<<B) || throw(ArgumentError("Bosonic exponent $count exceeds the $B-bit capacity $((1<<B)-1); select a wider word_type or reduce the model degree."))
    limb,lane=divrem(mode-1,64÷B)
    shift=B*lane; mask=((UInt64(1)<<B)-1)<<shift
    return ntuple(i->i==limb+1 ? (pack[i]&~mask)|(UInt64(count)<<shift) : pack[i],Val(K))
end
degree(w::W) where {W<:BosonWord} = sum(_boson_count(W,w.c,i)+_boson_count(W,w.a,i) for i in 1:_boson_modes(W))
star(w::W) where {W<:BosonWord} = (one(Phase),W(w.a,w.c))
ProductStyle(::Type{<:BosonWord}) = MultiTerm()

struct BosonSites{W}
    word::W
end
Base.eltype(::Type{<:BosonSites}) = Int
Base.length(s::BosonSites) = count(_->true,s)
function Base.iterate(s::BosonSites{W},mode=1) where W
    while mode<=_boson_modes(W)
        if _boson_count(W,s.word.c,mode)+_boson_count(W,s.word.a,mode)>0
            return mode,mode+1
        end
        mode+=1
    end
    return nothing
end
sites(w::BosonWord) = BosonSites(w)
function boson_letter(::Type{W},index::Integer) where {W<:BosonWord}
    1<=abs(index)<=_boson_modes(W) || throw(ArgumentError("Bosonic mode is outside the word capacity."))
    empty=one(W)
    pack=_boson_set(W,empty.c,abs(index),1)
    return index<0 ? W(pack,empty.a) : W(empty.c,pack)
end
word_letters(w::BosonWord) = WordLetters(w)
function Base.iterate(ls::WordLetters{W},state=(_boson_modes(W),0,true)) where {W<:BosonWord}
    mode,used,creating=state
    while true
        if creating && mode==0
            mode=1; used=0; creating=false
        end
        !creating && mode>_boson_modes(W) && return nothing
        count=_boson_count(W,creating ? ls.word.c : ls.word.a,mode)
        used<count && return boson_letter(W,creating ? -mode : mode),(mode,used+1,creating)
        mode+=creating ? -1 : 1
        used=0
    end
end
legacy_letters(w::W) where {W<:BosonWord} = [isone(W(l.c,one(W).a)) ? only(sites(l)) : -only(sites(l)) for l in word_letters(w)]
legacy_terms(w::BosonWord) = [(1.0+0im,legacy_letters(w))]
function word_from_legacy(::Type{W},letters) where {W<:BosonWord}
    out=one(W)
    issorted(letters) || throw(ArgumentError("Bosonic letters must be normal ordered; multiply raw letters with mul! first."))
    for index in letters
        l=boson_letter(W,index)
        out=_boson_add(out,l)
    end
    return out
end
function _boson_add(u::W,v::W) where {W<:BosonWord}
    c,a=u.c,u.a
    for i in 1:_boson_modes(W)
        c=_boson_set(W,c,i,_boson_count(W,u.c,i)+_boson_count(W,v.c,i))
        a=_boson_set(W,a,i,_boson_count(W,u.a,i)+_boson_count(W,v.a,i))
    end
    return W(c,a)
end
function _boson_scale(c::T,p,q,k) where {T<:Number}
    if T<:Integer || (T<:Complex && typeof(real(c))<:Integer)
        try
            I=typeof(real(c))
            F=promote_type(I,Int)
            factor=binomial(convert(F,p),k)
            for j in 0:k-1
                factor=Base.Checked.checked_mul(factor,convert(F,q-j))
            end
            factor=convert(I,factor)
            if T<:Integer
                return Base.Checked.checked_mul(c,factor)
            end
            return T(Base.Checked.checked_mul(real(c),factor),Base.Checked.checked_mul(imag(c),factor))
        catch error
            error isa Union{OverflowError,InexactError} || rethrow()
            throw(ArgumentError("Bosonic product coefficient exceeds $T; use a wider coefficient type."))
        end
    end
    out=c
    for j in 1:k
        out=convert(T,out*((p-j+1)*(q-j+1))/j)
        isfinite(out) || throw(ArgumentError("Bosonic product coefficient exceeds $T; use a wider coefficient type."))
    end
    return out
end
function _boson_contract!(sink,u::W,v::W,out::W,mode,c) where {W<:BosonWord}
    while mode<=_boson_modes(W)
        p,q=_boson_count(W,u.a,mode),_boson_count(W,v.c,mode)
        if min(p,q)>0
            nc,na=_boson_count(W,out.c,mode),_boson_count(W,out.a,mode)
            for k in 0:min(p,q)
                w=W(_boson_set(W,out.c,mode,nc-k),_boson_set(W,out.a,mode,na-k))
                _boson_contract!(sink,u,v,w,mode+1,_boson_scale(c,p,q,k))
            end
            return sink
        end
        mode+=1
    end
    push_term!(sink,c,out)
    return sink
end
function _mul!(::MultiTerm,sink,u::W,v::W,c::Number) where {W<:BosonWord}
    out=_boson_add(u,v) # Check every uncontracted exponent before emitting terms.
    return _boson_contract!(sink,u,v,out,1,c)
end
function _pack_word(::Type{W},w::BosonWord) where {W<:BosonWord}
    c,a=one(W).c,one(W).a
    for i in sites(w)
        i<=_boson_modes(W) || throw(ArgumentError("Bosonic repacking exceeds the destination mode capacity."))
        c=_boson_set(W,c,i,_boson_count(typeof(w),w.c,i))
        a=_boson_set(W,a,i,_boson_count(typeof(w),w.a,i))
    end
    return W(c,a)
end
function _extend_rows!(rows::Vector{W},::Bosonic,group,w::W,start,remaining) where {W<:BosonWord}
    remaining==0 && return rows
    n=length(group)
    for j in start:2n
        next=_boson_add(w,boson_letter(W,j<=n ? -group[n-j+1] : group[j-n]))
        push!(rows,next)
        _extend_rows!(rows,Bosonic(),group,next,j,remaining-1)
    end
    return rows
end
