"""
    LetterPermutation(a => b, b => a)
    LetterPermutation(pairs)

An optionally signed permutation of operator letters from one Free, Projector
or Unipotent registry, ordinary or functional. Omitted letters are fixed.
Sources have coefficient one; targets have coefficient +1 or -1. Negative
targets are allowed for Free and Unipotent, but not Projector: -P is not a
projector. The underlying letter map must be bijective and send all
letters of each site to one site, preserving the commuting-site partition.

For example, exchange two parties and reverse their setting labels:
```julia
registry, (A, B) = site_variables(Unipotent(), 2, 2)
exchange = LetterPermutation(A[1] => B[2], A[2] => B[1],
                             B[1] => A[2], B[2] => A[1])
problem = Problem(A[1]*B[2] + A[2]*B[1]; symmetry=(exchange,))
```

`compose(g, h)` applies `h` first, then `g`, as for the existing actions.
`LetterPermutation(W, g::SitePermutation)` converts a whole-site action for
the word type W without choosing a word-storage width for the action itself.
"""
struct LetterPermutation{A<:SiteRule} <: AbstractAction
    site_map::SitePermutation
    label_maps::Vector{SitePermutation}
    negative::Set{Tuple{Int,Int}}
    names::Dict{Symbol,Tuple{Int,Int}}
    function LetterPermutation{A}(site_map,label_maps,negative,names,::Val{:normal}) where {A<:SiteRule}
        while !isempty(label_maps) && isone(last(label_maps)); pop!(label_maps); end
        new{A}(site_map,label_maps,negative,names)
    end
end

_letter_rule(::Type{<:SiteWord{A}}) where {A} = A
_letter_rule(::Type{<:FunctionalWord{F,A}}) where {F,A} = A
_letter_location(w::W) where {A,U,B,S,W<:SiteWord{A,U,B,S}} = (_site(W,_letter(w,1)),Int(_letter(w,1)>>S))
_letter_location(w::FunctionalWord{F,A,B,S}) where {F,A,B,S} = (_fsite(_ftoken(w,1),S),Int(_ftoken(w,1)>>S))

function LetterPermutation(::Type{W},g::SitePermutation) where {W<:Union{SiteWord,FunctionalWord}}
    LetterPermutation{_letter_rule(W)}(g,SitePermutation[],Set{Tuple{Int,Int}}(),Dict{Symbol,Tuple{Int,Int}}(),Val(:normal))
end
@inline _label_image(g::LetterPermutation,s,l) = s<=length(g.label_maps) ? _image(g.label_maps[s],l) : l
@inline _letter_image(g::LetterPermutation,s,l) = (_image(g.site_map,s),_label_image(g,s,l))
@inline _letter_image(g::SitePermutation,s,l) = (_image(g,s),l)
@inline _letter_negative(g::LetterPermutation,s,l) = (s,l) in g.negative
@inline _letter_negative(::SitePermutation,s,l) = false
Base.one(g::LetterPermutation{A}) where {A} = LetterPermutation{A}(SitePermutation(Int[]),SitePermutation[],Set{Tuple{Int,Int}}(),g.names,Val(:normal))
Base.isone(g::LetterPermutation) = isone(g.site_map) && isempty(g.label_maps) && isempty(g.negative)
Base.:(==)(g::LetterPermutation,h::LetterPermutation) = typeof(g)===typeof(h) && g.site_map==h.site_map && g.label_maps==h.label_maps && g.negative==h.negative
Base.hash(g::LetterPermutation{A},h::UInt) where {A} = hash(A,hash(g.site_map,hash(g.label_maps,hash(g.negative,h))))
_action_sites(g::LetterPermutation) = max(length(g.site_map.images),length(g.label_maps),maximum(first,g.negative;init=0))
_label_count(g::LetterPermutation,s) = max(s<=length(g.label_maps) ? length(g.label_maps[s].images) : 0,
    maximum((l for (site,l) in g.negative if site==s);init=0))
function compose(g::LetterPermutation{A},h::LetterPermutation{A}) where {A}
    labels=SitePermutation[]; negative=Set{Tuple{Int,Int}}()
    for s in 1:max(_action_sites(g),_action_sites(h))
        hs=_image(h.site_map,s)
        n=max(_label_count(h,s),_label_count(g,hs))
        push!(labels,SitePermutation([_label_image(g,hs,_label_image(h,s,l)) for l in 1:n]))
        for l in 1:n
            xor(_letter_negative(h,s,l),_letter_negative(g,hs,_label_image(h,s,l))) && push!(negative,(s,l))
        end
    end
    # Actions built from one registry share names, including during group closure.
    names=g.names===h.names || isempty(h.names) ? g.names : isempty(g.names) ? h.names : merge(h.names,g.names)
    LetterPermutation{A}(compose(g.site_map,h.site_map),labels,negative,names,Val(:normal))
end
function Base.inv(g::LetterPermutation{A}) where {A}
    labels=[SitePermutation(Int[]) for _ in 1:_action_sites(g)]
    for s in eachindex(g.label_maps); labels[_image(g.site_map,s)]=inv(g.label_maps[s]); end
    negative=Set(_letter_image(g,s,l) for (s,l) in g.negative)
    LetterPermutation{A}(inv(g.site_map),labels,negative,g.names,Val(:normal))
end

function act(g::LetterPermutation{A},w::W) where {A,U,B,S,W<:SiteWord{A,U,B,S}}
    bits=zero(U); sign=false
    for i in 1:degree(w)
        token=_letter(w,i)
        sign=xor(sign,_letter_negative(g,_site(W,token),Int(token>>S)))
        site,label=_letter_image(g,_site(W,token),Int(token>>S))
        mapped=site_letter(W,site,label).bits
        j=1
        while j<i && _site(W,_letter(W(bits,Val(:normal)),j))<=site; j+=1; end
        shift=B*(j-1); low=(one(U)<<shift)-one(U)
        bits=(bits&low)|((bits&~low)<<B)|(mapped<<shift)
    end
    Phase(2sign),W(bits,Val(:normal))
end

# Stable insertion by site keeps the internal order of each mapped site word.
# Bijectivity preserves adjacent inequalities, so no base reductions are lost.
function _mapped_segment(g,w::W,lo,hi) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    raw=_fzero(W); n=0; sign=false
    for i in lo:hi
        token=_ftoken(w,i)
        sign=xor(sign,_letter_negative(g,_fsite(token,S),Int(token>>S)))
        site,label=_letter_image(g,_fsite(token,S),Int(token>>S))
        mapped=_ftoken(site_letter(W,site,label),1)
        j=n
        while j>0 && _fsite(_ftoken(raw,j,B),S)>site
            raw=_fsettoken(raw,j+1,B,_ftoken(raw,j,B)); j-=1
        end
        raw=_fsettoken(raw,j+1,B,mapped); n+=1
    end
    raw,n,sign
end
function _act_functional(g,w::W) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    raw,n,sign=_mapped_segment(g,w,1,_ftail(w))
    out=one(W)
    for i in 1:n; out=_fappend(out,_ftoken(raw,i,B)); end
    out=_flengths(out,n,n)
    i=_ftail(w)+1
    while i<=degree(w)
        j=_fstop(w,i)
        raw,n,factor_sign=_mapped_segment(g,w,i,j)
        sign=xor(sign,factor_sign)
        raw,n=_ffactor(raw,n,F,A,B,S)
        out=_finsert(out,raw,n); i=j+1
    end
    Phase(2sign),out
end
act(g::LetterPermutation{A},w::FunctionalWord{F,A}) where {F,A} = _act_functional(g,w)
act(g::SitePermutation,w::FunctionalWord) = _act_functional(g,w)

_letter_name(names,location) = something(findfirst(==(location),names),Symbol("site",location[1],"_",location[2]))
function Base.show(io::IO,g::LetterPermutation)
    print(io,"LetterPermutation(")
    if isempty(g.names) && !isone(g)
        print(io,"sites=",g.site_map.images,", settings=",[p.images for p in g.label_maps])
        isempty(g.negative) || print(io,", negative=",sort!(collect(g.negative)))
    else
        separator=""
        for (name,source) in sort!(collect(g.names);by=last)
            target=_letter_image(g,source...)
            negative=_letter_negative(g,source...)
            source==target && !negative && continue
            print(io,separator,name," => ",negative ? "-" : "",_letter_name(g.names,target)); separator=", "
        end
    end
    print(io,")")
end
