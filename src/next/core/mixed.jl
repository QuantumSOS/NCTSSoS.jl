"Tensor components with disjoint site ranges; distinct components commute."
struct MixedAlgebra{A,N}
    components::A
end

"A tensor word stores only its packed components; site offsets live in the type."
struct MixedWord{T<:Tuple,O} <: AbstractWord
    components::T
    function MixedWord{T,O}(parts::T) where {T<:Tuple,O}
        isbitstype(T) && all(w->w isa Union{PauliWord,FermiWord,BosonWord,SiteWord},parts) ||
            throw(ArgumentError("Mixed words need packed ordinary components; state/trace components are unsupported."))
        !isempty(parts) && length(O)==length(parts)+1 && first(O)==0 && all(i->O[i]<O[i+1],1:length(parts)) ||
            throw(ArgumentError("Mixed word offsets must describe positive component site counts."))
        for i in eachindex(parts)
            all(s->1<=s<=O[i+1]-O[i],sites(parts[i])) || throw(ArgumentError(
                "Mixed component $i has a site outside its declared range; use the component's variables."))
        end
        new{T,O}(parts)
    end
    MixedWord{T,O}(parts::T,::Val{:normal}) where {T<:Tuple,O} = new{T,O}(parts)
end
_mixed_offsets(sizes::Tuple) = (0,accumulate(+,sizes)...)
_real_structure(a::MixedAlgebra) = all(_real_structure,a.components)
_real_structure(::Type{MixedWord{T,O}}) where {T,O} = all(_real_structure,fieldtypes(T))
function word_type(a::MixedAlgebra{A,N},nsites::Integer,d::Integer) where {A,N}
    nsites==sum(N) || throw(ArgumentError("Mixed site count must equal the sum of its component counts."))
    types=map((algebra,n)->word_type(algebra,n,d),a.components,N)
    return MixedWord{Tuple{types...},_mixed_offsets(N)}
end
Base.one(::Type{MixedWord{T,O}}) where {T,O} = MixedWord{T,O}(map(one,fieldtypes(T)),Val(:normal))
Base.isone(w::MixedWord) = all(isone,w.components)
Base.:(==)(a::W,b::W) where {W<:MixedWord} = a.components==b.components
Base.hash(w::MixedWord,h::UInt) = hash(w.components,h)
degree(w::MixedWord) = sum(map(degree,w.components))
Base.isless(a::W,b::W) where {W<:MixedWord} =
    degree(a)==degree(b) ? isless(a.components,b.components) : degree(a)<degree(b)
function star(w::W) where {W<:MixedWord}
    images=map(star,w.components)
    return prod(map(first,images)),W(map(last,images),Val(:normal))
end
@generated function ProductStyle(::Type{MixedWord{T,O}}) where {T,O}
    return all(t->t<:Union{PauliWord,SiteWord},fieldtypes(T)) ? :(SingleTerm()) : :(MultiTerm())
end
function mul(a::W,b::W) where {W<:MixedWord}
    images=map(mul,a.components,b.components)
    return prod(map(first,images)),W(map(last,images),Val(:normal))
end

struct MixedProductSink{W,I,S}
    sink::S
    left::W
    right::W
    partial::W
end
@inline function _mixed_product!(sink,a::W,b::W,partial::W,::Val{I},c) where {T,O,W<:MixedWord{T,O},I}
    next=MixedProductSink{W,I,typeof(sink)}(sink,a,b,partial)
    _mul!(ProductStyle(fieldtype(T,I)),next,a.components[I],b.components[I],c)
    return sink
end
@inline function push_term!(s::MixedProductSink{W,I},c::Number,w::AbstractWord) where {T,O,W<:MixedWord{T,O},I}
    partial=W(ntuple(j->j==I ? w : s.partial.components[j],Val(fieldcount(T))),Val(:normal))
    if I==fieldcount(T)
        push_term!(s.sink,c,partial)
    else
        _mixed_product!(s.sink,s.left,s.right,partial,Val(I+1),c)
    end
    return s
end
mul!(sink,a::W,b::W,c::Number=1) where {W<:MixedWord} =
    _mixed_product!(sink,a,b,one(W),Val(1),c)

function _mixed_lift(::Type{W},::Val{I},part) where {T,O,W<:MixedWord{T,O},I}
    return W(ntuple(j->j==I ? part : one(fieldtype(T,j)),Val(fieldcount(T))),Val(:normal))
end
struct MixedSites{W}
    word::W
end
Base.eltype(::Type{<:MixedSites}) = Int
Base.length(s::MixedSites) = sum(map(w->length(sites(w)),s.word.components))
sites(w::MixedWord) = MixedSites(w)
word_letters(w::MixedWord) = WordLetters(w)
_mixed_iterator(::Val{:sites},w) = sites(w)
_mixed_iterator(::Val{:letters},w) = word_letters(w)
# Keep every component's iterator state in a fixed-type tuple, including empty
# components. A one-letter seed supplies its state type; fresh iteration ignores it.
_mixed_seed(kind,::W) where {W} = last(iterate(_mixed_iterator(kind,_mixed_letter(W,1,1))))
_mixed_initial(w,kind) = (1,true,map(part->_mixed_seed(kind,part),w.components))
@generated function _mixed_next(w::W,kind::Val{Kind},active,fresh,states) where {T,O,W<:MixedWord{T,O},Kind}
    steps=map(1:fieldcount(T)) do i
        value=Kind===:sites ? :(item[1]+$(O[i])) : :(_mixed_lift(W,Val($i),item[1]))
        updated=Expr(:tuple,(j==i ? :(item[2]) : :(states[$j]) for j in 1:fieldcount(T))...)
        quote
            if active<=$i
                iterator=_mixed_iterator(kind,w.components[$i])
                item=active==$i && !fresh ? iterate(iterator,states[$i]) : iterate(iterator)
                if !isnothing(item)
                    return $value,($i,false,$updated)
                end
            end
        end
    end
    return quote
        $(steps...)
        nothing
    end
end
Base.iterate(s::MixedSites,state=_mixed_initial(s.word,Val(:sites))) =
    _mixed_next(s.word,Val(:sites),state...)
Base.iterate(ls::WordLetters{<:MixedWord},state=_mixed_initial(ls.word,Val(:letters))) =
    _mixed_next(ls.word,Val(:letters),state...)
function _pack_word(::Type{MixedWord{T,O}},w::MixedWord{S,O}) where {T,S,O}
    return MixedWord{T,O}(map(_pack_word,fieldtypes(T),w.components),Val(:normal))
end

function _check_label(a::MixedAlgebra{A,N},site,label) where {A,N}
    offset=0
    for (component,n) in zip(a.components,N)
        site<=offset+n && return _check_label(component,label)
        offset+=n
    end
    throw(ArgumentError("Mixed operator site is outside the component ranges."))
end
function _mixed_letter(::Type{W},site,label) where {W<:PauliWord}
    pauli_letter(W,site,label)
end
_mixed_letter(::Type{W},site,label) where {W<:FermiWord} = fermi_letter(W,site*label)
_mixed_letter(::Type{W},site,label) where {W<:BosonWord} = boson_letter(W,site*label)
_mixed_letter(::Type{W},site,label) where {W<:SiteWord} = site_letter(W,site,label)
function word(r::Registry{<:MixedAlgebra},name::Symbol,::Type{W}) where {T,O,W<:MixedWord{T,O}}
    site,label=r[name]
    return W(ntuple(i->O[i]<site<=O[i+1] ? _mixed_letter(fieldtype(T,i),site-O[i],label) :
        one(fieldtype(T,i)),Val(fieldcount(T))))
end

function _mixed_rows!(out,::Type{W},rows::Tuple,prefix,remaining) where {W}
    for w in first(rows)
        d=degree(w)
        d>remaining && break
        _mixed_rows!(out,W,Base.tail(rows),(prefix...,w),remaining-d)
    end
    return out
end
_mixed_rows!(out,::Type{W},::Tuple{},prefix,remaining) where {W} = push!(out,W(prefix))
function basis_rows(::Type{W},a::MixedAlgebra{A,N},group,order::Integer) where {T,O,W<:MixedWord{T,O},A,N}
    order>=0 || throw(ArgumentError("Mixed basis order must be nonnegative."))
    selected=collect(Int,group)
    allunique(selected) && all(s->1<=s<=last(O),selected) ||
        throw(ArgumentError("Mixed basis sites must be distinct and inside the registry."))
    rows=ntuple(i->basis_rows(fieldtype(T,i),a.components[i],
        [s-O[i] for s in selected if O[i]<s<=O[i+1]],order),Val(fieldcount(T)))
    return sort!(_mixed_rows!(W[],W,rows,(),order))
end
