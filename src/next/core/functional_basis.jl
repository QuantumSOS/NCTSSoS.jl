function _pack_word(::Type{W},w::FunctionalWord{F,A,B,S}) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    degree(w)<=_fcap(W) || _foverflow(W,degree(w))
    out=one(W)
    for i in 1:degree(w)
        out=_fappend(out,_ftoken(w,i))
        _fend(w,i) && (out=_fmark(out,i))
    end
    _flengths(out,degree(w),_ftail(w))
end

function _functional_operators!(rows,letters,w::W,start,remaining) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    remaining==0 && return rows
    previous_site=0; group_start=start
    for j in start:length(letters)
        v=letters[j]; site=_fsite(_ftoken(v,1),S)
        if site!=previous_site; group_start=j; previous_site=site; end
        A!==Free && !isone(w) && _ftoken(w,degree(w))==_ftoken(v,1) && continue
        next=mul(w,v)[2]
        push!(rows,next)
        _functional_operators!(rows,letters,next,group_start,remaining-1)
    end
    rows
end
function _functional_rows!(rows,operators,symbols,scalar,start,order)
    for op in operators
        degree(scalar)+degree(op)>order && break
        push!(rows,mul(scalar,op)[2])
    end
    for i in start:length(symbols)
        degree(scalar)+degree(symbols[i])>order && break
        _functional_rows!(rows,operators,symbols,mul(scalar,symbols[i])[2],i,order)
    end
    rows
end
"""
Complete mixed rows, graded by all operator letters including those inside
symbols. `group` is either distinct site indices or distinct operator letters
of type W. A letter set is generated directly, without a global-basis filter.
"""
function basis_rows(::Type{W},algebra::FunctionalAlgebra,group,order::Integer) where {F,A,W<:FunctionalWord{F,A}}
    0<=order<=_fcap(W) || throw(ArgumentError("Row order must be nonnegative and within the functional word capacity $(_fcap(W))."))
    algebra isa FunctionalAlgebra{F,A} || throw(ArgumentError("Row letters and algebra must have the same functional domain and base rule."))
    selected=collect(group)
    if all(x->x isa Integer,selected)
        allunique(selected) && all(>(0),selected) || throw(ArgumentError("Row sites must be distinct positive integers."))
        selected=W[site_letter(W,s,l) for s in selected for l in 1:algebra.base.nletters]
    end
    all(w->w isa W && degree(w)==1 && _ftail(w)==1 && _letter_location(w)[2]<=algebra.base.nletters,selected) ||
        throw(ArgumentError("Row entries must be operator letters of type $W, not scalar symbols, identities or longer words."))
    allunique(selected) || throw(ArgumentError("Row operator letters must be distinct; remove repeated entries."))
    letters=sort!(collect(W,selected);by=_letter_location)
    operators=W[one(W)]
    _functional_operators!(operators,letters,one(W),1,Int(order))
    sort!(operators)
    symbols=sort!(unique!(W[first(moment_key(w)) for w in operators if !isone(w)]))
    filter!(!isone,symbols)
    rows=W[]
    _functional_rows!(rows,operators,symbols,one(W),1,Int(order))
    sort!(rows)
end
