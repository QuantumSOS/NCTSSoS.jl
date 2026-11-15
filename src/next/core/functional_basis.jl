function _pack_word(::Type{W},w::FunctionalWord{F,A,B,S}) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    degree(w)<=_fcap(W) || _foverflow(W,degree(w))
    out=one(W)
    for i in 1:degree(w)
        out=_fappend(out,_ftoken(w,i))
        _fend(w,i) && (out=_fmark(out,i))
    end
    _flengths(out,degree(w),_ftail(w))
end

function _functional_operators!(rows,algebra,group,w::W,start,remaining) where {F,A,B,S,K,W<:FunctionalWord{F,A,B,S,K}}
    remaining==0 && return rows
    for j in start:length(group),label in 1:algebra.base.nletters
        v=site_letter(W,group[j],label)
        A!==Free && !isone(w) && _ftoken(w,degree(w))==_ftoken(v,1) && continue
        next=mul(w,v)[2]
        push!(rows,next)
        _functional_operators!(rows,algebra,group,next,j,remaining-1)
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
"Complete mixed rows, graded by all operator letters including those inside symbols."
function basis_rows(::Type{W},algebra::FunctionalAlgebra,group,order::Integer) where {W<:FunctionalWord}
    0<=order<=_fcap(W) || throw(ArgumentError("Row order must be nonnegative and within the functional word capacity $(_fcap(W))."))
    selected=sort!(collect(Int,group))
    allunique(selected) && all(>(0),selected) || throw(ArgumentError("Row sites must be distinct positive integers."))
    operators=W[one(W)]
    _functional_operators!(operators,algebra,selected,one(W),1,Int(order))
    sort!(operators)
    symbols=sort!(unique!(W[first(moment_key(w)) for w in operators if !isone(w)]))
    filter!(!isone,symbols)
    rows=W[]
    _functional_rows!(rows,operators,symbols,one(W),1,Int(order))
    sort!(rows)
end
