struct Fermionic end

"Normal-ordered creators (descending modes), then annihilators (ascending modes)."
struct FermiWord{K} <: AbstractWord
    c::NTuple{K,UInt64}
    a::NTuple{K,UInt64}
    function FermiWord(c::NTuple{K,UInt64}, a::NTuple{K,UInt64}) where {K}
        K > 0 || throw(ArgumentError("A fermionic word needs at least one limb."))
        new{K}(c, a)
    end
end
function word_type(::Fermionic, nsites::Integer, max_letters::Integer)
    nsites >= 0 && max_letters >= 0 || throw(ArgumentError("Mode and letter counts must be nonnegative."))
    return FermiWord{max(1, cld(nsites, 64))}
end
Base.one(::Type{FermiWord{K}}) where {K} = FermiWord(ntuple(_->UInt64(0),Val(K)),ntuple(_->UInt64(0),Val(K)))
Base.isone(w::FermiWord) = all(iszero,w.c) && all(iszero,w.a)
Base.:(==)(a::FermiWord{K}, b::FermiWord{K}) where {K} = a.c == b.c && a.a == b.a
Base.hash(w::FermiWord,h::UInt) = hash(w.a,hash(w.c,h))
Base.isless(a::FermiWord{K},b::FermiWord{K}) where {K} = (degree(a),a.c,a.a) < (degree(b),b.c,b.a)
degree(w::FermiWord) = sum(count_ones,w.c) + sum(count_ones,w.a)
sites(w::FermiWord{K}) where {K} = MaskSites(ntuple(i->w.c[i] | w.a[i],Val(K)))
star(w::FermiWord) = (one(Phase),FermiWord(w.a,w.c))
ProductStyle(::Type{<:FermiWord}) = MultiTerm()

@inline function _count_above(m::NTuple{K,UInt64}, limb::Int, bit::UInt64) where {K}
    return sum(ntuple(i -> i > limb ? count_ones(m[i]) : i == limb ? count_ones(m[i] & ~(bit | (bit - UInt64(1)))) : 0,Val(K)))
end
@inline function _count_below(m::NTuple{K,UInt64}, limb::Int, bit::UInt64) where {K}
    return sum(ntuple(i -> i < limb ? count_ones(m[i]) : i == limb ? count_ones(m[i] & (bit - UInt64(1))) : 0,Val(K)))
end
@inline _replace(m::NTuple{K,UInt64}, limb::Int, value::UInt64) where {K} = ntuple(i->i==limb ? value : m[i],Val(K))

function _mul!(::MultiTerm, sink, l::FermiWord{K}, r::FermiWord{K}, coef::Number) where {K}
    limb = K
    while limb > 0 && iszero(r.c[limb])
        limb -= 1
    end
    if limb > 0
        bit = UInt64(1) << (63 - leading_zeros(r.c[limb]))
        rest = FermiWord(_replace(r.c,limb,r.c[limb] ⊻ bit),r.a)
        if !iszero(l.a[limb] & bit)
            sign = isodd(_count_above(l.a,limb,bit)) ? -coef : coef
            _mul!(MultiTerm(),sink,FermiWord(l.c,_replace(l.a,limb,l.a[limb] ⊻ bit)),rest,sign)
        end
        if iszero(l.c[limb] & bit)
            sign = isodd(sum(count_ones,l.a) + _count_below(l.c,limb,bit)) ? -coef : coef
            _mul!(MultiTerm(),sink,FermiWord(_replace(l.c,limb,l.c[limb] | bit),l.a),rest,sign)
        end
    elseif all(iszero,ntuple(i->l.a[i] & r.a[i],Val(K)))
        parity = 0
        for j in 1:K
            remaining = r.a[j]
            while !iszero(remaining)
                bit = remaining & -remaining
                parity += _count_above(l.a,j,bit)
                remaining &= remaining - UInt64(1)
            end
        end
        w = FermiWord(l.c,ntuple(i->l.a[i] | r.a[i],Val(K)))
        push_term!(sink,isodd(parity) ? -coef : coef,w)
    end
    return sink
end

"One fermionic operator; negative legacy index creates, positive index annihilates."
function fermi_letter(::Type{FermiWord{K}}, index::Integer) where {K}
    -64K <= index <= 64K && index != 0 || throw(ArgumentError("Fermionic index must be nonzero with absolute value at most $(64K)."))
    limb, offset = divrem(abs(Int(index))-1,64)
    bit = UInt64(1) << offset
    masks = ntuple(i->i==limb+1 ? bit : UInt64(0),Val(K))
    empty = ntuple(_->UInt64(0),Val(K))
    return index < 0 ? FermiWord(masks,empty) : FermiWord(empty,masks)
end
function legacy_letters(w::FermiWord)
    out = Int[]
    for site in MaskSites(w.c)
        push!(out,-site)
    end
    reverse!(out)
    for site in MaskSites(w.a)
        push!(out,site)
    end
    return out
end
function word_from_legacy(::Type{W}, letters::AbstractVector{<:Integer}) where {K,W<:FermiWord{K}}
    c, a, previous = one(W).c, one(W).a, typemin(Int)
    for l in letters
        l > previous || throw(ArgumentError("Fermionic normal form requires strictly increasing signed indices."))
        letter = fermi_letter(W,l)
        c = ntuple(i->c[i] | letter.c[i],Val(K))
        a = ntuple(i->a[i] | letter.a[i],Val(K))
        previous = Int(l)
    end
    return FermiWord(c,a)
end
