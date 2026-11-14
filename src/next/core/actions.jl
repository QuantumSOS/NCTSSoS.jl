"A finite site/mode permutation; sites beyond its support are fixed."
struct SitePermutation <: AbstractAction
    images::Vector{Int}
    function SitePermutation(images)
        p = collect(Int, images)
        sort(p) == collect(1:length(p)) || throw(ArgumentError("Site images must be a permutation of 1:n."))
        while !isempty(p) && p[end] == length(p)
            pop!(p)
        end
        new(p)
    end
end
@inline _image(g::SitePermutation, i::Int) = i <= length(g.images) ? g.images[i] : i
Base.one(::SitePermutation) = SitePermutation(Int[])
Base.isone(g::SitePermutation) = isempty(g.images)
Base.:(==)(g::SitePermutation,h::SitePermutation) = g.images == h.images
Base.hash(g::SitePermutation,h::UInt) = hash(g.images,h)
compose(g::SitePermutation,h::SitePermutation) = SitePermutation([_image(g,_image(h,i)) for i in 1:max(length(g.images),length(h.images))])
Base.inv(g::SitePermutation) = SitePermutation(invperm(g.images))

function act(g::SitePermutation,w::PauliWord{K}) where {K}
    out = one(typeof(w))
    for site in sites(w)
        limb, offset = divrem(site-1,64)
        bit = UInt64(1) << offset
        x,z = !iszero(w.x[limb+1]&bit), !iszero(w.z[limb+1]&bit)
        out = mul(out,pauli_letter(typeof(w),_image(g,site),x && z ? 1 : x ? 0 : 2))[2]
    end
    return one(Phase),out
end
function _permute_mask(g::SitePermutation,m::NTuple{K,UInt64}) where {K}
    out = ntuple(_->UInt64(0),Val(K))
    parity = 0
    for site in MaskSites(m)
        target = _image(g,site)
        target <= 64K || throw(ArgumentError("Permutation image exceeds the word's mode capacity."))
        limb,offset = divrem(target-1,64)
        bit = UInt64(1) << offset
        parity += _count_above(out,limb+1,bit)
        out = _replace(out,limb+1,out[limb+1]|bit)
    end
    return parity,out
end
function act(g::SitePermutation,w::FermiWord)
    pc,c = _permute_mask(g,w.c)
    pa,a = _permute_mask(g,w.a)
    return Phase(2(pc+pa)),FermiWord(c,a)
end
function act(g::SitePermutation,w::W) where {A,U,B,S,W<:SiteWord{A,U,B,S}}
    bits = zero(U)
    for i in 1:degree(w)
        l = _letter(w,i)
        site = _image(g,_site(W,l))
        mapped = site_letter(W,site,Int(l >> S)).bits
        j = 1
        while j < i && _site(W,_letter(W(bits,Val(:normal)),j)) <= site
            j += 1
        end
        shift = B*(j-1)
        low = (one(U)<<shift)-one(U)
        bits = (bits&low) | ((bits&~low)<<B) | (mapped<<shift)
    end
    return one(Phase),W(bits,Val(:normal))
end

@inline function _symplectic(a::PauliWord{K},b::PauliWord{K}) where {K}
    return isodd(sum(ntuple(i->count_ones(a.x[i]&b.z[i])+count_ones(a.z[i]&b.x[i]),Val(K))))
end
@inline _maskbit(m,site) = !iszero(m[(site-1)÷64+1] & (UInt64(1)<<((site-1)%64)))

"Canonical Pauli action: symplectic X/Z images, sign bits and a conjugation flag. Trailing fixed sites are omitted."
struct CliffordMap{W<:PauliWord} <: AbstractAction
    columns::Vector{W}
    negative::BitVector
    antiunitary::Bool
    function CliffordMap(columns::Vector{W},negative::AbstractVector{Bool}=falses(length(columns));antiunitary::Bool=false) where {W<:PauliWord}
        iseven(length(columns)) && length(columns)==length(negative) || throw(ArgumentError("Supply 2n Pauli columns and 2n sign bits."))
        n = length(columns)÷2
        for j in eachindex(columns)
            all(s->s<=n,sites(columns[j])) || throw(ArgumentError("Clifford columns must act within sites 1:n."))
            for i in 1:j
                _symplectic(columns[i],columns[j]) == (j==i+n) || throw(ArgumentError("Columns do not preserve the Pauli commutation relations."))
            end
        end
        active = n
        while active>0 && !negative[active] && !negative[n+active] &&
                columns[active]==pauli_letter(W,active,0) && columns[n+active]==pauli_letter(W,active,2)
            active-=1
        end
        keep = vcat(collect(1:active),collect(n+1:n+active))
        new{W}(columns[keep],BitVector(negative[keep]),antiunitary)
    end
end
is_antiunitary(g::CliffordMap) = g.antiunitary
function _identity_columns(::Type{W},n::Integer) where {W<:PauliWord}
    n >= 0 || throw(ArgumentError("Site count must be nonnegative."))
    return [pauli_letter(W,i,axis) for axis in (0,2) for i in 1:n]
end
Base.one(g::CliffordMap{W}) where {W} = CliffordMap(_identity_columns(W,length(g.columns)÷2))
function Base.isone(g::CliffordMap{W}) where {W}
    return isempty(g.columns) && !g.antiunitary
end
Base.:(==)(g::CliffordMap,h::CliffordMap) = g.columns==h.columns && g.negative==h.negative && g.antiunitary==h.antiunitary
Base.hash(g::CliffordMap,h::UInt) = hash(g.antiunitary,hash(g.negative,hash(g.columns,h)))
function act(g::CliffordMap{W},w::W) where {W}
    n = length(g.columns)÷2
    ny = sum(count_ones,w.x .& w.z)
    phase = Phase(g.antiunitary ? -ny : ny)
    out = one(W)
    for which in 0:1
        mask = which==0 ? w.x : w.z
        for site in MaskSites(mask)
            if site <= n
                j = site+which*n
                p,out = mul(out,g.columns[j])
                phase = phase*p*Phase(2Int(g.negative[j]))
            else
                p,out = mul(out,pauli_letter(W,site,2which))
                phase = phase*p
            end
        end
    end
    return phase,out
end
function compose(g::CliffordMap{W},h::CliffordMap{W}) where {W}
    n = max(length(g.columns),length(h.columns))÷2
    columns = _identity_columns(W,n)
    negative = falses(2n)
    for i in eachindex(columns)
        p,w = act(h,columns[i])
        q,columns[i] = act(g,w)
        negative[i] = (p*q).k==2
    end
    return CliffordMap(columns,negative;antiunitary=xor(g.antiunitary,h.antiunitary))
end
function Base.inv(g::CliffordMap{W}) where {W}
    n = length(g.columns)÷2
    columns = fill(one(W),2n)
    negative = falses(2n)
    for j in 1:2n
        site,which = mod1(j,n),j<=n ? 0 : 1
        w = one(W)
        for i in 1:n
            x = _maskbit(which==0 ? g.columns[n+i].z : g.columns[n+i].x,site)
            z = _maskbit(which==0 ? g.columns[i].z : g.columns[i].x,site)
            if x || z
                w = mul(w,pauli_letter(W,i,x && z ? 1 : x ? 0 : 2))[2]
            end
        end
        columns[j] = w
        negative[j] = act(g,w)[1].k==2
    end
    return CliffordMap(columns,negative;antiunitary=g.antiunitary)
end
function clifford_h(::Type{W},n::Integer,site::Integer) where {W<:PauliWord}
    1<=site<=n || throw(ArgumentError("H site must be in 1:n."))
    columns = _identity_columns(W,n)
    columns[site],columns[n+site] = columns[n+site],columns[site]
    return CliffordMap(columns)
end
function clifford_s(::Type{W},n::Integer,site::Integer) where {W<:PauliWord}
    1<=site<=n || throw(ArgumentError("S site must be in 1:n."))
    columns = _identity_columns(W,n)
    columns[site] = pauli_letter(W,site,1)
    return CliffordMap(columns)
end
function clifford_cnot(::Type{W},n::Integer,control::Integer,target::Integer) where {W<:PauliWord}
    1<=control<=n && 1<=target<=n && control!=target || throw(ArgumentError("CNOT needs two different sites in 1:n."))
    columns = _identity_columns(W,n)
    columns[control] = mul(columns[control],columns[target])[2]
    columns[n+target] = mul(columns[n+control],columns[n+target])[2]
    return CliffordMap(columns)
end
function clifford_swap(::Type{W},n::Integer,a::Integer,b::Integer) where {W<:PauliWord}
    1<=a<=n && 1<=b<=n && a!=b || throw(ArgumentError("SWAP needs two different sites in 1:n."))
    columns = _identity_columns(W,n)
    columns[a],columns[b] = columns[b],columns[a]
    columns[n+a],columns[n+b] = columns[n+b],columns[n+a]
    return CliffordMap(columns)
end

"Conjugation by a Pauli word, optionally after complex conjugation K."
struct PauliConjugation{W<:PauliWord} <: AbstractAction
    word::W
    antiunitary::Bool
end
AntiunitaryPauli(w::PauliWord) = PauliConjugation(w,true)
is_antiunitary(g::PauliConjugation) = g.antiunitary
Base.one(g::PauliConjugation{W}) where {W} = PauliConjugation(one(W),false)
Base.isone(g::PauliConjugation) = isone(g.word) && !g.antiunitary
Base.:(==)(g::PauliConjugation,h::PauliConjugation) = g.word==h.word && g.antiunitary==h.antiunitary
Base.hash(g::PauliConjugation,h::UInt) = hash(g.antiunitary,hash(g.word,h))
Base.inv(g::PauliConjugation) = g
compose(g::PauliConjugation{W},h::PauliConjugation{W}) where {W} = compose(CliffordMap(g),CliffordMap(h))
function act(g::PauliConjugation{W},w::W) where {W}
    sign = _symplectic(g.word,w)
    g.antiunitary && (sign=xor(sign,isodd(sum(count_ones,w.x .& w.z))))
    return Phase(2Int(sign)),w
end

"Convert a site permutation to the canonical action on Pauli words of type W."
function CliffordMap(::Type{W},n::Integer,g::SitePermutation) where {W<:PauliWord}
    n>=length(g.images) || throw(ArgumentError("Permutation acts beyond site $n; increase the conversion site count."))
    columns = _identity_columns(W,n)
    for i in eachindex(columns); columns[i]=act(g,columns[i])[2]; end
    return CliffordMap(columns)
end
CliffordMap(g::SitePermutation,::Type{W},n::Integer) where {W<:PauliWord} = CliffordMap(W,n,g)
function CliffordMap(g::PauliConjugation{W}) where {W}
    n = maximum(sites(g.word);init=0)
    columns = _identity_columns(W,n)
    negative = BitVector([act(g,w)[1].k==2 for w in columns])
    return CliffordMap(columns,negative;antiunitary=g.antiunitary)
end
_canonical_pauli(g::CliffordMap,::Type) = g
_canonical_pauli(g::PauliConjugation,::Type) = CliffordMap(g)
_canonical_pauli(g::SitePermutation,::Type{W}) where {W<:PauliWord} = CliffordMap(W,length(g.images),g)
function compose(g::Union{CliffordMap{W},PauliConjugation{W}},h::SitePermutation) where {W}
    return compose(_canonical_pauli(g,W),_canonical_pauli(h,W))
end
function compose(g::SitePermutation,h::Union{CliffordMap{W},PauliConjugation{W}}) where {W}
    return compose(_canonical_pauli(g,W),_canonical_pauli(h,W))
end
compose(g::CliffordMap{W},h::PauliConjugation{W}) where {W} = compose(g,CliffordMap(h))
compose(g::PauliConjugation{W},h::CliffordMap{W}) where {W} = compose(CliffordMap(g),h)
compose(g::AbstractAction,h::AbstractAction,k::AbstractAction,rest::AbstractAction...) =
    foldl(compose,rest;init=compose(compose(g,h),k))

"Infinitesimal spin rotation d(w) = -i/2 [Σ σ_axis, w] on the selected sites."
struct SpinRotation{K} <: AbstractDerivation
    axis::Int
    support::NTuple{K,UInt64}
    function SpinRotation(axis::Integer,support::NTuple{K,UInt64}) where {K}
        0<=axis<=2 || throw(ArgumentError("Spin axis must be 0 (X), 1 (Y), or 2 (Z)."))
        new{K}(axis,support)
    end
end
Base.:(==)(a::SpinRotation,b::SpinRotation) = a.axis==b.axis && a.support==b.support
Base.hash(d::SpinRotation,h::UInt) = hash(d.axis,hash(d.support,h))
function derive!(sink,d::SpinRotation{K},w::PauliWord{K},c::Number=1) where {K}
    for site in MaskSites(ntuple(i->d.support[i]&(w.x[i]|w.z[i]),Val(K)))
        generator = pauli_letter(typeof(w),site,d.axis)
        if _symplectic(generator,w)
            p,image = mul(generator,w)
            push_term!(sink,scale(Phase(3)*p,c),image)
        end
    end
    return sink
end
