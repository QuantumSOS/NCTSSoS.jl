"Pauli algebra descriptor. Width depends on the number of sites, not the degree."
struct Pauli end

"Hermitian Pauli strings, with 64 sites per limb and Y = iXZ."
struct PauliWord{K} <: AbstractWord
    x::NTuple{K,UInt64}
    z::NTuple{K,UInt64}
    function PauliWord(x::NTuple{K,UInt64}, z::NTuple{K,UInt64}) where {K}
        K > 0 || throw(ArgumentError("A Pauli word needs at least one limb."))
        new{K}(x, z)
    end
end

function word_type(::Pauli, nsites::Integer, max_letters::Integer)
    nsites >= 0 && max_letters >= 0 || throw(ArgumentError("Site and letter counts must be nonnegative."))
    return PauliWord{max(1, cld(nsites, 64))}
end

Base.one(::Type{PauliWord{K}}) where {K} = PauliWord(ntuple(_ -> UInt64(0), Val(K)), ntuple(_ -> UInt64(0), Val(K)))
Base.isone(w::PauliWord) = all(iszero, w.x) && all(iszero, w.z)
Base.:(==)(a::PauliWord{K}, b::PauliWord{K}) where {K} = a.x == b.x && a.z == b.z
Base.hash(w::PauliWord, h::UInt) = hash(w.z, hash(w.x, h))
Base.isless(a::PauliWord{K}, b::PauliWord{K}) where {K} = (degree(a), a.x, a.z) < (degree(b), b.x, b.z)
degree(w::PauliWord{K}) where {K} = sum(ntuple(i -> count_ones(w.x[i] | w.z[i]), Val(K)))
star(w::PauliWord) = (one(Phase), w)
ProductStyle(::Type{<:PauliWord}) = SingleTerm()

struct MaskSites{K}
    masks::NTuple{K,UInt64}
end
Base.eltype(::Type{<:MaskSites}) = Int
Base.length(s::MaskSites) = sum(count_ones, s.masks)
function Base.iterate(s::MaskSites{K}, state=(1, s.masks[1])) where {K}
    i, m = state
    while iszero(m)
        i += 1
        i > K && return nothing
        m = s.masks[i]
    end
    return (64 * (i - 1) + trailing_zeros(m) + 1, (i, m & (m - UInt64(1))))
end
sites(w::PauliWord{K}) where {K} = MaskSites(ntuple(i -> w.x[i] | w.z[i], Val(K)))

function mul(a::PauliWord{K}, b::PauliWord{K}) where {K}
    x = ntuple(i -> a.x[i] ⊻ b.x[i], Val(K))
    z = ntuple(i -> a.z[i] ⊻ b.z[i], Val(K))
    k = sum(ntuple(i -> count_ones(a.x[i] & a.z[i]) + count_ones(b.x[i] & b.z[i]) +
        2 * count_ones(a.z[i] & b.x[i]) - count_ones(x[i] & z[i]), Val(K)))
    return Phase(k), PauliWord(x, z)
end

"One Pauli letter; axis is 0 (X), 1 (Y), or 2 (Z)."
function pauli_letter(::Type{PauliWord{K}}, site::Integer, axis::Integer) where {K}
    1 <= site <= 64K || throw(ArgumentError("Site must be in 1:$(64K)."))
    0 <= axis <= 2 || throw(ArgumentError("Pauli axis must be 0 (X), 1 (Y), or 2 (Z)."))
    limb, offset = divrem(Int(site) - 1, 64)
    bit = UInt64(1) << offset
    return PauliWord(ntuple(i -> i == limb + 1 && axis != 2 ? bit : UInt64(0), Val(K)),
        ntuple(i -> i == limb + 1 && axis != 0 ? bit : UInt64(0), Val(K)))
end

function legacy_letters(w::PauliWord)
    out = Int[]
    for site in sites(w)
        i, offset = divrem(site - 1, 64)
        bit = UInt64(1) << offset
        x, z = !iszero(w.x[i + 1] & bit), !iszero(w.z[i + 1] & bit)
        push!(out, 3(site - 1) + (x && z ? 2 : x ? 1 : 3))
    end
    return out
end
function word_from_legacy(::Type{W}, letters::AbstractVector{<:Integer}) where {W<:PauliWord}
    out, previous = one(W), 0
    for l in letters
        site, axis = divrem(Int(l) - 1, 3)
        site += 1
        site > previous || throw(ArgumentError("Pauli normal form requires strictly increasing sites."))
        out = mul(out, pauli_letter(W, site, axis))[2]
        previous = site
    end
    return out
end
