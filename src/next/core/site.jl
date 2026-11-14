abstract type SiteRule end
struct Free <: SiteRule end
struct Projector <: SiteRule end
struct Unipotent <: SiteRule end

"Site algebra with operator labels 1:nletters (labels can occur at several sites)."
struct SiteAlgebra{A<:SiteRule}
    nletters::Int
    function SiteAlgebra{A}(nletters::Integer) where {A<:SiteRule}
        nletters > 0 || throw(ArgumentError("A site algebra needs at least one operator label."))
        new{A}(nletters)
    end
end
SiteAlgebra(::A, nletters::Integer) where {A<:SiteRule} = SiteAlgebra{A}(nletters)

"Packed letters (operator << S | site), lowest B bits first; zero is end-of-word."
struct SiteWord{A<:SiteRule,U<:Unsigned,B,S} <: AbstractWord
    bits::U
    function SiteWord{A,U,B,S}(bits::U) where {A<:SiteRule,U<:Unsigned,B,S}
        0 < S < B <= 8sizeof(U) || throw(ArgumentError("Require 0 < site bits < letter bits <= storage bits."))
        iszero(bits >> (B * (8sizeof(U) ÷ B))) || throw(ArgumentError("Unused high storage bits must be zero."))
        mask, smask = (one(U) << B) - one(U), (one(U) << S) - one(U)
        remaining, previous, previous_site = bits, zero(U), zero(U)
        while !iszero(remaining)
            letter = remaining & mask
            site = letter & smask
            site >= previous_site && site > 0 && letter >> S > 0 ||
                throw(ArgumentError("Site normal form needs nonzero labels and increasing site groups."))
            A !== Free && letter == previous && throw(ArgumentError("Adjacent equal letters are not in normal form."))
            previous, previous_site = letter, site
            remaining >>= B
        end
        new{A,U,B,S}(bits)
    end
    # Only internal kernels call this after maintaining the normal-form invariant.
    SiteWord{A,U,B,S}(bits::U, ::Val{:normal}) where {A<:SiteRule,U<:Unsigned,B,S} = new{A,U,B,S}(bits)
end

function word_type(a::SiteAlgebra{A}, nsites::Integer, max_letters::Integer) where {A}
    nsites >= 0 && max_letters >= 0 || throw(ArgumentError("Site and letter counts must be nonnegative."))
    s = max(1, ndigits(nsites; base=2))
    b = s + ndigits(a.nletters; base=2)
    max_letters <= 128 ÷ b || throw(ArgumentError("Site words have a hard maximum of $(128 ÷ b) letters for $b-bit letters (128 storage bits)."))
    needed = max(b, b * max_letters)
    for U in (UInt8, UInt16, UInt32, UInt64, UInt128)
        needed <= 8sizeof(U) && return SiteWord{A,U,b,s}
    end
    throw(ArgumentError("Site letters have a hard maximum of 128 bits."))
end
Base.one(::Type{W}) where {A,U,B,S,W<:SiteWord{A,U,B,S}} = W(zero(U), Val(:normal))
Base.isone(w::SiteWord) = iszero(w.bits)
Base.:(==)(a::W, b::W) where {W<:SiteWord} = a.bits == b.bits
Base.hash(w::SiteWord, h::UInt) = hash(w.bits, h)
Base.isless(a::W, b::W) where {W<:SiteWord} = (degree(a), a.bits) < (degree(b), b.bits)
degree(w::SiteWord{A,U,B}) where {A,U,B} = cld(8sizeof(U) - leading_zeros(w.bits), B)
ProductStyle(::Type{<:SiteWord}) = SingleTerm()
@inline _letter(w::SiteWord{A,U,B}, i::Int) where {A,U,B} = (w.bits >> (B * (i - 1))) & ((one(U) << B) - one(U))
@inline _site(::Type{<:SiteWord{A,U,B,S}}, letter) where {A,U,B,S} = Int(letter & ((one(U) << S) - one(U)))

struct WordSites{W<:SiteWord}
    word::W
end
Base.eltype(::Type{<:WordSites}) = Int
function Base.iterate(s::WordSites{W}, state=(1, 0)) where {W}
    i, previous = state
    n = degree(s.word)
    while i <= n
        site = _site(W, _letter(s.word, i))
        i += 1
        site != previous && return site, (i, site)
    end
    return nothing
end
Base.length(s::WordSites) = count(_ -> true, s)
sites(w::SiteWord) = WordSites(w)

function site_letter(::Type{W}, site::Integer, label::Integer) where {A,U,B,S,W<:SiteWord{A,U,B,S}}
    1 <= site <= (one(U) << S) - 1 || throw(ArgumentError("Site does not fit $S bits."))
    1 <= label <= (one(U) << (B - S)) - 1 || throw(ArgumentError("Operator label does not fit $(B-S) bits."))
    return W((U(label) << S) | U(site), Val(:normal))
end

function mul(a::W, b::W) where {A,U,B,S,W<:SiteWord{A,U,B,S}}
    i = j = 1
    na, nb = degree(a), degree(b)
    bits, n = zero(U), 0
    mask = (one(U) << B) - one(U)
    while i <= na || j <= nb
        if j > nb || (i <= na && _site(W, _letter(a, i)) <= _site(W, _letter(b, j)))
            letter = _letter(a, i); i += 1
        else
            letter = _letter(b, j); j += 1
        end
        previous = n == 0 ? zero(U) : (bits >> (B * (n - 1))) & mask
        if A !== Free && previous == letter
            if A === Unipotent
                n -= 1
                bits &= ~(mask << (B * n))
            end
        else
            n < 8sizeof(U) ÷ B || throw(ArgumentError("Product exceeds the word capacity of $(8sizeof(U) ÷ B) letters; select a wider word_type."))
            bits |= letter << (B * n)
            n += 1
        end
    end
    return one(Phase), W(bits, Val(:normal))
end

function star(w::W) where {A,U,B,S,W<:SiteWord{A,U,B,S}}
    bits, n, start = zero(U), degree(w), 1
    while start <= n
        stop = start
        site = _site(W, _letter(w, start))
        while stop < n && _site(W, _letter(w, stop + 1)) == site
            stop += 1
        end
        for i in start:stop
            bits |= _letter(w, stop - i + start) << (B * (i - 1))
        end
        start = stop + 1
    end
    return one(Phase), W(bits, Val(:normal))
end

# The temporary bridge fixes the old encoding to UInt32: 8 site bits.
function legacy_letters(w::W) where {A,U,B,S,W<:SiteWord{A,U,B,S}}
    out = Int[]
    for i in 1:degree(w)
        letter = _letter(w, i)
        site, label = _site(W, letter), Int(letter >> S)
        site <= 255 && label < 1 << 24 || throw(ArgumentError("Legacy UInt32 encoding supports 255 sites and 2^24-1 labels."))
        push!(out, (label << 8) | site)
    end
    return out
end
function word_from_legacy(::Type{W}, letters::AbstractVector{<:Integer}) where {A,U,B,S,W<:SiteWord{A,U,B,S}}
    length(letters) <= 8sizeof(U) ÷ B || throw(ArgumentError("Legacy word exceeds packed word capacity."))
    bits = zero(U)
    for (i, l) in enumerate(letters)
        0 < l <= typemax(UInt32) || throw(ArgumentError("Expected UInt32 legacy encoding."))
        bits |= site_letter(W, l & 255, l >> 8).bits << (B * (i - 1))
    end
    return W(bits)
end
