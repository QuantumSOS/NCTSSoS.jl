# A small word type for the tests of the contracts and of the reference
# builders. It is not the package word type: it holds at most 64 sites and is
# written for clarity. The core layer supplies the real `PauliWord`.

"Pauli string on at most 64 sites: `prod(im^(x*z) X^x Z^z)` over the sites. X = (1,0), Y = (1,1), Z = (0,1)."
struct ToyPauli <: Nx.AbstractWord
    x::UInt64
    z::UInt64
end

"One Pauli letter. `axis` is 0 for X, 1 for Y and 2 for Z, as in the current package."
function toy_letter(site::Integer, axis::Integer)
    bit = UInt64(1) << (site - 1)
    return ToyPauli(axis != 2 ? bit : UInt64(0), axis != 0 ? bit : UInt64(0))
end

Base.one(::Type{ToyPauli}) = ToyPauli(0, 0)
Base.isone(w::ToyPauli) = (w.x | w.z) == 0
Base.hash(w::ToyPauli, h::UInt) = hash(w.z, hash(w.x, h))
Base.isless(a::ToyPauli, b::ToyPauli) = (Nx.degree(a), a.x, a.z) < (Nx.degree(b), b.x, b.z)

Nx.degree(w::ToyPauli) = count_ones(w.x | w.z)

struct BitSites
    mask::UInt64
end
Base.eltype(::Type{BitSites}) = Int
Base.length(s::BitSites) = count_ones(s.mask)
Base.iterate(s::BitSites, m::UInt64=s.mask) = m == 0 ? nothing : (trailing_zeros(m) + 1, m & (m - 1))
Nx.sites(w::ToyPauli) = BitSites(w.x | w.z)

Nx.star(w::ToyPauli) = (one(Nx.Phase), w)
Nx.ProductStyle(::Type{ToyPauli}) = Nx.SingleTerm()
function Nx.mul(a::ToyPauli, b::ToyPauli)
    x, z = a.x ⊻ b.x, a.z ⊻ b.z
    k = count_ones(a.x & a.z) + count_ones(b.x & b.z) + 2 * count_ones(a.z & b.x) - count_ones(x & z)
    return Nx.Phase(k), ToyPauli(x, z)
end

function Nx.legacy_letters(w::ToyPauli)
    letters = Int[]
    for site in Nx.sites(w)
        bit = UInt64(1) << (site - 1)
        x, z = w.x & bit != 0, w.z & bit != 0
        push!(letters, 3 * (site - 1) + (x && z ? 1 : x ? 0 : 2) + 1)
    end
    return letters
end
function Nx.word_from_legacy(::Type{ToyPauli}, letters::AbstractVector{<:Integer})
    w = one(ToyPauli)
    for l in letters
        w = Nx.mul(w, toy_letter((Int(l) - 1) ÷ 3 + 1, (Int(l) - 1) % 3))[2]   # one letter for each site: no phase
    end
    return w
end

"All toy words with at most `d` letters on sites `1:n`, the identity first."
function toy_words(n::Integer, d::Integer)
    out = [one(ToyPauli)]
    function extend(current, start, remaining)
        for site in start:n, axis in 0:2
            next = Nx.mul(current, toy_letter(site, axis))[2]
            push!(out, next)
            remaining > 1 && extend(next, site + 1, remaining - 1)
        end
    end
    extend(one(ToyPauli), 1, d)
    return out
end

"Cyclic shift of the sites `1:n` by one site."
struct ToyShift <: Nx.AbstractAction
    n::Int
end
function Nx.act(g::ToyShift, w::ToyPauli)
    mask = (UInt64(1) << g.n) - 1
    rotate(m) = ((m << 1) | (m >> (g.n - 1))) & mask
    return one(Nx.Phase), ToyPauli(rotate(w.x), rotate(w.z))
end
