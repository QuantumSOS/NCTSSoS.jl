"Binary SympleQ tableau with owning packed Pauli words, in detector term order."
struct SymplecticTableau{W<:PauliWord,C<:Number}
    coeffs::Vector{C}
    paulis::Matrix{UInt8}
    eta::Vector{Int8}
    monomials::Vector{W}
    nqubits::Int
end
Base.length(tab::SymplecticTableau) = length(tab.monomials)

@inline _gf2(x::Integer) = UInt8(isodd(x))

function _symplectic_row_to_pauli_word(row::AbstractVector{<:Integer}, ::Type{W}) where {W<:PauliWord}
    iseven(length(row)) || throw(ArgumentError("Symplectic rows must have even length."))
    n = length(row) ÷ 2
    word = one(W)
    eta = 0
    for i in 1:n
        x, z = isodd(row[i]), isodd(row[n+i])
        (x || z) || continue
        word = mul(word, pauli_letter(W,i,x ? (z ? 1 : 0) : 2))[2]
        eta += x && z
    end
    return word, Int8(mod(eta,4))
end

function SymplecticTableau(poly::Polynomial{W}) where {W<:PauliWord}
    # The search's tie breaks use graded lexicographic site/axis order, not
    # packed-mask order. Keep that order without constructing legacy words.
    ordered = sort!(collect(terms(poly)); by=t->(degree(t[2]),collect(_locations(t[2]))))
    n = maximum((maximum(sites(w);init=0) for (_,w) in ordered);init=0)
    paulis = zeros(UInt8,length(ordered),2n)
    eta = zeros(Int8,length(ordered))
    for (i,(_,w)) in enumerate(ordered)
        for (site,axis) in _locations(w)
            paulis[i,site] = axis != 2
            paulis[i,n+site] = axis != 0
            eta[i] = mod(eta[i] + (axis == 1),4)
        end
    end
    return SymplecticTableau(ComplexF64[c for (c,_) in ordered],paulis,eta,
        W[w for (_,w) in ordered],n)
end
