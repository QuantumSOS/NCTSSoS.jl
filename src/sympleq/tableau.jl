# =============================================================================
# SympleQ Steps 1--2: Pauli polynomial -> binary symplectic tableau
# =============================================================================

"""
    SymplecticTableau

Pauli Hamiltonian tableau used by the find-side SympleQ pipeline.

Rows correspond to nonzero Pauli terms. `paulis[i, :]` is the binary vector
`(x₁,…,xₙ,z₁,…,zₙ)` for term `i`; `eta[i]` is the canonical Pauli-word phase in
`ℤ₄` from writing `Y = iXZ`. Polynomial coefficients are stored separately and
are used as graph colours.
"""
struct SymplecticTableau{T<:Integer,C<:Number}
    coeffs::Vector{C}
    paulis::Matrix{UInt8}
    eta::Vector{Int8}
    monomials::Vector{NormalMonomial{PauliAlgebra,T}}
    nqubits::Int
end

Base.length(tab::SymplecticTableau) = length(tab.monomials)
Base.size(tab::SymplecticTableau) = size(tab.paulis)

function Base.show(io::IO, tab::SymplecticTableau)
    print(io, "SymplecticTableau($(length(tab)) terms, $(tab.nqubits) qubits)")
end

@inline _gf2(x::Integer) = UInt8(isodd(x) ? 1 : 0)

function _pauli_nqubits(poly::Polynomial{PauliAlgebra,T,C}) where {T<:Integer,C<:Number}
    max_site = 0
    for (_, mono) in poly.terms, idx in mono.word
        max_site = max(max_site, _pauli_site(idx))
    end
    return max_site
end

function _pauli_word_to_symplectic(
    mono::NormalMonomial{PauliAlgebra,T},
    nqubits::Integer,
) where {T<:Integer}
    row = zeros(UInt8, 2 * nqubits)
    eta = Int8(0)

    for idx in mono.word
        site = _pauli_site(idx)
        site <= nqubits || throw(ArgumentError(
            "Pauli word contains site $site but tableau has only $nqubits qubits."
        ))
        typ = _pauli_type(idx)
        if typ == 0          # X
            row[site] = 1
        elseif typ == 1      # Y = iXZ
            row[site] = 1
            row[nqubits + site] = 1
            eta = Int8(mod(Int(eta) + 1, 4))
        elseif typ == 2      # Z
            row[nqubits + site] = 1
        else
            throw(ArgumentError("Invalid Pauli type $typ for index $idx."))
        end
    end

    return row, eta
end

function _symplectic_row_to_pauli_word(
    row::AbstractVector{<:Integer},
    ::Type{T}=UInt8,
) where {T<:Integer}
    iseven(length(row)) || throw(ArgumentError("Symplectic rows must have even length."))
    nqubits = length(row) ÷ 2
    word = T[]
    eta = Int8(0)

    for site in 1:nqubits
        x = isodd(row[site])
        z = isodd(row[nqubits + site])
        if x && z
            push!(word, T(_pauli_index(site, 1)))
            eta = Int8(mod(Int(eta) + 1, 4))
        elseif x
            push!(word, T(_pauli_index(site, 0)))
        elseif z
            push!(word, T(_pauli_index(site, 2)))
        end
    end

    return word, eta
end

function _single_pauli_letter_from_symplectic(
    row::AbstractVector{<:Integer},
    ::Type{T}=UInt8,
) where {T<:Integer}
    word, eta = _symplectic_row_to_pauli_word(row, T)
    length(word) == 1 || return nothing
    return only(word), eta
end

"""
    SymplecticTableau(H::Polynomial{PauliAlgebra})

Construct the binary symplectic tableau for a Pauli polynomial. Identity terms
are allowed; they become all-zero rows and are usually irrelevant for graph
automorphisms.
"""
function SymplecticTableau(poly::Polynomial{PauliAlgebra,T,C}) where {T<:Integer,C<:Number}
    nqubits = _pauli_nqubits(poly)
    m = length(poly.terms)
    paulis = zeros(UInt8, m, 2 * nqubits)
    eta = Vector{Int8}(undef, m)
    coeffs = Vector{C}(undef, m)
    monos = Vector{NormalMonomial{PauliAlgebra,T}}(undef, m)

    for (i, (coef, mono)) in enumerate(poly.terms)
        row, row_eta = _pauli_word_to_symplectic(mono, nqubits)
        if !isempty(row)
            paulis[i, :] .= row
        end
        eta[i] = row_eta
        coeffs[i] = coef
        monos[i] = mono
    end

    return SymplecticTableau{T,C}(coeffs, paulis, eta, monos, nqubits)
end
