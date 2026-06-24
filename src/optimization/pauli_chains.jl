"""
    pauli_contiguous_chain_basis(registry, degree; periodic=true)

Build the sparse Pauli half-basis with contiguous support on a one-dimensional
chain.  The basis contains the identity and every Pauli word supported on a
contiguous window of length `1:degree`.  Periodic wrapping is enabled by default.

This is the basis used for large translation-invariant spin-chain relaxations:
for `degree < nqubits` and `periodic=true` its size is
`1 + nqubits * sum(3^k for k in 1:degree)` instead of the full NPA basis size.
"""
function pauli_contiguous_chain_basis(
    registry::VariableRegistry{PauliAlgebra,T},
    degree::Integer;
    periodic::Bool=true,
) where {T<:Integer}
    d = Int(degree)
    d >= 0 || throw(ArgumentError("`degree` must be nonnegative; got $d."))

    nqubits = _pauli_chain_nqubits(registry)
    basis = Set{NormalMonomial{PauliAlgebra,T}}()
    push!(basis, one(NormalMonomial{PauliAlgebra,T}))
    d == 0 && return sort!(collect(basis))

    max_width = min(d, nqubits)
    sizehint!(basis, _pauli_contiguous_chain_basis_size_hint(nqubits, max_width; periodic))

    for width in 1:max_width
        starts = periodic ? (1:nqubits) : (1:(nqubits - width + 1))
        for start in starts
            sites = if periodic
                [mod1(start + offset, nqubits) for offset in 0:(width - 1)]
            else
                collect(start:(start + width - 1))
            end
            for code in 0:(3^width - 1)
                push!(basis, _pauli_chain_word(T, sites, code))
            end
        end
    end

    return sort!(collect(basis))
end

function _pauli_chain_nqubits(registry::VariableRegistry{PauliAlgebra,T}) where {T<:Integer}
    !isempty(registry.idx_to_variables) || throw(ArgumentError("Pauli chain basis needs a non-empty registry."))
    nqubits = maximum(_pauli_site(idx) for idx in keys(registry.idx_to_variables))
    for site in 1:nqubits, pauli_type in (_PAULI_X_TYPE, _PAULI_Y_TYPE, _PAULI_Z_TYPE)
        idx = convert(T, _pauli_index(site, pauli_type))
        haskey(registry.idx_to_variables, idx) || throw(ArgumentError(
            "Pauli chain basis needs a complete site-contiguous Pauli registry; missing index $idx for site $site."
        ))
    end
    length(registry.idx_to_variables) == 3 * nqubits || throw(ArgumentError(
        "Pauli chain basis needs a complete site-contiguous Pauli registry with 3 variables per site."
    ))
    return nqubits
end

function _pauli_contiguous_chain_basis_size_hint(nqubits::Integer, degree::Integer; periodic::Bool)
    if periodic
        return 1 + Int(nqubits) * sum(3^width for width in 1:Int(degree); init=0)
    end
    return 1 + sum((Int(nqubits) - width + 1) * 3^width for width in 1:Int(degree); init=0)
end

function _pauli_chain_word(::Type{T}, sites::AbstractVector{<:Integer}, code::Integer) where {T<:Integer}
    raw = Vector{T}(undef, length(sites))
    value = Int(code)
    @inbounds for i in eachindex(sites)
        pauli_type = value % 3
        value ÷= 3
        raw[i] = convert(T, _pauli_index(sites[i], pauli_type))
    end

    word, phase = simplify(PauliAlgebra, raw)
    phase == UInt8(0) || throw(ArgumentError(
        "Internal Pauli chain basis construction produced a phased word; this should not happen for distinct sites."
    ))
    return NormalMonomial{PauliAlgebra,T}(word)
end

"""
    pauli_sign_symmetry(nqubits; integer_type=Int)

Construct the global Pauli sign symmetry induced by conjugation with
`∏ᵢ σᵢᶻ`: `Xᵢ ↦ -Xᵢ`, `Yᵢ ↦ -Yᵢ`, and `Zᵢ ↦ Zᵢ`.
"""
function pauli_sign_symmetry(nqubits::Integer; integer_type::Type{T}=Int) where {T<:Integer}
    n = _clifford_validate_nqubits(nqubits)
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    sizehint!(images, 2 * n)
    for site in 1:n
        x = _pauli_letter(T, site, _PAULI_X_TYPE)
        y = _pauli_letter(T, site, _PAULI_Y_TYPE)
        images[x] = (-1, x)
        images[y] = (-1, y)
    end
    return CliffordSymmetry(images; nqubits=n)
end
