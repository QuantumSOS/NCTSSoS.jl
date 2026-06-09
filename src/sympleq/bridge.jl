# =============================================================================
# SympleQ generator -> CliffordSymmetry bridge
# =============================================================================

function _pauli_letter_symplectic_row(
    site::Integer,
    pauli_type::Integer,
    nqubits::Integer,
)
    row = zeros(UInt8, 2 * nqubits)
    eta = Int8(0)

    if pauli_type == 0          # X
        row[site] = 1
    elseif pauli_type == 1      # Y = iXZ
        row[site] = 1
        row[nqubits + site] = 1
        eta = Int8(1)
    elseif pauli_type == 2      # Z
        row[nqubits + site] = 1
    else
        throw(ArgumentError("Invalid Pauli type $pauli_type."))
    end

    return row, eta
end

function _clifford_image_sign(
    source_eta::Integer,
    source_row::AbstractVector{<:Integer},
    target_eta::Integer,
    phase::PhaseVector,
)
    k = mod(Int(source_eta) + _phase_dot(source_row, phase) - Int(target_eta), 4)
    k == 0 && return 1
    k == 2 && return -1
    throw(ArgumentError(
        "Clifford phase i^$k cannot be represented as a Hermitian Pauli-word sign ±1."
    ))
end

function _clifford_to_clifford_symmetry(
    S::SymplecticMatrix,
    phase::PhaseVector,
    ::Type{T},
    nqubits::Integer,
) where {T<:Integer}
    size(S.data, 1) == 2 * nqubits || throw(DimensionMismatch(
        "Symplectic matrix dimension $(size(S.data, 1)) does not match $nqubits Pauli qubits."
    ))
    length(phase) == 2 * nqubits || throw(DimensionMismatch(
        "Phase vector length $(length(phase)) does not match $nqubits Pauli qubits."
    ))

    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    for site in 1:nqubits, pauli_type in 0:2
        source_idx = convert(T, _pauli_index(site, pauli_type))
        source = NormalMonomial{PauliAlgebra,T}(T[source_idx])
        source_row, source_eta = _pauli_letter_symplectic_row(site, pauli_type, nqubits)

        target_row = _gf2_matvec(source_row, S.data)
        target_word, target_eta = _symplectic_row_to_pauli_word(target_row, T)
        target = NormalMonomial{PauliAlgebra,T}(target_word)
        sign = _clifford_image_sign(source_eta, source_row, target_eta, phase)

        if sign != 1 || target != source
            images[source] = (sign, target)
        end
    end

    return CliffordSymmetry(images; nqubits)
end

"""
    clifford_to_clifford_symmetry(S, φ; nqubits, integer_type=Int)
    clifford_to_clifford_symmetry(generator; nqubits, integer_type=Int)

Convert a SympleQ binary Clifford candidate into NCTSSoS's native
[`CliffordSymmetry`](@ref) action.

The conversion preserves multi-qubit Pauli-word images, unlike the old
letter-permutation bridge. It still rejects non-Hermitian phase lifts because
`CliffordSymmetry` represents conjugation of Hermitian Pauli words by a real
sign and a Pauli word.
"""
function clifford_to_clifford_symmetry(
    S::SymplecticMatrix,
    phase::PhaseVector;
    nqubits::Integer,
    integer_type::Type{T}=Int,
) where {T<:Integer}
    return _clifford_to_clifford_symmetry(S, phase, T, nqubits)
end

function clifford_to_clifford_symmetry(
    generator::SympleQGenerator;
    nqubits::Integer,
    integer_type::Type{T}=Int,
) where {T<:Integer}
    return clifford_to_clifford_symmetry(generator.S, generator.phase; nqubits, integer_type=T)
end

function _is_identity_on_domain(g::CliffordSymmetry{T}, domain::AbstractVector{T}) where {T<:Integer}
    return all(idx -> _clifford_action_signature(g, idx) == (1, (idx,)), domain)
end

function _valid_clifford_symmetry_for_polynomial(
    g::CliffordSymmetry,
    H::Polynomial{PauliAlgebra,T,C},
) where {T<:Integer,C<:Number}
    try
        return _act_polynomial(g, H) == H
    catch
        return false
    end
end

"""
    sympleq_symmetry_spec(H; cycle_strategy=:cycle_basis, ga_backend=:backtracking)

Discover SympleQ find-side Pauli symmetries and return a `SymmetrySpec` using
NCTSSoS's native [`CliffordSymmetry`](@ref) representation.

This implements the reproducible find side only: tableau construction,
anticommutation/cycle graph construction, colour-preserving automorphism
enumeration, binary symplectic synthesis, conservative phase recovery, and
conversion to `CliffordSymmetry`. Candidate generators are still checked against
`H`; unsupported or non-invariant phase lifts are skipped rather than faked.
"""
function sympleq_symmetry_spec(
    H::Polynomial{PauliAlgebra,T,C};
    cycle_strategy::Symbol=:cycle_basis,
    ga_backend::Symbol=:backtracking,
    max_automorphisms::Int=100_000,
) where {T<:Integer,C<:Number}
    tab = SymplecticTableau(H)
    generators = sympleq_generators(
        H;
        cycle_strategy,
        ga_backend,
        max_automorphisms,
    )

    clifford = CliffordSymmetry{T}[]
    seen = Set{Any}()
    domain = T.(1:(3 * tab.nqubits))

    for generator in generators
        !generator.phase_verified && @warn(
            "Including a SympleQ generator with phase_verified=false; downstream invariance checks remain enabled."
        )

        g = try
            clifford_to_clifford_symmetry(generator; nqubits=tab.nqubits, integer_type=T)
        catch err
            @warn "Skipping SympleQ generator not representable as a Hermitian CliffordSymmetry." exception=(err, catch_backtrace())
            continue
        end

        _is_identity_on_domain(g, domain) && continue
        if !_valid_clifford_symmetry_for_polynomial(g, H)
            @warn "Skipping SympleQ generator whose Clifford action does not leave the Hamiltonian invariant."
            continue
        end

        key = _symmetry_key(g, domain)
        key in seen && continue
        push!(seen, key)
        push!(clifford, g)
    end

    isempty(clifford) && throw(ArgumentError(
        "SympleQ did not find any nontrivial symmetry representable by the current CliffordSymmetry seam."
    ))

    return SymmetrySpec(clifford; check_invariance=true)
end
