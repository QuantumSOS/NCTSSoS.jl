# =============================================================================
# Phase 6: Clifford/SympleQ generator -> existing SymmetrySpec seam
# =============================================================================

function _pauli_registry_nqubits(registry::VariableRegistry{PauliAlgebra,T}) where {T<:Integer}
    isempty(registry.idx_to_variables) && return 0
    max_idx = maximum(Int.(keys(registry.idx_to_variables)))
    max_idx % 3 == 0 || throw(ArgumentError(
        "Pauli registry should contain complete X/Y/Z triples; max index is $max_idx."
    ))
    return max_idx ÷ 3
end

function _validate_pauli_letter_registry(registry::VariableRegistry{PauliAlgebra,T}) where {T<:Integer}
    nqubits = _pauli_registry_nqubits(registry)
    for idx in one(T):T(3 * nqubits)
        haskey(registry.idx_to_variables, idx) || throw(ArgumentError(
            "Pauli registry is not the standard letter-level layout: missing index $idx."
        ))
    end
    return nqubits
end

function _phase_sign_for_row(row::AbstractVector{<:Integer}, phase::PhaseVector)
    k = _phase_dot(row, phase)
    k == 0 && return 1
    k == 2 && return -1
    throw(ArgumentError(
        "Clifford phase i^$k cannot be represented by SignedPermutation signs ±1."
    ))
end

function _clifford_to_signed_permutation(
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

    images = Dict{T,Tuple{Int,T}}()
    for site in 1:nqubits, typ in 0:2
        src = T(_pauli_index(site, typ))
        row = zeros(UInt8, 2 * nqubits)
        if typ == 0          # X
            row[site] = 1
        elseif typ == 1      # Y = X + Z at the binary level
            row[site] = 1
            row[nqubits + site] = 1
        else                 # Z
            row[nqubits + site] = 1
        end

        image_row = _gf2_matvec(row, S.data)
        image = _single_pauli_letter_from_symplectic(image_row, T)
        isnothing(image) && throw(ArgumentError(
            "Clifford image of Pauli letter index $src is a multi-qubit Pauli word. " *
            "The existing SymmetrySpec/SignedPermutation seam is letter-level; " *
            "word-level Clifford actions are deliberately not faked."
        ))
        dst, _eta = image
        sign = _phase_sign_for_row(row, phase)
        (sign != 1 || dst != src) && (images[src] = (sign, dst))
    end

    return SignedPermutation(images)
end

"""
    clifford_to_signed_permutation(S, φ, registry)

Convert a supported Clifford action to the existing letter-level
`SignedPermutation` representation for a standard Pauli `VariableRegistry`.

This intentionally rejects entangling Clifford images such as `X₁ ↦ X₁X₂`.
Those are real Clifford symmetries, but the current `SymmetrySpec` seam cannot
represent them as letter permutations. Lying here would corrupt the relaxation.
"""
function clifford_to_signed_permutation(
    S::SymplecticMatrix,
    phase::PhaseVector,
    registry::VariableRegistry{PauliAlgebra,T},
) where {T<:Integer}
    nqubits = _validate_pauli_letter_registry(registry)
    return _clifford_to_signed_permutation(S, phase, T, nqubits)
end

function clifford_to_signed_permutation(
    generator::SympleQGenerator,
    registry::VariableRegistry{PauliAlgebra,T},
) where {T<:Integer}
    return clifford_to_signed_permutation(generator.S, generator.phase, registry)
end

function _signed_permutation_key(g::SignedPermutation{T}, domain::AbstractVector{T}) where {T<:Integer}
    return Tuple(_signed_image(g, idx) for idx in domain)
end

function _infer_pauli_registry(::Type{T}, nqubits::Integer) where {T<:Integer}
    idx_to_vars = Dict{T,Symbol}()
    vars_to_idx = Dict{Symbol,T}()
    for site in 1:nqubits
        for (typ, prefix) in enumerate(("σx", "σy", "σz"))
            idx = T(_pauli_index(site, typ - 1))
            sym = Symbol(prefix * _subscript_string(site))
            idx_to_vars[idx] = sym
            vars_to_idx[sym] = idx
        end
    end
    return VariableRegistry{PauliAlgebra,T}(idx_to_vars, vars_to_idx)
end

function _valid_signed_symmetry_for_polynomial(g::SignedPermutation{T}, H::Polynomial{PauliAlgebra,T,C}) where {T<:Integer,C<:Number}
    try
        return _act_polynomial(g, H) == H
    catch
        return false
    end
end

"""
    sympleq_symmetry_spec(H; cycle_strategy=:cycle_basis, ga_backend=:bliss)

Discover supported SympleQ Pauli symmetries and return a ready-to-feed
`SymmetrySpec` for dense/no-sparsity relaxations.

Only letter-level Clifford actions are converted today because that is what the
existing `SignedPermutation` pipeline can represent for Pauli registries.
"""
function sympleq_symmetry_spec(
    H::Polynomial{PauliAlgebra,T,C};
    cycle_strategy::Symbol=:cycle_basis,
    ga_backend::Symbol=:bliss,
    max_automorphisms::Int=100_000,
) where {T<:Integer,C<:Number}
    tab = SymplecticTableau(H)
    registry = _infer_pauli_registry(T, tab.nqubits)
    generators = sympleq_generators(
        H;
        cycle_strategy,
        ga_backend,
        max_automorphisms,
    )

    signed = SignedPermutation{T}[]
    seen = Set{Any}()
    domain = T.(1:(3 * tab.nqubits))

    for generator in generators
        !generator.phase_verified && @warn(
            "Including a SympleQ generator with phase_verified=false; downstream invariance checks remain enabled."
        )
        sp = try
            clifford_to_signed_permutation(generator, registry)
        catch err
            @warn "Skipping SympleQ generator not representable as a Pauli letter SignedPermutation." exception=(err, catch_backtrace())
            continue
        end
        isone_on_domain = all(idx -> _signed_image(sp, idx) == (1, idx), domain)
        isone_on_domain && continue
        if !_valid_signed_symmetry_for_polynomial(sp, H)
            @warn "Skipping SympleQ generator whose letter-level action does not leave the Hamiltonian invariant."
            continue
        end
        key = _signed_permutation_key(sp, domain)
        key in seen && continue
        push!(seen, key)
        push!(signed, sp)
    end

    isempty(signed) && throw(ArgumentError(
        "SympleQ did not find any nontrivial symmetry representable by the current Pauli letter-level SignedPermutation seam."
    ))

    return SymmetrySpec(signed; check_invariance=true)
end
