# =============================================================================
# SympleQ generator -> current CliffordSymmetry / SymmetrySpec seam
# =============================================================================

function _sympleq_check_dimensions(
    S::SymplecticMatrix,
    phase::PhaseVector,
    nqubits::Integer,
)
    size(S.data, 1) == 2 * nqubits || throw(DimensionMismatch(
        "Symplectic matrix dimension $(size(S.data, 1)) does not match $nqubits Pauli qubits."
    ))
    length(phase) == 2 * nqubits || throw(DimensionMismatch(
        "Phase vector length $(length(phase)) does not match $nqubits Pauli qubits."
    ))
    return nothing
end

function _sympleq_real_phase_sign(phase_k::Integer; context::AbstractString="SympleQ Clifford image")
    k = mod(Int(phase_k), 4)
    k == 0 && return 1
    k == 2 && return -1
    throw(ArgumentError(
        "$context produced phase i^$k, which cannot be represented as a Hermitian Pauli image with sign ±1."
    ))
end

function _sympleq_generator_row(nqubits::Integer, site::Integer, pauli_type::Integer)
    row = zeros(UInt8, 2 * nqubits)
    if pauli_type == _PAULI_X_TYPE
        row[site] = 1
    elseif pauli_type == _PAULI_Z_TYPE
        row[nqubits + site] = 1
    else
        throw(ArgumentError("Expected X/Z generator type, got Pauli type $pauli_type."))
    end
    return row
end

function _sympleq_generator_image(
    S::SymplecticMatrix,
    phase::PhaseVector,
    ::Type{T},
    nqubits::Integer,
    site::Integer,
    pauli_type::Integer,
) where {T<:Integer}
    row = _sympleq_generator_row(nqubits, site, pauli_type)
    image_row = _gf2_matvec(row, S.data)
    image_word, image_eta = _symplectic_row_to_pauli_word(image_row, T)
    isempty(image_word) && throw(ArgumentError(
        "SympleQ Clifford image of a non-identity Pauli generator became identity."
    ))

    # Xᵢ/Zᵢ have source η = 0. `phase` is the phase of the binary generator
    # image; converting from binary X/Z convention back to canonical Pauli words
    # subtracts the target word's η (Y = iXZ contributes η = 1).
    sign = _sympleq_real_phase_sign(
        _phase_dot(row, phase) - Int(image_eta);
        context="SympleQ image of Pauli generator on site $site",
    )
    return sign, NormalMonomial{PauliAlgebra,T}(image_word)
end

function _sympleq_y_image_from_xz(
    x_image::Tuple{Int,NormalMonomial{PauliAlgebra,T}},
    z_image::Tuple{Int,NormalMonomial{PauliAlgebra,T}},
) where {T<:Integer}
    x_sign, x_mono = x_image
    z_sign, z_mono = z_image
    product_word, product_phase = simplify(PauliAlgebra, vcat(x_mono.word, z_mono.word))

    # Y = i X Z. The Pauli simplifier returns X'Z' = i^product_phase P'.
    total_phase = Int(product_phase) + 1 + (x_sign * z_sign == -1 ? 2 : 0)
    y_sign = _sympleq_real_phase_sign(total_phase; context="SympleQ derived image of Pauli Y")
    return y_sign, NormalMonomial{PauliAlgebra,T}(product_word)
end

"""
    sympleq_clifford_symmetry(S, φ; integer_type=Int)
    sympleq_clifford_symmetry(generator; integer_type=Int)

Convert a SympleQ binary symplectic generator and phase vector into the
[`CliffordSymmetry`](@ref) action used by `NCTSSoS`'s Pauli symmetry pipeline.

The conversion keeps the honest contract: if the phase lift would map a
Hermitian Pauli generator to an `±im` multiple of a Pauli word, it errors rather
than smuggling a non-Clifford action into the SDP reduction.
"""
function sympleq_clifford_symmetry(
    S::SymplecticMatrix,
    phase::PhaseVector;
    integer_type::Type{T}=Int,
) where {T<:Integer}
    nqubits = size(S.data, 1) ÷ 2
    return sympleq_clifford_symmetry(S, phase, T, nqubits)
end

function sympleq_clifford_symmetry(
    generator::SympleQGenerator;
    integer_type::Type{T}=Int,
) where {T<:Integer}
    return sympleq_clifford_symmetry(generator.S, generator.phase; integer_type=T)
end

function sympleq_clifford_symmetry(
    S::SymplecticMatrix,
    phase::PhaseVector,
    ::Type{T},
    nqubits::Integer,
) where {T<:Integer}
    _sympleq_check_dimensions(S, phase, nqubits)

    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    for site in 1:nqubits
        source_x = _pauli_letter(T, site, _PAULI_X_TYPE)
        source_y = _pauli_letter(T, site, _PAULI_Y_TYPE)
        source_z = _pauli_letter(T, site, _PAULI_Z_TYPE)

        x_image = _sympleq_generator_image(S, phase, T, nqubits, site, _PAULI_X_TYPE)
        z_image = _sympleq_generator_image(S, phase, T, nqubits, site, _PAULI_Z_TYPE)
        y_image = _sympleq_y_image_from_xz(x_image, z_image)

        for (source, (sign, image)) in ((source_x, x_image), (source_y, y_image), (source_z, z_image))
            if sign != 1 || image != source
                images[source] = (sign, image)
            end
        end
    end

    return CliffordSymmetry(images; nqubits)
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

function _sympleq_base_phase_values(S::SymplecticMatrix, nqubits::Integer)
    values = zeros(Int8, 2 * nqubits)
    for site in 1:nqubits
        x_row = _sympleq_generator_row(nqubits, site, _PAULI_X_TYPE)
        z_row = _sympleq_generator_row(nqubits, site, _PAULI_Z_TYPE)
        _, x_eta = _symplectic_row_to_pauli_word(_gf2_matvec(x_row, S.data), Int)
        _, z_eta = _symplectic_row_to_pauli_word(_gf2_matvec(z_row, S.data), Int)
        values[site] = x_eta
        values[nqubits + site] = z_eta
    end
    return values
end

function _sympleq_gf2_solution_space(A::AbstractMatrix{<:Integer}, b::AbstractVector{<:Integer})
    size(A, 1) == length(b) || throw(DimensionMismatch("GF(2) solve has incompatible right-hand side."))
    nvars = size(A, 2)
    if isempty(b)
        return zeros(UInt8, nvars), _gf2_nullspace_basis(zeros(UInt8, 0, nvars))
    end

    augmented = hcat(UInt8.(mod.(A, 2)), reshape(UInt8.(mod.(b, 2)), :, 1))
    rref, pivots = _gf2_rref(augmented)
    rhs_col = nvars + 1

    for r in axes(rref, 1)
        if all(rref[r, c] == 0 for c in 1:nvars) && rref[r, rhs_col] == 1
            throw(ArgumentError("SympleQ phase equations are inconsistent for this graph automorphism."))
        end
    end

    x = zeros(UInt8, nvars)
    for (r, pivot_col) in enumerate(pivots)
        pivot_col <= nvars || continue
        x[pivot_col] = rref[r, rhs_col]
    end
    return x, _gf2_nullspace_basis(A)
end

function _sympleq_sign_bit(sign::Integer)
    sign == 1 && return UInt8(0)
    sign == -1 && return UInt8(1)
    throw(ArgumentError("Expected sign ±1, got $sign."))
end

function _sympleq_coeff_ratio_sign(source, target)
    target == source && return 1
    target == -source && return -1
    throw(ArgumentError(
        "SympleQ graph automorphism maps coefficient $(repr(source)) to incompatible coefficient $(repr(target))."
    ))
end

function _sympleq_phase_from_sign_bits(base_phase::PhaseVector, sign_bits::AbstractVector{<:Integer})
    values = copy(base_phase.values)
    for i in eachindex(values)
        values[i] = Int8(mod(Int(values[i]) + 2 * Int(sign_bits[i]), 4))
    end
    return PhaseVector(values, true)
end

function _sympleq_active_phase_indices(tab::SymplecticTableau)
    sites = Set{Int}()
    for mono in tab.monomials, idx in mono.word
        push!(sites, _pauli_site(idx))
    end
    active = Int[]
    for site in sort!(collect(sites))
        push!(active, site)
        push!(active, tab.nqubits + site)
    end
    return active
end

function _sympleq_phase_vectors_for_hamiltonian(
    tab::SymplecticTableau{T,C},
    generator::SympleQGenerator,
) where {T<:Integer,C<:Number}
    nqubits = tab.nqubits
    _sympleq_check_dimensions(generator.S, generator.phase, nqubits)

    base_phase = PhaseVector(_sympleq_base_phase_values(generator.S, nqubits), true)
    base_clifford = sympleq_clifford_symmetry(generator.S, base_phase, T, nqubits)

    rhs = UInt8[]
    sizehint!(rhs, length(tab))
    for i in 1:length(tab)
        j = generator.permutation.images[i]
        sign0, image0 = _act_monomial(base_clifford, tab.monomials[i])
        image0 == tab.monomials[j] || throw(ArgumentError(
            "SympleQ phase solve expected term $i to map to term $j, got $image0 instead of $(tab.monomials[j])."
        ))

        desired_sign = _sympleq_coeff_ratio_sign(tab.coeffs[i], tab.coeffs[j])
        push!(rhs, _sympleq_sign_bit(sign0 * desired_sign))
    end

    particular, kernel = _sympleq_gf2_solution_space(tab.paulis, rhs)
    phases = PhaseVector[_sympleq_phase_from_sign_bits(base_phase, particular)]
    active_phase_indices = _sympleq_active_phase_indices(tab)
    for direction in kernel
        any(idx -> direction[idx] == 1, active_phase_indices) || continue
        push!(phases, _sympleq_phase_from_sign_bits(base_phase, particular .⊻ direction))
    end
    return phases
end

function _sympleq_clifford_key(g::CliffordSymmetry{T}, nqubits::Integer) where {T<:Integer}
    return _symmetry_key(g, _clifford_letter_domain(T, nqubits))
end

function _sympleq_clifford_is_identity(g::CliffordSymmetry{T}, nqubits::Integer) where {T<:Integer}
    domain = _clifford_letter_domain(T, nqubits)
    return all(idx -> _clifford_action_signature(g, idx) == (1, (idx,)), domain)
end

"""
    sympleq_symmetry_spec(H; cycle_strategy=:cycle_basis, ga_backend=:backtracking, max_automorphisms=100_000)

Recognize Clifford symmetries of a Pauli Hamiltonian using SympleQ's find-side
pipeline and return a ready-to-use [`SymmetrySpec`](@ref).

The implemented pipeline is the qubit/Pauli version of SympleQ: build a binary
symplectic tableau, encode commutation and linear-dependency structure as a
coloured graph, enumerate colour-preserving automorphisms, synthesize candidate
symplectic Clifford actions, convert them to `CliffordSymmetry`, and keep only
candidates that leave `H` invariant under the existing Pauli action.

This is a recognizer/bridge, not the SympleQ exploit-side qubit-tapering
machinery. Dense/no-sparsity SDP reductions are handled by the existing
`SolverConfig(symmetry=...)` path.
"""
function sympleq_symmetry_spec(
    H::Polynomial{PauliAlgebra,T,C};
    cycle_strategy::Symbol=:cycle_basis,
    ga_backend::Symbol=:backtracking,
    max_automorphisms::Int=100_000,
) where {T<:Integer,C<:Number}
    tab = SymplecticTableau(H)
    tab.nqubits > 0 || throw(ArgumentError("SympleQ needs a non-constant Pauli Hamiltonian."))

    generators = sympleq_generators(
        H;
        cycle_strategy,
        ga_backend,
        max_automorphisms,
    )
    identity_perm = TermPermutation(collect(1:length(tab)))
    identity_S = SymplecticMatrix(_gf2_identity(2 * tab.nqubits))
    push!(
        generators,
        SympleQGenerator(identity_perm, identity_S, PhaseVector(zeros(Int8, 2 * tab.nqubits), true), true),
    )

    cliffords = CliffordSymmetry{T}[]
    seen = Set{Any}()

    for generator in generators
        phases = try
            _sympleq_phase_vectors_for_hamiltonian(tab, generator)
        catch err
            @warn "Skipping SympleQ generator whose phase lift is not representable as a CliffordSymmetry." exception=(err, catch_backtrace())
            continue
        end

        for phase in phases
            clifford = try
                sympleq_clifford_symmetry(generator.S, phase, T, tab.nqubits)
            catch err
                @warn "Skipping SympleQ phase solution that does not define a valid CliffordSymmetry." exception=(err, catch_backtrace())
                continue
            end

            _sympleq_clifford_is_identity(clifford, tab.nqubits) && continue
            _valid_clifford_symmetry_for_polynomial(clifford, H) || continue

            key = _sympleq_clifford_key(clifford, tab.nqubits)
            key in seen && continue
            push!(seen, key)
            push!(cliffords, clifford)
        end
    end

    isempty(cliffords) && throw(ArgumentError(
        "SympleQ did not find any nontrivial Clifford symmetry that preserves this Pauli Hamiltonian."
    ))

    return SymmetrySpec(clifford_generators=cliffords; check_invariance=true)
end
