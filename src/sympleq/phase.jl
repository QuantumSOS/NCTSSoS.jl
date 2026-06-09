# =============================================================================
# SympleQ Step 7: best-effort phase verification
# =============================================================================

"""
    PhaseVector

`ℤ₄` phase corrections for images of the `2n` binary Pauli generators. This MVP
keeps the vector explicit and marks whether it was verified against the tableau.
"""
struct PhaseVector
    values::Vector{Int8}
    verified::Bool

    function PhaseVector(values::AbstractVector{<:Integer}, verified::Bool)
        return new(Int8[mod(Int(v), 4) for v in values], verified)
    end
end

Base.length(φ::PhaseVector) = length(φ.values)
Base.getindex(φ::PhaseVector, i::Integer) = φ.values[i]

"""
    SympleQGenerator

A candidate Clifford symmetry found by SympleQ's find side.

`phase_verified=false` is not swept under the rug. It means the binary
symplectic action was found, but the companion-paper phase lift was not proven
from the available public specification.
"""
struct SympleQGenerator
    permutation::TermPermutation
    S::SymplecticMatrix
    phase::PhaseVector
    phase_verified::Bool
end

function _phase_dot(row::AbstractVector{<:Integer}, phase::PhaseVector)
    length(row) == length(phase) || throw(DimensionMismatch("Phase vector length does not match symplectic row."))
    acc = 0
    for i in eachindex(row)
        isodd(row[i]) && (acc += Int(phase[i]))
    end
    return mod(acc, 4)
end

function _verify_zero_phase_lift(tab::SymplecticTableau, perm::TermPermutation, S::SymplecticMatrix)
    _gf2_matmul(tab.paulis, S.data) == tab.paulis[perm.images, :] || return false

    # This is intentionally conservative. Public SympleQ text does not pin down
    # Step 7's phase recovery. A zero phase vector is verified only when the
    # canonical Pauli-word phases and coefficients already match the term
    # permutation. More exotic Clifford phases are returned as unverified.
    for i in 1:length(tab)
        j = perm.images[i]
        tab.eta[i] == tab.eta[j] || return false
        tab.coeffs[i] == tab.coeffs[j] || return false
    end
    return true
end

"""
    recover_phase_vector(tab, perm, S)

Best-effort Step-7 phase recovery for the raw generator list. It verifies the
zero phase lift and otherwise returns `phase_verified=false`; the public
`sympleq_symmetry_spec` bridge performs a Hamiltonian-specific GF(2) phase solve
before accepting a Clifford action.
"""
function recover_phase_vector(tab::SymplecticTableau, perm::TermPermutation, S::SymplecticMatrix)
    values = zeros(Int8, size(S.data, 1))
    verified = _verify_zero_phase_lift(tab, perm, S)
    return PhaseVector(values, verified)
end

"""
    sympleq_generators(H; cycle_strategy=:cycle_basis, ga_backend=:backtracking)

Run SympleQ find-side Steps 1--7 and return candidate binary Clifford
generators.
"""
function sympleq_generators(
    H::Polynomial{PauliAlgebra,T,C};
    cycle_strategy::Symbol=:cycle_basis,
    ga_backend::Symbol=:backtracking,
    max_automorphisms::Int=100_000,
) where {T<:Integer,C<:Number}
    tab = SymplecticTableau(H)
    graph = cycle_augmented_graph(tab; cycle_strategy)
    perms = automorphism_generators(graph; backend=ga_backend, max_automorphisms)
    generators = SympleQGenerator[]

    for perm in perms
        S = try
            symplectic_matrix_from_permutation(tab, perm)
        catch err
            @warn "Skipping graph automorphism that is not a symplectic tableau action." exception=(err, catch_backtrace())
            continue
        end
        phase = recover_phase_vector(tab, perm, S)
        push!(generators, SympleQGenerator(perm, S, phase, phase.verified))
    end

    return generators
end
