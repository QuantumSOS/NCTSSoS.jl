# =============================================================================
# SympleQ Step 6: synthesize binary symplectic S from a term permutation
# =============================================================================

"""
    SymplecticMatrix

Binary row-convention symplectic matrix. A Pauli row vector `p` maps to `p*S`.
"""
struct SymplecticMatrix
    data::Matrix{UInt8}

    function SymplecticMatrix(data::AbstractMatrix{<:Integer})
        size(data, 1) == size(data, 2) || throw(ArgumentError("SymplecticMatrix must be square."))
        iseven(size(data, 1)) || throw(ArgumentError("SymplecticMatrix dimension must be even."))
        mat = UInt8.(mod.(data, 2))
        return new(mat)
    end
end

Base.size(S::SymplecticMatrix) = size(S.data)
Base.getindex(S::SymplecticMatrix, i::Integer, j::Integer) = S.data[i, j]
Base.Matrix(S::SymplecticMatrix) = copy(S.data)

function Base.show(io::IO, S::SymplecticMatrix)
    print(io, "SymplecticMatrix($(size(S.data, 1))×$(size(S.data, 2)))")
end

function _gf2_matmul(A::AbstractMatrix{<:Integer}, B::AbstractMatrix{<:Integer})
    size(A, 2) == size(B, 1) || throw(DimensionMismatch("GF(2) matrix product dimension mismatch."))
    C = zeros(UInt8, size(A, 1), size(B, 2))
    @inbounds for i in axes(A, 1), k in axes(A, 2)
        isodd(A[i, k]) || continue
        for j in axes(B, 2)
            C[i, j] ⊻= UInt8(isodd(B[k, j]))
        end
    end
    return C
end

function _gf2_matvec(row::AbstractVector{<:Integer}, S::AbstractMatrix{<:Integer})
    length(row) == size(S, 1) || throw(DimensionMismatch("GF(2) row/matrix product dimension mismatch."))
    out = zeros(UInt8, size(S, 2))
    @inbounds for k in eachindex(row)
        isodd(row[k]) || continue
        for j in axes(S, 2)
            out[j] ⊻= UInt8(isodd(S[k, j]))
        end
    end
    return out
end

function _gf2_identity(n::Integer)
    I = zeros(UInt8, n, n)
    for i in 1:n
        I[i, i] = 1
    end
    return I
end

function _gf2_inverse(A::AbstractMatrix{<:Integer})
    n, m = size(A)
    n == m || throw(ArgumentError("GF(2) inverse requires a square matrix."))
    B = hcat(UInt8.(mod.(A, 2)), _gf2_identity(n))
    row = 1

    for col in 1:n
        pivot = 0
        for r in row:n
            if B[r, col] == 1
                pivot = r
                break
            end
        end
        pivot == 0 && throw(ArgumentError("Matrix is singular over GF(2)."))

        if pivot != row
            B[row, :], B[pivot, :] = copy(B[pivot, :]), copy(B[row, :])
        end

        for r in 1:n
            if r != row && B[r, col] == 1
                @inbounds for c in col:(2n)
                    B[r, c] ⊻= B[row, c]
                end
            end
        end
        row += 1
    end

    return B[:, (n + 1):(2n)]
end

function symplectic_form(nqubits::Integer)
    Ω = zeros(UInt8, 2 * nqubits, 2 * nqubits)
    for i in 1:nqubits
        Ω[i, nqubits + i] = 1
        Ω[nqubits + i, i] = 1
    end
    return Ω
end

function _gf2_pairing(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer})
    return _symplectic_product_rows(a, b)
end

function _gf2_append_row(rows::AbstractMatrix{<:Integer}, row::AbstractVector{<:Integer})
    return vcat(UInt8.(mod.(rows, 2)), reshape(UInt8.(mod.(row, 2)), 1, :))
end

function _all_gf2_vectors(dim::Integer)
    dim <= 16 || throw(ArgumentError(
        "SympleQ rank-deficient synthesis currently enumerates GF(2)^$dim; refusing dimensions > 16."
    ))
    vectors = Vector{UInt8}[]
    sizehint!(vectors, 2^dim)
    for mask in 0:(UInt(1) << dim) - UInt(1)
        v = zeros(UInt8, dim)
        for bit in 1:dim
            ((mask >> (bit - 1)) & UInt(1)) == UInt(1) && (v[bit] = 1)
        end
        push!(vectors, v)
    end
    return vectors
end

function _standard_gf2_vectors(dim::Integer)
    vectors = Vector{UInt8}[]
    for i in 1:dim
        v = zeros(UInt8, dim)
        v[i] = 1
        push!(vectors, v)
    end
    return vectors
end

function _choose_domain_extension(D::AbstractMatrix{<:Integer})
    dim = size(D, 2)
    for v in _standard_gf2_vectors(dim)
        _gf2_row_in_span(v, D) || return v
    end
    for v in _all_gf2_vectors(dim)
        _gf2_row_in_span(v, D) || return v
    end
    throw(ArgumentError("No independent GF(2) domain extension vector exists."))
end

function _choose_target_extension(
    domain_vector::AbstractVector{<:Integer},
    D::AbstractMatrix{<:Integer},
    R::AbstractMatrix{<:Integer},
)
    dim = size(R, 2)
    for candidate in _all_gf2_vectors(dim)
        _gf2_row_in_span(candidate, R) && continue
        ok = true
        for i in axes(D, 1)
            if _gf2_pairing(domain_vector, view(D, i, :)) != _gf2_pairing(candidate, view(R, i, :))
                ok = false
                break
            end
        end
        ok && return candidate
    end
    throw(ArgumentError(
        "Failed to extend the partial Pauli isometry to a full symplectic basis. " *
        "The graph automorphism may not preserve the tableau's linear relations."
    ))
end

function _extend_symplectic_isometry(
    domain_rows::AbstractMatrix{<:Integer},
    target_rows::AbstractMatrix{<:Integer},
)
    size(domain_rows) == size(target_rows) || throw(DimensionMismatch("Domain and target bases must have the same shape."))
    dim = size(domain_rows, 2)
    D = UInt8.(mod.(domain_rows, 2))
    R = UInt8.(mod.(target_rows, 2))

    _gf2_rank(D) == size(D, 1) || throw(ArgumentError("Domain rows are not independent over GF(2)."))
    _gf2_rank(R) == size(R, 1) || throw(ArgumentError("Target rows are not independent over GF(2)."))

    while size(D, 1) < dim
        d = _choose_domain_extension(D)
        r = _choose_target_extension(d, D, R)
        D = _gf2_append_row(D, d)
        R = _gf2_append_row(R, r)
    end

    return D, R
end

"""
    is_symplectic_matrix(S::SymplecticMatrix)

Check `S Ω Sᵀ = Ω` over GF(2), matching the row-vector convention used by this
subsystem.
"""
function is_symplectic_matrix(S::SymplecticMatrix)
    dim = size(S.data, 1)
    Ω = symplectic_form(dim ÷ 2)
    return _gf2_matmul(_gf2_matmul(S.data, Ω), transpose(S.data)) == Ω
end

"""
    symplectic_matrix_from_permutation(tab, perm)

Synthesize a row-convention binary symplectic matrix `S` satisfying
`tab.paulis * S == tab.paulis[perm.images, :]` over GF(2) whenever the supplied
term permutation induces a Pauli isometry.

Rank-deficient tableaux are extended by a small finite-field Witt-extension
search. That is the explicit, testable patch for the Step-6 corner case.
"""
function symplectic_matrix_from_permutation(tab::SymplecticTableau, perm::TermPermutation)
    length(perm) == length(tab) || throw(ArgumentError("Permutation length does not match tableau term count."))
    dim = 2 * tab.nqubits
    dim == 0 && return SymplecticMatrix(zeros(UInt8, 0, 0))

    P = UInt8.(tab.paulis)
    Q = P[perm.images, :]
    basis_indices = _gf2_independent_row_indices(P)

    if isempty(basis_indices)
        return SymplecticMatrix(_gf2_identity(dim))
    end

    D0 = P[basis_indices, :]
    R0 = Q[basis_indices, :]
    D, R = _extend_symplectic_isometry(D0, R0)
    S = SymplecticMatrix(_gf2_matmul(_gf2_inverse(D), R))

    is_symplectic_matrix(S) || throw(ArgumentError("Synthesized matrix is not symplectic over GF(2)."))
    _gf2_matmul(P, S.data) == Q || throw(ArgumentError(
        "Synthesized symplectic matrix does not reproduce the requested term permutation."
    ))

    return S
end
