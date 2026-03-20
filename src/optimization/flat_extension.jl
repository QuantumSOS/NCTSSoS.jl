"""
    FlatnessResult

Summary of a numerical flatness check for a Hankel matrix partitioned into a
principal block and its extension block.
"""
struct FlatnessResult
    is_flat::Bool
    rank_principal::Int
    rank_full::Int
    err_flat::Float64
end

function _flat_extension_projection(
    Htilde::AbstractMatrix,
    B::AbstractMatrix;
    atol::Real=1e-8,
)
    factor = svd(Matrix(Htilde))
    rank_principal, cutoff = _gns_rank(factor.S; atol=atol, rtol=0.0)
    rank_principal == 0 && throw(
        ArgumentError(
            "Principal Hankel block has numerical rank 0 at cutoff $cutoff; cannot build a flat extension.",
        ),
    )

    U = factor.U[:, 1:rank_principal]
    Vt = factor.Vt[1:rank_principal, :]
    S_inv = Diagonal(inv.(factor.S[1:rank_principal]))
    Z = Vt' * (S_inv * (U' * Matrix(B)))
    projected_B = Matrix(Htilde) * Z
    C_flat = Z' * projected_B
    B_residual = norm(Matrix(B) - projected_B) / (1 + norm(B) + norm(projected_B))

    return (; Z, C_flat, B_residual, rank_principal)
end

function _flat_extension_partition(
    hankel::AbstractMatrix{T},
    full_basis::Vector{M},
    basis::Vector{M},
) where {T<:Number,M}
    size(hankel, 1) == size(hankel, 2) ||
        throw(ArgumentError("Hankel matrix must be square, got size $(size(hankel))."))
    length(full_basis) == size(hankel, 1) || throw(
        ArgumentError(
            "Hankel matrix size $(size(hankel, 1)) does not match full basis length $(length(full_basis)).",
        ),
    )

    principal = _gns_basis_indices(full_basis, basis)
    is_principal = falses(length(full_basis))
    is_principal[principal] .= true
    extension = findall(!, is_principal)
    perm = vcat(principal, extension)
    permuted = Matrix(hankel[perm, perm])

    n_principal = length(principal)
    Htilde = permuted[1:n_principal, 1:n_principal]
    B = permuted[1:n_principal, n_principal + 1:end]
    C = permuted[n_principal + 1:end, n_principal + 1:end]

    return (; principal, extension, perm, n_principal, Htilde, B, C)
end

"""
    test_flatness(hankel, full_basis, basis; atol=1e-8)

Check whether `hankel` is numerically flat over the principal block indexed by
`basis`.
"""
function test_flatness(
    hankel::AbstractMatrix,
    full_basis::Vector,
    basis::Vector;
    atol::Real=1e-8,
)
    part = _flat_extension_partition(hankel, full_basis, basis)

    principal_svals = svdvals(Matrix(part.Htilde))
    rank_principal, _ = _gns_rank(principal_svals; atol=atol, rtol=0.0)

    full_svals = svdvals(Matrix(hankel))
    rank_full, _ = _gns_rank(full_svals; atol=atol, rtol=0.0)

    if isempty(part.extension)
        return FlatnessResult(rank_full == rank_principal, rank_principal, rank_full, 0.0)
    end

    projection = _flat_extension_projection(part.Htilde, part.B; atol=atol)
    C_flat = projection.C_flat
    err_flat = norm(part.C - C_flat) / (1 + norm(part.C) + norm(C_flat))

    return FlatnessResult(
        rank_full == rank_principal && projection.B_residual < atol && err_flat < atol,
        rank_principal,
        rank_full,
        Float64(err_flat),
    )
end

"""
    flat_extend(hankel, full_basis, basis; atol=1e-8)

Replace the extension block of `hankel` by its canonical flat projection
`Z' * Htilde * Z`, where `Htilde * Z = B`.
"""
function flat_extend(
    hankel::AbstractMatrix{T},
    full_basis::Vector{M},
    basis::Vector{M};
    atol::Real=1e-8,
) where {T<:Number,M}
    part = _flat_extension_partition(hankel, full_basis, basis)
    isempty(part.extension) && return Matrix(hankel)

    projection = _flat_extension_projection(part.Htilde, part.B; atol=atol)
    C_flat = projection.C_flat
    flat_permuted = Matrix(hankel[part.perm, part.perm])
    flat_permuted[part.n_principal + 1:end, part.n_principal + 1:end] .= C_flat

    inv_perm = invperm(part.perm)
    return flat_permuted[inv_perm, inv_perm]
end
