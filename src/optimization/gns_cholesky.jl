function _gns_cholesky_basis_indices(
    hankel_block::AbstractMatrix{T};
    atol::Real=1e-8,
    rtol::Real=0.0,
) where {T<:Number}
    qr_factor = qr(Matrix(hankel_block), ColumnNorm())
    principal_svals = svdvals(Matrix(hankel_block))
    rank_basis, cutoff = _gns_rank(principal_svals; atol=atol, rtol=rtol)
    rank_basis == 0 && throw(
        ArgumentError(
            "No singular values exceed the numerical cutoff $cutoff. Consider decreasing `atol`/`rtol`.",
        ),
    )

    pivots = collect(qr_factor.p[1:rank_basis])
    if 1 ∉ pivots
        pivots = vcat(1, filter(!=(1), pivots))[1:rank_basis]
    end

    selected_block = Matrix(hankel_block[pivots, pivots])
    selected_svals = svdvals(selected_block)
    minimum(selected_svals) > cutoff || throw(
        ArgumentError(
            "Failed to select a numerically independent GNS basis containing the identity word.",
        ),
    )

    return pivots, principal_svals
end

function _gns_finalize_cholesky(
    hankel::AbstractMatrix{T},
    full_basis::Vector{M},
    basis::Vector{M},
    registry::VariableRegistry{A,TI};
    atol::Real=1e-8,
    rtol::Real=0.0,
) where {T<:Number,A<:AlgebraType,TI<:Integer,M<:NormalMonomial{A,TI}}
    size(hankel, 1) == size(hankel, 2) ||
        throw(ArgumentError("Hankel matrix must be square, got size $(size(hankel))."))
    length(full_basis) == size(hankel, 1) || throw(
        ArgumentError(
            "Hankel matrix size $(size(hankel, 1)) does not match full basis length $(length(full_basis)).",
        ),
    )

    hankel_indices = _gns_basis_indices(full_basis, basis)
    hankel_block = Matrix(hankel[hankel_indices, hankel_indices])
    selected, principal_svals = _gns_cholesky_basis_indices(
        hankel_block;
        atol=atol,
        rtol=rtol,
    )

    rank_basis = length(selected)
    full_svals = svdvals(Matrix(hankel))
    rank_full, _ = _gns_rank(full_svals; atol=atol, rtol=rtol)

    if rank_full != rank_basis
        @warn "Flatness condition violated: rank(H) = $rank_full ≠ rank(hankel_block) = $rank_basis. The full Hankel matrix is not a flat extension of the principal block."
    end

    gns_basis = basis[selected]
    columns = hankel_block[:, selected]
    principal_subblock = Matrix(hankel_block[selected, selected])
    chol = cholesky(Hermitian((principal_subblock + principal_subblock') / 2))
    G = Matrix(chol.U)
    Ginv = inv(G)
    Trec = promote_type(eltype(hankel), eltype(G), Float64)

    matrices = Dict{TI,Matrix{Trec}}()
    for var_idx in indices(registry)
        K = _gns_localizing_from_hankel(hankel, full_basis, var_idx, basis)
        localizing_columns = Matrix(K[:, selected])
        Abar = columns \ localizing_columns
        X = G * Abar * Ginv
        matrices[var_idx] = _gns_postprocess_operator(A, X)
    end

    e1 = zeros(Trec, rank_basis)
    e1[1] = one(Trec)
    xi = G * e1

    return GNSResult{Trec,eltype(principal_svals),A,TI,M}(
        matrices,
        Vector{Trec}(xi),
        collect(gns_basis),
        full_basis,
        collect(principal_svals),
        rank_basis,
        rank_full,
    )
end
