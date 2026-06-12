# =============================================================================
# Ground-state and local physical moment constraints
# =============================================================================

"""
Provenance tag for physical PSD blocks injected outside the automatic
moment/localizing-matrix construction.
"""
struct PhysicalPSDOrigin <: BlockOrigin
    family::Symbol
    block_key::Tuple
end

function _assert_physical_psd_matrix_symmetry(cone::Symbol, matrix::AbstractMatrix)
    if cone == :HPSD
        for i in axes(matrix, 1)
            iszero(matrix[i, i] - adjoint(matrix[i, i])) || throw(ArgumentError(
                "Physical HPSD block has a non-Hermitian diagonal entry at ($i, $i)."
            ))
        end
    end
    for j in axes(matrix, 2), i in first(axes(matrix, 1)):(j - 1)
        diff = cone == :HPSD ? matrix[i, j] - adjoint(matrix[j, i]) : matrix[i, j] - matrix[j, i]
        iszero(diff) || throw(ArgumentError(
            "Physical $(cone) block is not $(cone == :HPSD ? "Hermitian" : "symmetric") at entries ($i, $j) and ($j, $i)."
        ))
    end
    return nothing
end

"""
    PhysicalPSDConstraint(cone, matrix; row_labels=nothing, origin=PhysicalPSDOrigin(:user, ()))

Polynomial-valued PSD/HPSD block to append to a moment relaxation. Entries are
interpreted as moment linear forms `ℓ(matrix[i,j])`; no localizing multipliers
are added.
"""
struct PhysicalPSDConstraint{P<:Polynomial,M<:NormalMonomial}
    cone::Symbol
    matrix::Matrix{P}
    row_labels::Vector{M}
    origin::BlockOrigin

    function PhysicalPSDConstraint{P,M}(
        cone::Symbol,
        matrix::Matrix{P},
        row_labels::Vector{M},
        origin::BlockOrigin,
    ) where {P<:Polynomial,M<:NormalMonomial}
        cone in (:PSD, :HPSD) || throw(ArgumentError(
            "Physical PSD constraint cone must be :PSD or :HPSD, got $(repr(cone))"
        ))
        size(matrix, 1) == size(matrix, 2) || throw(DimensionMismatch(
            "Physical PSD constraint matrix must be square, got $(size(matrix))"
        ))
        _assert_physical_psd_matrix_symmetry(cone, matrix)
        length(row_labels) == size(matrix, 1) || throw(DimensionMismatch(
            "Physical PSD row label count $(length(row_labels)) does not match block size $(size(matrix, 1))"
        ))
        return new{P,M}(cone, matrix, row_labels, origin)
    end
end

function PhysicalPSDConstraint(
    cone::Symbol,
    matrix::Matrix{P};
    row_labels=nothing,
    origin::BlockOrigin=PhysicalPSDOrigin(:user, ()),
) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C}}
    M = NormalMonomial{A,T}
    labels = row_labels === nothing ? [one(M) for _ in 1:size(matrix, 1)] : M[collect(row_labels)...]
    return PhysicalPSDConstraint{P,M}(cone, matrix, labels, origin)
end

Base.size(block::PhysicalPSDConstraint) = size(block.matrix)
Base.size(block::PhysicalPSDConstraint, dim::Integer) = size(block.matrix, dim)

@inline _physical_psd_cone(::Type{A}) where {A<:AlgebraType} = _is_complex_problem(A) ? :HPSD : :PSD

"""
    commutator_constraints(H, basis) -> Vector{Polynomial}

Return scalar moment equalities for the eigenstate condition `ℓ([H,u]) = 0`.
Pass the result to `polyopt(...; scalar_moment_eq_constraints=...)`, not to
`moment_eq_constraints`.
"""
function commutator_constraints(
    H::Polynomial{A,T,C},
    basis::AbstractVector{M},
) where {A<:AlgebraType,T<:Integer,C<:Number,M<:NormalMonomial{A,T}}
    Cout = promote_type(C, coeff_type(A))
    Pout = Polynomial{A,T,Cout}
    constraints = Pout[]
    for u in basis
        comm = convert(Pout, H * u - u * H)
        iszero(comm) || push!(constraints, comm)
    end
    return constraints
end

"""
    curvature_block(H, basis) -> PhysicalPSDConstraint

Build the ground-state optimality/curvature PSD block with entries

`ℓ(v† H w - 1/2 * (H v† w + v† w H))`.

This condition is valid for ground-state minimization problems. It is deliberately
not enabled automatically; solvers can be sensitive to this block's conditioning.
"""
function curvature_block(
    H::Polynomial{A,T,C},
    basis::AbstractVector{M},
) where {A<:AlgebraType,T<:Integer,C<:Number,M<:NormalMonomial{A,T}}
    _is_hermitian_polynomial(H) || throw(ArgumentError("curvature_block requires a Hermitian Hamiltonian"))

    Cout = promote_type(C, coeff_type(A), Float64)
    Pout = Polynomial{A,T,Cout}
    n = length(basis)
    entries = Matrix{Pout}(undef, n, n)
    half = convert(Cout, 0.5)

    for (i, v) in enumerate(basis), (j, w) in enumerate(basis)
        vstar = adjoint(v)
        vstar_w = vstar * w
        entries[i, j] = convert(Pout, vstar * H * w - half * (H * vstar_w + vstar_w * H))
    end

    return PhysicalPSDConstraint(
        _physical_psd_cone(A),
        entries;
        row_labels=collect(basis),
        origin=PhysicalPSDOrigin(:ground_state_curvature, (n,)),
    )
end

function _check_pauli_rdm_sites(registry::VariableRegistry{PauliAlgebra,T}, sites) where {T<:Unsigned}
    isempty(sites) && throw(ArgumentError("rdm_block requires at least one site"))
    length(unique(sites)) == length(sites) || throw(ArgumentError("rdm_block sites must be unique"))

    for site in sites
        site >= 1 || throw(ArgumentError("Pauli site indices are 1-based; got $site"))
        for pauli_type in 0:2
            idx = convert(T, _pauli_index(site, pauli_type))
            haskey(registry.idx_to_variables, idx) || throw(ArgumentError(
                "registry does not contain Pauli operator type $pauli_type at internal site $site"
            ))
        end
    end
    return nothing
end

@inline function _pauli_local_matrix_entry(op::Int, row_bit::Int, col_bit::Int)
    if op == 0       # I
        return row_bit == col_bit ? 1.0 + 0.0im : 0.0 + 0.0im
    elseif op == 1   # X
        return row_bit != col_bit ? 1.0 + 0.0im : 0.0 + 0.0im
    elseif op == 2   # Y
        row_bit == 0 && col_bit == 1 && return 0.0 - 1.0im
        row_bit == 1 && col_bit == 0 && return 0.0 + 1.0im
        return 0.0 + 0.0im
    elseif op == 3   # Z
        return row_bit == col_bit ? (row_bit == 0 ? 1.0 + 0.0im : -1.0 + 0.0im) : 0.0 + 0.0im
    else
        throw(ArgumentError("invalid Pauli code $op"))
    end
end

function _pauli_tensor_matrix_entry(ops::AbstractVector{Int}, row::Int, col::Int)
    k = length(ops)
    coef = 1.0 + 0.0im
    for pos in 1:k
        shift = k - pos
        row_bit = ((row - 1) >> shift) & 0x01
        col_bit = ((col - 1) >> shift) & 0x01
        local_coef = _pauli_local_matrix_entry(ops[pos], Int(row_bit), Int(col_bit))
        iszero(local_coef) && return 0.0 + 0.0im
        coef *= local_coef
    end
    return coef
end

function _pauli_ops_from_code(code::Integer, k::Int)
    ops = Vector{Int}(undef, k)
    value = Int(code)
    for pos in k:-1:1
        ops[pos] = value % 4
        value ÷= 4
    end
    return ops
end

function _pauli_cluster_word(
    ::Type{T},
    sites::AbstractVector{<:Integer},
    ops::AbstractVector{Int},
) where {T<:Unsigned}
    P = Polynomial{PauliAlgebra,T,ComplexF64}
    raw = T[]
    for (site, op) in zip(sites, ops)
        op == 0 && continue
        push!(raw, convert(T, _pauli_index(site, op - 1)))
    end
    isempty(raw) && return one(P)
    return convert(P, Polynomial(simplify(PauliAlgebra, raw)))
end

"""
    rdm_block(registry::VariableRegistry{PauliAlgebra}, sites) -> PhysicalPSDConstraint

Build the spin `k`-RDM positivity block for the given Pauli sites:

`ρ_S(ℓ) = 2^{-k} Σ_a ℓ(σ_{s1}^{a1}⋯σ_{sk}^{ak}) σ^{a1}⊗⋯⊗σ^{ak} ⪰ 0`.

`sites` are the registry's internal 1-based Pauli site positions.
"""
function rdm_block(
    registry::VariableRegistry{PauliAlgebra,T},
    sites::AbstractVector{<:Integer},
) where {T<:Unsigned}
    cluster = Int.(collect(sites))
    _check_pauli_rdm_sites(registry, cluster)

    k = length(cluster)
    dim = 1 << k
    nwords = 4^k
    P = Polynomial{PauliAlgebra,T,ComplexF64}
    entries = [zero(P) for _ in 1:dim, _ in 1:dim]
    scale = 1.0 / dim

    for code in 0:(nwords - 1)
        ops = _pauli_ops_from_code(code, k)
        word_poly = _pauli_cluster_word(T, cluster, ops)
        for row in 1:dim, col in 1:dim
            coef = scale * _pauli_tensor_matrix_entry(ops, row, col)
            iszero(coef) || (entries[row, col] = entries[row, col] + coef * word_poly)
        end
    end

    return PhysicalPSDConstraint(
        :HPSD,
        entries;
        origin=PhysicalPSDOrigin(:spin_rdm, (Tuple(cluster),)),
    )
end

"""
    rdm_blocks(registry, clusters)

Build one spin-RDM PSD block per cluster.
"""
rdm_blocks(registry::VariableRegistry{PauliAlgebra}, clusters) =
    [rdm_block(registry, collect(cluster)) for cluster in clusters]
