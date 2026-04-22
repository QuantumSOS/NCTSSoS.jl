using LinearAlgebra

# Types from NCTSSoS main module (no submodule needed)

"""
    GNSResult{T,RT,A,TI,M}

Result of a GNS (Gelfand-Naimark-Segal) reconstruction.

# Fields
- `matrices::Dict{TI,Matrix{T}}`: Reconstructed multiplication operators for each variable
- `xi::Vector{T}`: Distinguished vector corresponding to the class of the identity word
- `basis::Vector{M}`: Basis words used for the quotient space / principal Hankel block
- `full_basis::Vector{M}`: Basis words used for the full Hankel matrix
- `singular_values::Vector{RT}`: Singular values of the principal Hankel block
- `rank::Int`: Numerical rank of the principal Hankel block
- `full_rank::Int`: Numerical rank of the full Hankel matrix

The matrices are expressed in an orthonormal basis obtained from the principal
Hankel block via SVD. This is the numerically robust analogue of the column-space
construction from the notes: the resulting operators are equivalent up to an
orthogonal/unitary change of basis.
"""
struct GNSResult{T,RT<:Real,A<:AlgebraType,TI<:Integer,M<:NormalMonomial{A,TI}}
    matrices::Dict{TI,Matrix{T}}
    xi::Vector{T}
    basis::Vector{M}
    full_basis::Vector{M}
    singular_values::Vector{RT}
    rank::Int
    full_rank::Int
end

"""
    _gns_extract_monomials_from_basis(basis_polys::Vector{Polynomial{A,T,C}}) where {A,T,C}

Extract monomials from a vector of single-term polynomials for GNS reconstruction.

For NonCommutative, Pauli, Projector, and Unipotent algebras, each basis polynomial
has exactly one term. This function extracts the monomial from each polynomial.

For Fermionic/Bosonic algebras, basis polynomials may have multiple terms due to
normal-ordering corrections. In those cases, we keep the leading monomial.

This is a GNS-internal helper distinct from basis extraction used elsewhere in the
sparsity pipeline.
"""
function _gns_extract_monomials_from_basis(basis_monos::Vector{NormalMonomial{A,T}}) where {A,T}
    return basis_monos
end

function _gns_extract_monomials_from_basis(basis_polys::Vector{Polynomial{A,T,C}}) where {A<:AlgebraType,T<:Integer,C<:Number}
    result = NormalMonomial{A,T}[]
    for p in basis_polys
        isempty(p.terms) && continue
        push!(result, p.terms[1].monomial)
    end
    return result
end

@inline function _gns_zero_value(monomap::AbstractDict)
    if valtype(monomap) === Any
        isempty(monomap) && return 0.0
        return zero(first(values(monomap)))
    end
    return zero(valtype(monomap))
end

@inline function _gns_basis_indices(full_basis::Vector{M}, basis::Vector{M}) where {M}
    full_basis_to_idx = Dict(m => i for (i, m) in enumerate(full_basis))
    indices = Int[]
    sizehint!(indices, length(basis))
    for mono in basis
        haskey(full_basis_to_idx, mono) ||
            throw(ArgumentError("Basis word $mono is missing from the full Hankel basis."))
        push!(indices, full_basis_to_idx[mono])
    end
    return indices
end

@inline function _gns_rank(svals::AbstractVector{RT}; atol::Real=1e-8, rtol::Real=0.0) where {RT<:Real}
    isempty(svals) && return (0, zero(RT))
    cutoff = max(RT(atol), RT(rtol) * maximum(svals))
    return (count(>(cutoff), svals), cutoff)
end

function _gns_validate_degrees(H_deg::Int, hankel_deg::Int)
    H_deg <= 0 && throw(ArgumentError("`H_deg` must be positive, got $H_deg."))
    hankel_deg < 0 && throw(ArgumentError("`hankel_deg` must be non-negative, got $hankel_deg."))
    hankel_deg >= H_deg && throw(ArgumentError("`hankel_deg` must be strictly smaller than `H_deg`; got hankel_deg=$hankel_deg, H_deg=$H_deg."))
    return nothing
end

function _gns_basis_pair(registry::VariableRegistry{A,TI}, H_deg::Int, hankel_deg::Int) where {A<:AlgebraType,TI<:Integer}
    full_basis = _gns_extract_monomials_from_basis(get_ncbasis(registry, H_deg))
    basis = _gns_extract_monomials_from_basis(get_ncbasis(registry, hankel_deg))
    return full_basis, basis
end

@inline _gns_moment_key(mono::NormalMonomial) = symmetric_canon(expval(mono))

@inline function _gns_moment_value(monomap::AbstractDict, mono::NormalMonomial)
    key = _gns_moment_key(mono)
    if haskey(monomap, key)
        return monomap[key]
    end
    return _gns_zero_value(monomap)
end

function _gns_moment_from_word(
    monomap::AbstractDict,
    ::Type{A},
    word::Vector{T},
) where {A<:AlgebraType,T<:Integer}
    total = _gns_zero_value(monomap)
    for (coef, mono) in _simplified_to_terms(A, simplify(A, word), T)
        total += coef * _gns_moment_value(monomap, mono)
    end
    return total
end

function _gns_missing_moment_keys(
    monomap::AbstractDict,
    basis::Vector{M},
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T}}
    missing = Any[]
    seen = Set{Any}()

    for row_mono in basis, col_mono in basis
        for (_coef, mono) in _simplified_to_terms(A, simplify(A, neat_dot(row_mono.word, col_mono.word)), T)
            key = _gns_moment_key(mono)
            if !haskey(monomap, key) && !(key in seen)
                push!(seen, key)
                push!(missing, key)
            end
        end
    end

    return missing
end

function _gns_validate_monomap_coverage(
    monomap::AbstractDict,
    basis::Vector{M},
    H_deg::Int,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T}}
    missing = _gns_missing_moment_keys(monomap, basis)
    isempty(missing) && return nothing

    preview_count = min(length(missing), 5)
    preview = join(repr.(missing[1:preview_count]), ", ")
    suffix = length(missing) > preview_count ? ", ..." : ""
    throw(ArgumentError(
        "`monomap` is missing $(length(missing)) moment value(s) required for GNS reconstruction up to H_deg=$H_deg. " *
        "First missing keys: $preview$suffix. Pass a dense Hankel matrix instead, or provide a `monomap` covering all moments up to the requested degree."
    ))
end

"""
    hankel_matrix(monomap::AbstractDict, basis::Vector{<:NormalMonomial})
    hankel_matrix(monomap::AbstractDict, registry::VariableRegistry, degree::Int)

Build a Hankel / moment matrix from solved moment values.

The entry at `(i, j)` is `L(basis[i]' * basis[j])`, where `L` is encoded by
`monomap`.
"""
function hankel_matrix(monomap::AbstractDict, basis::Vector{M}) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T}}
    return [
        _gns_moment_from_word(monomap, A, neat_dot(row_mono.word, col_mono.word))
        for row_mono in basis, col_mono in basis
    ]
end

function hankel_matrix(monomap::AbstractDict, registry::VariableRegistry{A,TI}, degree::Int) where {A<:AlgebraType,TI<:Integer}
    basis = _gns_extract_monomials_from_basis(get_ncbasis(registry, degree))
    return hankel_matrix(monomap, basis)
end

"""
    construct_localizing_matrix(monomap::AbstractDict, var_idx::TI, basis::Vector{M})

Construct the localizing matrix for left multiplication by a variable using a
moment lookup dictionary.

The entry `(i, j)` equals `L(basis[i]' * x_var * basis[j])`.
"""
function construct_localizing_matrix(
    monomap::AbstractDict,
    var_idx::TI,
    basis::Vector{M},
) where {A<:AlgebraType,TI<:Integer,T<:Integer,M<:NormalMonomial{A,T}}
    var_word = T[T(var_idx)]
    return [
        _gns_moment_from_word(monomap, A, _neat_dot3(row_mono.word, var_word, col_mono.word))
        for row_mono in basis, col_mono in basis
    ]
end

"""
    hankel_entries_dict(hankel::Matrix{T}, basis::Vector{<:NormalMonomial}) where {T <: Number}

Create a lookup dictionary mapping canonical monomials to Hankel matrix entries.

This helper is mainly useful for tests and debugging. If multiple basis pairs
produce the same canonical monomial, the first encountered entry is kept; for
valid Hankel data these entries should agree.
"""
function hankel_entries_dict(hankel::Matrix{T}, basis::Vector{M}) where {A<:MonoidAlgebra,TM<:Integer,M<:NormalMonomial{A,TM},T<:Number}
    size(hankel, 1) == size(hankel, 2) ||
        throw(ArgumentError("Hankel matrix must be square, got size $(size(hankel))."))
    length(basis) == size(hankel, 1) || throw(
        ArgumentError(
            "Basis length $(length(basis)) must match Hankel size $(size(hankel, 1)).",
        ),
    )

    dict = Dict{M,T}()
    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            key = M(simplify(A, neat_dot(row_mono.word, col_mono.word)))
            get!(dict, key, hankel[i, j])
        end
    end

    return dict
end

"""
    construct_localizing_matrix(hankel_dict::Dict{M,T}, var_idx::TI, basis::Vector{M})

Construct the localizing matrix for left multiplication by a variable using a
precomputed Hankel dictionary.

This method is currently defined for monoid algebras, where products reduce to a
single canonical monomial.
"""
function construct_localizing_matrix(
    hankel_dict::Dict{M,T},
    var_idx::TI,
    basis::Vector{M},
) where {A<:MonoidAlgebra,TM<:Integer,M<:NormalMonomial{A,TM},T<:Number,TI<:Integer}
    n = length(basis)
    K = zeros(T, n, n)
    var_word = TM[TM(var_idx)]

    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            key = M(simplify(A, _neat_dot3(row_mono.word, var_word, col_mono.word)))
            K[i, j] = get(hankel_dict, key, zero(T))
        end
    end

    return K
end

function _gns_localizing_from_hankel(
    hankel::AbstractMatrix{T},
    full_basis::Vector{M},
    var_idx::TI,
    basis::Vector{M},
) where {A<:AlgebraType,TM<:Integer,TI<:Integer,M<:NormalMonomial{A,TM},T<:Number}
    basis_to_idx = Dict(m => i for (i, m) in enumerate(full_basis))
    row_indices = Int[]
    sizehint!(row_indices, length(basis))
    for mono in basis
        haskey(basis_to_idx, mono) ||
            throw(ArgumentError("Basis word $mono is missing from the full Hankel basis."))
        push!(row_indices, basis_to_idx[mono])
    end
    var_mono = NormalMonomial{A,TM}(TM[TM(var_idx)])

    K = zeros(T, length(basis), length(basis))
    for (i, row_idx) in enumerate(row_indices)
        for (j, col_mono) in enumerate(basis)
            entry = zero(T)
            prod = var_mono * col_mono
            for (coef, mono) in terms(prod)
                haskey(basis_to_idx, mono) ||
                    throw(ArgumentError("Product $(var_mono) * $(col_mono) is missing from the full Hankel basis."))
                entry += coef * hankel[row_idx, basis_to_idx[mono]]
            end
            K[i, j] = entry
        end
    end

    return K
end

@inline _gns_postprocess_operator(::Type{A}, X::AbstractMatrix{T}) where {A<:AlgebraType,T} = Matrix{T}(X)
@inline _gns_postprocess_operator(::Type{NonCommutativeAlgebra}, X::AbstractMatrix{T}) where {T} = Matrix{T}((X + X') / 2)
@inline _gns_postprocess_operator(::Type{PauliAlgebra}, X::AbstractMatrix{T}) where {T} = Matrix{T}((X + X') / 2)
@inline _gns_postprocess_operator(::Type{ProjectorAlgebra}, X::AbstractMatrix{T}) where {T} = Matrix{T}((X + X') / 2)
@inline _gns_postprocess_operator(::Type{UnipotentAlgebra}, X::AbstractMatrix{T}) where {T} = Matrix{T}((X + X') / 2)

function _gns_finalize(
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

    F = svd(hankel_block)
    rank_basis, cutoff = _gns_rank(F.S; atol=atol, rtol=rtol)
    rank_basis == 0 && throw(
        ArgumentError(
            "No singular values exceed the numerical cutoff $cutoff. Consider decreasing `atol`/`rtol`.",
        ),
    )

    full_svals = svdvals(Matrix(hankel))
    rank_full, _ = _gns_rank(full_svals; atol=atol, rtol=rtol)

    if rank_full != rank_basis
        @warn "Flatness condition violated: rank(H) = $rank_full ≠ rank(hankel_block) = $rank_basis. The full Hankel matrix is not a flat extension of the principal block."
    end

    U = F.U[:, 1:rank_basis]
    S = F.S[1:rank_basis]
    scale = Diagonal(inv.(sqrt.(S)))
    Trec = promote_type(eltype(hankel), eltype(U), eltype(scale))

    matrices = Dict{TI,Matrix{Trec}}()
    for var_idx in indices(registry)
        K = _gns_localizing_from_hankel(hankel, full_basis, var_idx, basis)
        X = scale * (U' * K * U) * scale
        matrices[var_idx] = _gns_postprocess_operator(A, X)
    end

    e1 = zeros(Trec, length(basis))
    e1[1] = one(Trec)
    xi = Diagonal(sqrt.(S)) * (U' * e1)

    return GNSResult{Trec,eltype(F.S),A,TI,M}(
        matrices,
        Vector{Trec}(xi),
        basis,
        full_basis,
        collect(F.S),
        rank_basis,
        rank_full,
    )
end

"""
    gns_reconstruct(hankel::AbstractMatrix, registry::VariableRegistry, H_deg::Int;
                    method::Symbol=:svd, hankel_deg::Int=H_deg-1, atol::Real=1e-8, rtol::Real=0.0)

Perform GNS reconstruction from a full Hankel matrix.

`H_deg` is the degree of the full Hankel matrix (use `H_deg = order` where
`order` is the SDP relaxation order). `hankel_deg` is the degree of its
principal block used to define the quotient space (standard choice:
`hankel_deg = H_deg - 1`).

The result contains the reconstructed matrices, the distinguished vector `ξ`,
and rank information.

!!! tip "Recommended workflow"
    Raw solver Hankel matrices are typically not flat. For best results,
    pre-process with [`test_flatness`](@ref) and [`flat_extend`](@ref):

    ```julia
    hankel_flat = flat_extend(hankel, full_basis, basis; atol=1e-6)
    gns = gns_reconstruct(hankel_flat, registry, H_deg; ...)
    ```

    The `monomap` overload (below) skips this step, which is fine for
    `:svd` but may degrade `:cholesky` accuracy.

!!! note "`atol` guidance"
    The default `atol=1e-8` is well-tuned for standard SDP solvers
    (COSMO, Mosek). Setting `atol` below the solver noise floor
    (~1e-9) will inflate the detected rank with spurious directions.
"""
function gns_reconstruct(
    hankel::AbstractMatrix{T},
    registry::VariableRegistry{A,TI},
    H_deg::Int;
    method::Symbol=:svd,
    hankel_deg::Int=H_deg - 1,
    atol::Real=1e-8,
    rtol::Real=0.0,
) where {T<:Number,A<:AlgebraType,TI<:Integer}
    _gns_validate_degrees(H_deg, hankel_deg)
    full_basis, basis = _gns_basis_pair(registry, H_deg, hankel_deg)
    if method == :svd
        return _gns_finalize(hankel, full_basis, basis, registry; atol=atol, rtol=rtol)
    elseif method == :cholesky
        return _gns_finalize_cholesky(hankel, full_basis, basis, registry; atol=atol, rtol=rtol)
    end
    throw(ArgumentError("Unknown GNS method: $method. Use :svd or :cholesky."))
end

"""
    gns_reconstruct(monomap::AbstractDict, registry::VariableRegistry, H_deg::Int;
                    method::Symbol=:svd, hankel_deg::Int=H_deg-1, atol::Real=1e-8, rtol::Real=0.0)

Perform GNS reconstruction directly from solved moment values.

This is the most convenient entry point after `solve_moment_problem`: pass
`result.monomap`, the registry, and a full Hankel degree (typically `H_deg =
order` where `order` is the SDP relaxation order). The full Hankel matrix is
assembled automatically.

The `monomap` keys are `StateSymbol` objects of the form
`symmetric_canon(expval(mono))` — the same format returned by
`solve_moment_problem(...).monomap`.

The input `monomap` must contain every moment needed to build the dense Hankel
matrix; sparse/custom bases should use the matrix overload instead. In
particular, symmetry-reduced `monomap`s from `SolverConfig(; symmetry=...)`
are not directly compatible today because they omit many dense moments needed
by GNS reconstruction.
"""
function gns_reconstruct(
    monomap::AbstractDict,
    registry::VariableRegistry{A,TI},
    H_deg::Int;
    method::Symbol=:svd,
    hankel_deg::Int=H_deg - 1,
    atol::Real=1e-8,
    rtol::Real=0.0,
) where {A<:AlgebraType,TI<:Integer}
    _gns_validate_degrees(H_deg, hankel_deg)
    full_basis, basis = _gns_basis_pair(registry, H_deg, hankel_deg)
    _gns_validate_monomap_coverage(monomap, full_basis, H_deg)
    hankel = hankel_matrix(monomap, full_basis)
    if method == :svd
        return _gns_finalize(hankel, full_basis, basis, registry; atol=atol, rtol=rtol)
    elseif method == :cholesky
        return _gns_finalize_cholesky(hankel, full_basis, basis, registry; atol=atol, rtol=rtol)
    end
    throw(ArgumentError("Unknown GNS method: $method. Use :svd or :cholesky."))
end

@inline _gns_abs_index(idx::Integer) = abs(Int(idx))

function _gns_word_in_clique(word::AbstractVector{<:Integer}, clique_set::Set{Int})
    for idx in word
        _gns_abs_index(idx) in clique_set || return false
    end
    return true
end

function _gns_project_word_to_clique(word::AbstractVector{TI}, clique_set::Set{Int}) where {TI<:Integer}
    projected = TI[]
    sizehint!(projected, length(word))
    for idx in word
        _gns_abs_index(idx) in clique_set && push!(projected, idx)
    end
    return projected
end

function _gns_filter_monomap_to_clique(monomap::AbstractDict, clique_set::Set{Int})
    filtered = Dict{keytype(monomap),valtype(monomap)}()
    for (key, value) in pairs(monomap)
        if key isa NormalMonomial
            _gns_word_in_clique(key.word, clique_set) || continue
        elseif key isa AbstractVector{<:Integer}
            _gns_word_in_clique(key, clique_set) || continue
        else
            continue
        end
        filtered[key] = value
    end
    return filtered
end

function _gns_pairwise_disjoint_cliques(cliques)
    seen = Set{Int}()
    for clique in cliques
        clique_set = Set(Int.(clique))
        isempty(intersect(seen, clique_set)) || return false
        union!(seen, clique_set)
    end
    return true
end

function _gns_validate_sparse_registry_coverage(corr_sparsity)
    registry_vars = Set(_gns_abs_index.(indices(corr_sparsity.registry)))
    covered_vars = Set{Int}()
    for clique in corr_sparsity.cliques
        union!(covered_vars, _gns_abs_index.(clique))
    end

    inactive_vars = sort!(collect(setdiff(registry_vars, covered_vars)))
    isempty(inactive_vars) || throw(ArgumentError(
        "Sparse GNS reconstruction requires every registry variable to appear in at least one clique. Inactive variable indices: $(inactive_vars)."
    ))
    return nothing
end

function _gns_sparse_clique_data(
    monomap::AbstractDict,
    corr_sparsity::CorrelativeSparsity{A,TI,P,M,Nothing},
    H_deg::Int,
) where {A<:AlgebraType,TI<:Integer,C<:Number,P<:Polynomial{A,TI,C},M<:NormalMonomial{A,TI}}
    data = NamedTuple[]

    for clique in corr_sparsity.cliques
        clique_set = Set(Int.(clique))
        clique_monomap = _gns_filter_monomap_to_clique(monomap, clique_set)
        clique_registry = _basis_subregistry(corr_sparsity.registry, clique)
        clique_full_basis = _gns_extract_monomials_from_basis(get_ncbasis(clique_registry, H_deg))
        _gns_validate_monomap_coverage(clique_monomap, clique_full_basis, H_deg)
        push!(data, (clique_set=clique_set, monomap=clique_monomap))
    end

    return data
end

function _gns_sparse_disjoint_moment_value(
    clique_data,
    ::Type{A},
    mono::NormalMonomial{A,TI},
) where {A<:AlgebraType,TI<:Integer}
    isempty(clique_data) && return one(Float64)

    first_data = first(clique_data)
    first_word = _gns_project_word_to_clique(mono.word, first_data.clique_set)
    total = _gns_moment_from_word(first_data.monomap, A, first_word)

    for data in Iterators.drop(clique_data, 1)
        projected = _gns_project_word_to_clique(mono.word, data.clique_set)
        total *= _gns_moment_from_word(data.monomap, A, projected)
    end

    return total
end

function _gns_sparse_disjoint_moment_from_word(
    clique_data,
    ::Type{A},
    word::Vector{TI},
) where {A<:AlgebraType,TI<:Integer}
    isempty(clique_data) && return 0.0

    total = _gns_zero_value(clique_data[1].monomap)
    for (coef, mono) in _simplified_to_terms(A, simplify(A, word), TI)
        total += coef * _gns_sparse_disjoint_moment_value(clique_data, A, mono)
    end
    return total
end

function _gns_sparse_disjoint_hankel(clique_data, full_basis::Vector{M}) where {A<:AlgebraType,TI<:Integer,M<:NormalMonomial{A,TI}}
    return [
        _gns_sparse_disjoint_moment_from_word(clique_data, A, neat_dot(row_mono.word, col_mono.word))
        for row_mono in full_basis, col_mono in full_basis
    ]
end

"""
    gns_reconstruct(monomap::AbstractDict, sparsity::SparsityResult, H_deg::Int;
                    hankel_deg::Int=H_deg-1, atol::Real=1e-8, rtol::Real=0.0)

Perform sparse GNS reconstruction from solved sparse moments and their sparsity
metadata.

Current support is intentionally conservative:
- a single clique falls back to dense [`gns_reconstruct`](@ref), and
- multiple cliques are supported only when they are pairwise disjoint and there
  are no global constraints.

In the disjoint-clique case, missing cross-clique moments are completed by the
product functional induced by the clique-local sparse moments, and dense GNS is
then applied to the resulting Hankel matrix. Overlapping-clique amalgamation
from SparseGNS / Theorem 4.2 is not implemented yet.
"""
function gns_reconstruct(
    monomap::AbstractDict,
    sparsity::SparsityResult{A,TI,P,M,Nothing},
    H_deg::Int;
    hankel_deg::Int=H_deg - 1,
    atol::Real=1e-8,
    rtol::Real=0.0,
) where {A<:AlgebraType,TI<:Integer,C<:Number,P<:Polynomial{A,TI,C},M<:NormalMonomial{A,TI}}
    _gns_validate_degrees(H_deg, hankel_deg)

    corr_sparsity = sparsity.corr_sparsity
    registry = corr_sparsity.registry
    n_cliques = length(corr_sparsity.cliques)

    n_cliques == 0 && throw(ArgumentError("Sparse GNS reconstruction requires at least one clique; got none."))
    _gns_validate_sparse_registry_coverage(corr_sparsity)
    n_cliques == 1 && return gns_reconstruct(
        monomap,
        registry,
        H_deg;
        hankel_deg=hankel_deg,
        atol=atol,
        rtol=rtol,
    )

    isempty(corr_sparsity.global_cons) || throw(ArgumentError(
        "Sparse GNS reconstruction currently requires `global_cons` to be empty for multi-clique problems; got $(length(corr_sparsity.global_cons)) global constraint(s)."
    ))
    _gns_pairwise_disjoint_cliques(corr_sparsity.cliques) || throw(ArgumentError(
        "Sparse GNS reconstruction currently supports only pairwise-disjoint cliques. Overlapping-clique amalgamation from SparseGNS / Theorem 4.2 is not implemented yet."
    ))

    clique_data = _gns_sparse_clique_data(monomap, corr_sparsity, H_deg)
    full_basis, basis = _gns_basis_pair(registry, H_deg, hankel_deg)
    hankel = _gns_sparse_disjoint_hankel(clique_data, full_basis)
    return _gns_finalize(hankel, full_basis, basis, registry; atol=atol, rtol=rtol)
end

"""
    reconstruct(H::AbstractMatrix, registry::VariableRegistry, H_deg::Int;
                hankel_deg::Int=H_deg-1, atol::Real=1e-8, rtol::Real=0.0)

Compatibility wrapper returning only the reconstructed operator matrices.

For the full GNS output, including the distinguished vector `ξ`, use
[`gns_reconstruct`](@ref).
"""
function reconstruct(
    hankel::AbstractMatrix,
    registry::VariableRegistry,
    H_deg::Int;
    hankel_deg::Int=H_deg - 1,
    atol::Real=1e-8,
    rtol::Real=0.0,
)
    return gns_reconstruct(
        hankel,
        registry,
        H_deg;
        hankel_deg=hankel_deg,
        atol=atol,
        rtol=rtol,
    ).matrices
end

function reconstruct(
    monomap::AbstractDict,
    sparsity::SparsityResult{A,TI,P,M,Nothing},
    H_deg::Int;
    hankel_deg::Int=H_deg - 1,
    atol::Real=1e-8,
    rtol::Real=0.0,
) where {A<:AlgebraType,TI<:Integer,C<:Number,P<:Polynomial{A,TI,C},M<:NormalMonomial{A,TI}}
    return gns_reconstruct(
        monomap,
        sparsity,
        H_deg;
        hankel_deg=hankel_deg,
        atol=atol,
        rtol=rtol,
    ).matrices
end
