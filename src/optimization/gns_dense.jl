"""
    DenseMomentSolution

Dense, non-sparse moment data for the phase-1 GNS path.

`DenseMomentSolution` packages the dense basis, dense moment values, and dense
Hankel matrix used by the phase-1 GNS reconstruction path.

# Fields
- `registry::VariableRegistry{A,TI}`: Variable registry used to generate the dense basis and interpret monomials or variable indices.
- `order::Int`: Relaxation order of the dense moment data. The dense basis runs up to this degree, while the stored moment dictionary covers moments needed through degree `2 * order`.
- `basis::Vector{M}`: Dense monomial basis up to `order`. Row `i` and column `j` of `hankel` are indexed by `basis[i]` and `basis[j]`.
- `moments::Dict{M,C}`: Dense moment dictionary keyed by monomials. These values are used to rebuild Hankel and shifted Gram matrices and to verify reconstructed moments.
- `hankel::Matrix{C}`: Dense Hankel (moment) matrix on `basis`, where entry `(i, j)` is the moment associated with `basis[i]' * basis[j]`.

# Example
See the worked example source at `docs/src/examples/literate/dense_gns_construction.jl`,
rendered as [Dense GNS Reconstruction](@ref dense-gns-construction).
"""
struct DenseMomentSolution{A<:MonoidAlgebra,TI<:Integer,M<:NormalMonomial{A,TI},C<:Number}
    registry::VariableRegistry{A,TI}
    order::Int
    basis::Vector{M}
    moments::Dict{M,C}
    hankel::Matrix{C}
end

"""
    GNSFlatnessReport

Machine-readable diagnostics for dense GNS reconstruction.

`GNSFlatnessReport` summarizes the numerical rank test behind flatness and the
moment-reproduction error measured for a reconstructed dense GNS model.

# Fields
- `order::Int`: Order of the source dense moment data.
- `quotient_order::Int`: Lower degree used to form the quotient Hankel block and select quotient representatives.
- `rank_full::Int`: Numerical rank of the full dense Hankel matrix.
- `rank_quotient::Int`: Numerical rank of the lower-degree Hankel block at `quotient_order`.
- `flat::Bool`: Whether the dense moment data is flat, i.e. `rank_full == rank_quotient`.
- `atol::Float64`: Absolute tolerance used in the numerical rank test.
- `rtol::Float64`: Relative tolerance used in the numerical rank test.
- `verified_degree::Int`: Highest monomial degree checked when comparing reconstructed moments against the source data.
- `max_moment_error::Float64`: Maximum absolute moment mismatch up to `verified_degree`. This is `NaN` on the initial flatness report before `verify_gns` fills in the verification error.

# Example
See the worked example source at `docs/src/examples/literate/dense_gns_construction.jl`,
rendered as [Dense GNS Reconstruction](@ref dense-gns-construction).
"""
struct GNSFlatnessReport
    order::Int
    quotient_order::Int
    rank_full::Int
    rank_quotient::Int
    flat::Bool
    atol::Float64
    rtol::Float64
    verified_degree::Int
    max_moment_error::Float64
end

"""
    GNSModel

Finite-dimensional dense GNS model with cyclic vector and verification report.

`GNSModel` is the finite-dimensional operator model recovered from dense moment
data. It stores both the quotient-space representation and the diagnostics used
to check that the reconstructed operators reproduce the source moments.

# Fields
- `registry::VariableRegistry{A,TI}`: Variable registry for the reconstructed generators.
- `order::Int`: Order of the source dense moment data.
- `quotient_order::Int`: Degree used to choose quotient representatives from the lower Hankel block.
- `basis::Vector{M}`: Full dense basis from the source moment data.
- `quotient_basis::Vector{M}`: Selected quotient representatives spanning the reconstructed Hilbert space.
- `gram::Matrix{C}`: Gram matrix of the selected `quotient_basis` vectors.
- `shifted_grams::Dict{TI,Matrix{C}}`: Shifted Gram matrices, keyed by variable index, obtained by multiplying the quotient basis by one generator.
- `operators::Dict{TI,Matrix{C}}`: Reconstructed multiplication operators, keyed by variable index, acting in the `quotient_basis` coordinates.
- `cyclic_vector::Vector{C}`: Canonical coordinates of the identity representative in the quotient basis.
- `report::GNSFlatnessReport`: Flatness and verification diagnostics attached to this reconstructed model.

# Example
See the worked example source at `docs/src/examples/literate/dense_gns_construction.jl`,
rendered as [Dense GNS Reconstruction](@ref dense-gns-construction).
"""
struct GNSModel{A<:MonoidAlgebra,TI<:Integer,M<:NormalMonomial{A,TI},C<:Number}
    registry::VariableRegistry{A,TI}
    order::Int
    quotient_order::Int
    basis::Vector{M}
    quotient_basis::Vector{M}
    gram::Matrix{C}
    shifted_grams::Dict{TI,Matrix{C}}
    operators::Dict{TI,Matrix{C}}
    cyclic_vector::Vector{C}
    report::GNSFlatnessReport
end

_numeric_value(x::Number) = x
_numeric_value(x) = value(x)

function _rank_tolerance(svals::AbstractVector, atol::Real, rtol::Real)
    isempty(svals) && return Float64(atol)
    return Float64(max(atol, rtol * first(svals)))
end

function _numerical_rank(mat::AbstractMatrix; atol::Real=1e-8, rtol::Real=1e-8)
    svals = svdvals(mat)
    tol = _rank_tolerance(svals, atol, rtol)
    return count(>(tol), svals), tol
end

_lower_basis_indices(basis, quotient_order::Int) = findall(m -> degree(m) <= quotient_order, basis)

function _dense_hankel_key(
    row::NormalMonomial{A,T},
    col::NormalMonomial{A,T},
) where {A<:MonoidAlgebra,T<:Integer}
    return NormalMonomial{A,T}(simplify(A, neat_dot(row.word, col.word)))
end

function _dense_hankel_matrix(
    basis::Vector{M},
    moments::Dict{M,C},
) where {A<:MonoidAlgebra,TI<:Integer,M<:NormalMonomial{A,TI},C<:Number}
    hankel = Matrix{C}(undef, length(basis), length(basis))
    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            hankel[i, j] = get(moments, _dense_hankel_key(row_mono, col_mono), zero(C))
        end
    end
    return hankel
end

function _dense_moment_dictionary(
    registry::VariableRegistry{A,TI},
    order::Int,
    raw_monomap,
) where {A<:MonoidAlgebra,TI<:Integer}
    basis = get_ncbasis(registry, 2 * order)
    C = eltype(values(raw_monomap))
    moments = Dict{eltype(basis),C}()
    for mono in basis
        key = symmetric_canon(expval(mono))
        moments[mono] = get(raw_monomap, key, zero(C))
    end
    return moments
end

function _dense_moment_solution_from_hankel(
    hankel::Matrix{C},
    registry::VariableRegistry{A,TI},
    order::Int,
) where {C<:Number,A<:MonoidAlgebra,TI<:Integer}
    basis = get_ncbasis(registry, order)
    size(hankel, 1) == length(basis) == size(hankel, 2) ||
        throw(ArgumentError("Hankel matrix size $(size(hankel)) does not match dense basis length $(length(basis)) for order $order."))
    moments = hankel_entries_dict(hankel, basis)
    return DenseMomentSolution(registry, order, basis, moments, Matrix(hankel))
end

function _supports_dense_gns(sparsity::SparsityResult)
    length(sparsity.corr_sparsity.cliques) == 1 || return false
    isempty(sparsity.corr_sparsity.global_cons) || return false
    length(sparsity.cliques_term_sparsities) == 1 || return false
    return all(ts -> length(ts.block_bases) == 1, only(sparsity.cliques_term_sparsities))
end

function _dense_moment_solution_from_solve(
    sparsity::SparsityResult{A,TI,P,M,Nothing},
    relaxation_spec,
    monomap,
) where {A<:MonoidAlgebra,TI<:Integer,P,M<:NormalMonomial{A,TI}}
    monomap === nothing && return nothing
    relaxation_spec isa Int || return nothing
    relaxation_spec > 0 || return nothing
    _supports_dense_gns(sparsity) || return nothing

    raw_monomap = Dict(copy(key) => _numeric_value(val) for (key, val) in pairs(monomap))
    registry = sparsity.corr_sparsity.registry
    basis = get_ncbasis(registry, relaxation_spec)
    moments = _dense_moment_dictionary(registry, relaxation_spec, raw_monomap)
    hankel = _dense_hankel_matrix(basis, moments)
    return DenseMomentSolution(registry, relaxation_spec, basis, moments, hankel)
end

_dense_moment_solution_from_solve(::SparsityResult, _relaxation_spec, _monomap) = nothing

function dense_moment_solution(result::PolyOptResult)
    isnothing(result.dense_moment_solution) || return result.dense_moment_solution
    throw(ArgumentError(
        "Dense moment data is unavailable. Use `cs_nctssos(...; dualize=false)` with a single dense moment block over `NonCommutativeAlgebra`, `ProjectorAlgebra`, or `UnipotentAlgebra`.",
    ))
end

function flatness_report(
    sol::DenseMomentSolution;
    quotient_order::Int=sol.order - 1,
    atol::Real=1e-8,
    rtol::Real=1e-8,
)
    0 <= quotient_order < sol.order || throw(ArgumentError(
        "`quotient_order` must satisfy 0 <= quotient_order < order=$(sol.order).",
    ))

    lower_idx = _lower_basis_indices(sol.basis, quotient_order)
    lower_hankel = sol.hankel[lower_idx, lower_idx]
    rank_full, _ = _numerical_rank(sol.hankel; atol, rtol)
    rank_quotient, _ = _numerical_rank(lower_hankel; atol, rtol)

    return GNSFlatnessReport(
        sol.order,
        quotient_order,
        rank_full,
        rank_quotient,
        rank_full == rank_quotient,
        Float64(atol),
        Float64(rtol),
        min(2 * quotient_order + 1, 2 * sol.order),
        NaN,
    )
end

function _select_quotient_basis(
    lower_hankel::AbstractMatrix,
    target_rank::Int,
    tol::Real,
)
    selected = Int[]
    for col_idx in axes(lower_hankel, 2)
        residual = if isempty(selected)
            lower_hankel[col_idx, col_idx]
        else
            gram = lower_hankel[selected, selected]
            column = lower_hankel[selected, col_idx]
            lower_hankel[col_idx, col_idx] - dot(column, gram \ column)
        end

        if abs(residual) > tol
            push!(selected, col_idx)
        end
        length(selected) == target_rank && break
    end

    length(selected) == target_rank || throw(ArgumentError(
        "Could not select $target_rank independent quotient representatives from the lower-degree Hankel block.",
    ))
    return selected
end

function _cyclic_vector(
    quotient_basis::Vector{M},
    ::Type{C},
) where {A<:MonoidAlgebra,TI<:Integer,M<:NormalMonomial{A,TI},C<:Number}
    first(quotient_basis) == one(first(quotient_basis)) || throw(ArgumentError(
        "The quotient basis must start with the identity representative so the cyclic vector can be fixed canonically.",
    ))
    cyclic_vector = zeros(C, length(quotient_basis))
    cyclic_vector[1] = one(C)
    return cyclic_vector
end

function _evaluate_moment(
    model::GNSModel{A,TI,M,C},
    mono::NormalMonomial{A,TI},
) where {A<:MonoidAlgebra,TI<:Integer,M<:NormalMonomial{A,TI},C<:Number}
    vector = copy(model.cyclic_vector)
    for var_idx in reverse(mono.word)
        vector = model.operators[var_idx] * vector
    end
    return dot(model.cyclic_vector, model.gram * vector)
end

function verify_gns(
    model::GNSModel{A,TI,M,C},
    sol::DenseMomentSolution{A,TI,M,C};
    max_degree::Int=min(2 * model.quotient_order + 1, 2 * sol.order),
) where {A<:MonoidAlgebra,TI<:Integer,M<:NormalMonomial{A,TI},C<:Number}
    0 <= max_degree <= 2 * sol.order || throw(ArgumentError(
        "`max_degree` must satisfy 0 <= max_degree <= $(2 * sol.order).",
    ))

    max_error = 0.0
    for mono in get_ncbasis(sol.registry, max_degree)
        expected = get(sol.moments, mono, zero(C))
        actual = _evaluate_moment(model, mono)
        max_error = max(max_error, Float64(abs(actual - expected)))
    end

    report = model.report
    return GNSFlatnessReport(
        report.order,
        report.quotient_order,
        report.rank_full,
        report.rank_quotient,
        report.flat,
        report.atol,
        report.rtol,
        max_degree,
        max_error,
    )
end

function gns_reconstruct(
    sol::DenseMomentSolution{A,TI,M,C};
    quotient_order::Int=sol.order - 1,
    atol::Real=1e-8,
    rtol::Real=1e-8,
    require_flat::Bool=false,
) where {A<:MonoidAlgebra,TI<:Integer,M<:NormalMonomial{A,TI},C<:Number}
    report = flatness_report(sol; quotient_order, atol, rtol)
    require_flat && !report.flat && throw(ArgumentError(
        "Flatness failed: rank(H)=$(report.rank_full), rank(H_≤$(report.quotient_order))=$(report.rank_quotient).",
    ))
    report.rank_quotient > 0 || throw(ArgumentError("The lower-degree Hankel block has rank zero."))

    lower_idx = _lower_basis_indices(sol.basis, quotient_order)
    lower_basis = sol.basis[lower_idx]
    lower_hankel = sol.hankel[lower_idx, lower_idx]
    _, tol = _numerical_rank(lower_hankel; atol, rtol)
    selected = _select_quotient_basis(lower_hankel, report.rank_quotient, tol)

    quotient_basis = lower_basis[selected]
    gram = lower_hankel[selected, selected]
    hankel_dict = hankel_entries_dict(sol.hankel, sol.basis)

    shifted_grams = Dict{TI,Matrix{C}}()
    operators = Dict{TI,Matrix{C}}()
    for var_idx in indices(sol.registry)
        shifted = construct_localizing_matrix(hankel_dict, var_idx, quotient_basis)
        shifted_grams[var_idx] = shifted
        operators[var_idx] = gram \ shifted
    end

    cyclic_vector = _cyclic_vector(quotient_basis, C)
    model = GNSModel(
        sol.registry,
        sol.order,
        quotient_order,
        sol.basis,
        quotient_basis,
        gram,
        shifted_grams,
        operators,
        cyclic_vector,
        report,
    )

    verified_report = verify_gns(model, sol)
    return GNSModel(
        model.registry,
        model.order,
        model.quotient_order,
        model.basis,
        model.quotient_basis,
        model.gram,
        model.shifted_grams,
        model.operators,
        model.cyclic_vector,
        verified_report,
    )
end

function gns_reconstruct(
    hankel::Matrix{C},
    registry::VariableRegistry{A,TI},
    order::Int;
    kwargs...,
) where {C<:Number,A<:MonoidAlgebra,TI<:Integer}
    return gns_reconstruct(_dense_moment_solution_from_hankel(hankel, registry, order); kwargs...)
end

gns_reconstruct(result::PolyOptResult; kwargs...) = gns_reconstruct(dense_moment_solution(result); kwargs...)

function reconstruct(
    hankel::Matrix{C},
    registry::VariableRegistry{A,TI},
    order::Int;
    atol::Float64=1e-3,
    rtol::Float64=1e-8,
    require_flat::Bool=false,
) where {C<:Number,A<:MonoidAlgebra,TI<:Integer}
    return gns_reconstruct(hankel, registry, order; atol, rtol, require_flat).operators
end
