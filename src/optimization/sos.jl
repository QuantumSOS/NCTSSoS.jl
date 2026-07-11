# =============================================================================
# SOS Problem Type
# =============================================================================

"""
    SOSProblem{T}

A Sum of Squares (SOS) optimization problem in JuMP form.

This is the dual formulation of a moment problem. The SOS problem is typically
smaller and faster to solve than the primal moment problem.

# Type Parameters
- `T`: Coefficient type (Float64 for real problems)

# Fields
- `model::GenericModel{T}`: The JuMP model ready to be optimized
- `n_unique_elements::Int`: Number of unique moment variables after canonicalization
- `psd_dual_blocks`: PSD dual variables with their symbolic block metadata
- `zero_duals`: scalar dual variables for zero constraints

# Usage
```julia
# From moment problem
mp = moment_relax(pop, corr_sparsity, term_sparsities)
sos = sos_dualize(mp)

# Solve
set_optimizer(sos.model, Clarabel.Optimizer)
optimize!(sos.model)
obj = objective_value(sos.model)
```

See also: [`sos_dualize`](@ref), [`MomentProblem`](@ref)
"""
struct SOSDualBlock{M}
    meta::BlockMeta{M}
    variable::Any
    representation::Symbol

    function SOSDualBlock{M}(meta::BlockMeta{M}, variable, representation::Symbol) where {M}
        representation in (:psd, :hermitian_lift, :hermitian) || throw(ArgumentError(
            "Unsupported SOS dual block representation $(repr(representation))."
        ))
        return new{M}(meta, variable, representation)
    end
end

SOSDualBlock(meta::BlockMeta{M}, variable, representation::Symbol) where {M} =
    SOSDualBlock{M}(meta, variable, representation)

struct SOSProblem{T}
    model::GenericModel{T}
    n_unique_elements::Int
    psd_dual_blocks::Vector{Any}
    zero_duals::Vector{Any}
end

SOSProblem(model::GenericModel{T}, n_unique_elements::Integer) where {T} =
    SOSProblem{T}(model, Int(n_unique_elements), Any[], Any[])

SOSProblem(
    model::GenericModel{T},
    n_unique_elements::Integer,
    psd_dual_blocks::AbstractVector,
    zero_duals::AbstractVector,
) where {T} = SOSProblem{T}(
    model,
    Int(n_unique_elements),
    Any[block for block in psd_dual_blocks],
    Any[dual for dual in zero_duals],
)

"""
    sos_dual_blocks(sos)

Return PSD dual variable handles together with their `BlockMeta` provenance.
For Hermitian SOS dualization, `representation == :hermitian_lift` means the
stored variable is the real lifted PSD matrix, while `:hermitian` keeps the
native Hermitian cone variable.
"""
sos_dual_blocks(sos::SOSProblem) = copy(sos.psd_dual_blocks)

"""
    sos_zero_duals(moment_data, sos)

Return zero/equality dual variable handles together with their origin label,
label feature fields, coefficient-domain metadata, and linear-form provenance.
`moment_data` may be either a symbolic `MomentProblem` or finalized
`MomentLinearData`.
"""
function sos_zero_duals(mp::MomentProblem, sos::SOSProblem)
    return sos_zero_duals(mp.linear, sos)
end

function sos_zero_duals(L::MomentLinearData, sos::SOSProblem)
    _assert_sos_dual_metadata_matches(L, sos)
    return [
        let label = _sos_origin_field(zc.origin, :label, nothing)
            (
                origin=zc.origin,
                label=label,
                feature=_sos_label_field(label, :feature, nothing),
                decomposition=_sos_label_field(label, :decomposition, nothing),
                reason=_sos_label_field(label, :reason, nothing),
                coefficient_domain=_sos_label_field(label, :coefficient_domain, nothing),
                exact_coefficient_domain=_sos_label_field(
                    label,
                    :exact_coefficient_domain,
                    nothing,
                ),
                form=zc.form,
                kind=zc.kind,
                term_count=length(zc.form),
                variable=dual,
            )
        end
        for (zc, dual) in zip(L.zero_constraints, sos.zero_duals)
    ]
end

"""
    sos_dual_block_values(sos)

Extract numeric PSD dual block values from a solved SOS model, preserving each
block's metadata and representation.
"""
_sos_origin_field(origin, field::Symbol, default) =
    hasproperty(origin, field) ? getproperty(origin, field) : default

_sos_transform_field(transform, field::Symbol, default) =
    transform === nothing ? default :
    hasproperty(transform, field) ? getproperty(transform, field) : default

_sos_label_field(label, field::Symbol, default) =
    label === nothing ? default :
    hasproperty(label, field) ? getproperty(label, field) : default

function _sos_native_dual_block_value(block::SOSDualBlock, value_matrix)
    if block.representation == :psd
        return value_matrix
    elseif block.representation == :hermitian
        return value_matrix
    elseif block.representation == :hermitian_lift
        n2 = size(value_matrix, 1)
        size(value_matrix, 2) == n2 && iseven(n2) || throw(DimensionMismatch(
            "Hermitian lifted SOS dual block must be an even square matrix; got $(size(value_matrix))."
        ))
        n = n2 ÷ 2
        native = Matrix{ComplexF64}(undef, n, n)
        for j in 1:n, i in 1:n
            x1 = value_matrix[i, j] + value_matrix[n + i, n + j]
            x2 = value_matrix[n + i, j] - value_matrix[i, n + j]
            native[i, j] = complex(Float64(real(x1)), Float64(real(x2)))
        end
        return native
    end
    throw(ArgumentError("Unsupported SOS dual block representation $(repr(block.representation))."))
end

function _sos_dual_block_value_record(block::SOSDualBlock)
    value_matrix = value.(block.variable)
    transform = _sos_origin_field(block.meta.origin, :transform, nothing)
    return (
        meta=block.meta,
        label=_sos_origin_field(block.meta.origin, :label, nothing),
        logical_row_labels=_sos_origin_field(block.meta.origin, :logical_row_labels, Any[]),
        transform=transform,
        transform_family=_sos_transform_field(transform, :family, nothing),
        coefficient_domain=_sos_transform_field(transform, :coefficient_domain, nothing),
        exact_coefficient_domain=_sos_transform_field(transform, :exact_coefficient_domain, nothing),
        representation=block.representation,
        value=value_matrix,
        native_value=_sos_native_dual_block_value(block, value_matrix),
    )
end

function sos_dual_block_values(sos::SOSProblem)
    return [
        _sos_dual_block_value_record(block)
        for block in sos.psd_dual_blocks
    ]
end

function sos_dual_block_values(mp::MomentProblem, sos::SOSProblem)
    return sos_dual_block_values(mp.linear, sos)
end

function sos_dual_block_values(L::MomentLinearData, sos::SOSProblem)
    _assert_sos_dual_metadata_matches(L, sos)
    return sos_dual_block_values(sos)
end

function _sos_hermitian_residual(matrix)
    isempty(matrix) && return 0.0
    return Float64(norm(matrix - matrix', Inf))
end

function _sos_min_eigenvalue(matrix)
    isempty(matrix) && return Inf
    hermitian_part = (matrix + matrix') / 2
    return Float64(real(eigmin(Hermitian(hermitian_part))))
end

function _sos_dual_block_diagnostic_record(block)
    return (
        meta=block.meta,
        label=block.label,
        logical_row_labels=block.logical_row_labels,
        transform=block.transform,
        transform_family=block.transform_family,
        coefficient_domain=block.coefficient_domain,
        exact_coefficient_domain=block.exact_coefficient_domain,
        representation=block.representation,
        value_size=size(block.value),
        native_value_size=size(block.native_value),
        native_hermitian_residual=_sos_hermitian_residual(block.native_value),
        native_min_eigenvalue=_sos_min_eigenvalue(block.native_value),
        lifted_min_eigenvalue=_sos_min_eigenvalue(block.value),
    )
end

"""
    sos_dual_block_diagnostics(block_values)
    sos_dual_block_diagnostics(sos)
    sos_dual_block_diagnostics(moment_data, sos)

Return per-block numerical diagnostics for solved SOS dual PSD variables,
including the native Hermitian block reconstructed from any real lift.
"""
function sos_dual_block_diagnostics(block_values::AbstractVector{<:NamedTuple})
    return [
        _sos_dual_block_diagnostic_record(block)
        for block in block_values
    ]
end

sos_dual_block_diagnostics(sos::SOSProblem) =
    sos_dual_block_diagnostics(sos_dual_block_values(sos))

function sos_dual_block_diagnostics(mp::MomentProblem, sos::SOSProblem)
    return sos_dual_block_diagnostics(mp.linear, sos)
end

function sos_dual_block_diagnostics(L::MomentLinearData, sos::SOSProblem)
    _assert_sos_dual_metadata_matches(L, sos)
    return sos_dual_block_diagnostics(sos)
end

"""
    sos_zero_dual_values(moment_data, sos)

Extract numeric zero/equality dual values from a solved SOS model, preserving
the origin label, label feature fields, coefficient-domain metadata, and linear
form of each zero constraint.
"""
function sos_zero_dual_values(mp::MomentProblem, sos::SOSProblem)
    return sos_zero_dual_values(mp.linear, sos)
end

function sos_zero_dual_values(L::MomentLinearData, sos::SOSProblem)
    return [
        (
            origin=dual.origin,
            label=dual.label,
            feature=dual.feature,
            decomposition=dual.decomposition,
            reason=dual.reason,
            coefficient_domain=dual.coefficient_domain,
            exact_coefficient_domain=dual.exact_coefficient_domain,
            form=dual.form,
            kind=dual.kind,
            term_count=dual.term_count,
            value=value(dual.variable),
        )
        for dual in sos_zero_duals(L, sos)
    ]
end

function _sos_form_max_abs_coefficient(form)
    entries = collect(form)
    isempty(entries) && return 0.0
    return Float64(maximum(pair -> abs(pair.second), entries))
end

function _sos_zero_dual_diagnostic_record(zero)
    label = _sos_origin_field(zero.origin, :label, nothing)
    return (
        origin=zero.origin,
        label=label,
        feature=_sos_label_field(label, :feature, nothing),
        decomposition=_sos_label_field(label, :decomposition, nothing),
        reason=_sos_label_field(label, :reason, nothing),
        coefficient_domain=_sos_label_field(label, :coefficient_domain, nothing),
        exact_coefficient_domain=_sos_label_field(label, :exact_coefficient_domain, nothing),
        kind=zero.kind,
        value=zero.value,
        abs_value=Float64(abs(zero.value)),
        term_count=length(zero.form),
        max_abs_coefficient=_sos_form_max_abs_coefficient(zero.form),
    )
end

"""
    sos_zero_dual_diagnostics(zero_values)
    sos_zero_dual_diagnostics(moment_data, sos)

Return per-row diagnostics for solved zero/equality dual variables, preserving
the row origin, coefficient-domain metadata, linear-form size, and multiplier
magnitude.
"""
function sos_zero_dual_diagnostics(zero_values::AbstractVector)
    return [
        _sos_zero_dual_diagnostic_record(zero)
        for zero in zero_values
    ]
end

function sos_zero_dual_diagnostics(mp::MomentProblem, sos::SOSProblem)
    return sos_zero_dual_diagnostics(mp.linear, sos)
end

function sos_zero_dual_diagnostics(L::MomentLinearData, sos::SOSProblem)
    return sos_zero_dual_diagnostics(sos_zero_dual_values(L, sos))
end

"""
    sos_dual_certificate_diagnostics(moment_data, sos)

Return the solved SOS dual certificate diagnostics in one record: coefficient
residuals, PSD block diagnostics, and zero/equality dual diagnostics.
This is a numerical helper for small-N checks and certificate extraction
plumbing; it is not a rigorous interval certificate.
"""
function sos_dual_certificate_diagnostics(mp::MomentProblem, sos::SOSProblem)
    return sos_dual_certificate_diagnostics(mp.linear, sos)
end

function sos_dual_certificate_diagnostics(L::MomentLinearData, sos::SOSProblem)
    residual = sos_dual_certificate_residual(L, sos)
    psd_blocks = sos_dual_block_diagnostics(L, sos)
    zero_duals = sos_zero_dual_diagnostics(L, sos)
    return (
        residual=residual,
        psd_blocks=psd_blocks,
        zero_duals=zero_duals,
        moment_count=residual.moment_count,
        psd_block_count=length(psd_blocks),
        zero_dual_count=length(zero_duals),
        max_abs_residual=residual.max_abs_residual,
        identity_residual=residual.identity_residual,
        max_residual_moment=residual.max_residual_moment,
        max_residual_value=residual.max_residual_value,
    )
end

function _assert_sos_dual_metadata_matches(mp::MomentProblem, sos::SOSProblem)
    return _assert_sos_dual_metadata_matches(mp.linear, sos)
end

function _assert_sos_dual_metadata_matches(L::MomentLinearData, sos::SOSProblem)
    length(sos.psd_dual_blocks) == length(L.psd_blocks_lin) || throw(ArgumentError(
        "SOS problem has $(length(sos.psd_dual_blocks)) PSD dual blocks but moment problem has $(length(L.psd_blocks_lin))."
    ))
    length(sos.zero_duals) == length(L.zero_constraints) || throw(ArgumentError(
        "SOS problem has $(length(sos.zero_duals)) zero duals but moment problem has $(length(L.zero_constraints))."
    ))
    return nothing
end

function _sos_residual_result(L::MomentLinearData, residuals::AbstractVector)
    residual_by_moment = Dict(L.moments[i] => residuals[i] for i in eachindex(L.moments))
    if isempty(residuals)
        max_abs_residual = 0.0
        max_residual_moment = nothing
        max_residual_value = nothing
    else
        max_idx = argmax(abs.(residuals))
        max_abs_residual = abs(residuals[max_idx])
        max_residual_moment = L.moments[max_idx]
        max_residual_value = residuals[max_idx]
    end
    return (
        moment_count=length(L.moments),
        identity_moment=L.identity,
        max_abs_residual=max_abs_residual,
        max_residual_moment=max_residual_moment,
        max_residual_value=max_residual_value,
        identity_residual=residual_by_moment[L.identity],
        residual_by_moment=residual_by_moment,
    )
end

"""
    sos_dual_certificate_residual(moment_data, sos)

Evaluate the coefficient residual of a solved SOS dual certificate:
`objective - A'(dual_psd) - B'(dual_zero) - objective_value(sos.model)`.
The returned `max_abs_residual` is a small-N numerical certificate check, not
a rigorous interval certificate. `moment_data` may be either a symbolic
`MomentProblem` or finalized `MomentLinearData`.
"""
function sos_dual_certificate_residual(
    mp::MomentProblem{A,T,M,P},
    sos::SOSProblem,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P<:Polynomial{A,T}}
    _assert_sos_dual_metadata_matches(mp, sos)
    if _is_real_moment_problem(mp)
        return _sos_dual_certificate_residual_real(mp.linear, sos)
    elseif _is_complex_problem(A)
        return _sos_dual_certificate_residual_hermitian(mp.linear, sos)
    else
        return _sos_dual_certificate_residual_real(mp.linear, sos)
    end
end

function sos_dual_certificate_residual(L::MomentLinearData, sos::SOSProblem)
    _assert_sos_dual_metadata_matches(L, sos)
    if _sos_coeff_type(L) <: Real
        return _sos_dual_certificate_residual_real(L, sos)
    else
        return _sos_dual_certificate_residual_hermitian(L, sos)
    end
end


# =============================================================================
# Coefficient Extraction for Constraints
# =============================================================================

"""
    get_Cαj(unsymmetrized_basis::Vector{M}, localizing_mtx::Matrix{P}) where {T,M,P<:AbstractPolynomial{T}}

Extract coefficient matrix for dualization from a symbolic polynomial constraint matrix.

# Arguments
- `unsymmetrized_basis`: Sorted vector of monomials as the basis
- `localizing_mtx`: Polynomial-valued constraint matrix

# Returns
- Dictionary mapping (basis_idx, row, col) to coefficient values

This function decomposes the constraint matrix into a sparse representation
where each non-zero coefficient is indexed by the basis monomial position
and the matrix position (row, col).

Note: Monomials not found in the basis (e.g., odd-parity fermionic monomials
filtered by superselection) are skipped since they have expectation value 0.
"""
function get_Cαj(unsymmetrized_basis::Vector{M}, localizing_mtx::Matrix{P}) where {T,M,P<:AbstractPolynomial{T}}
    dim = size(localizing_mtx, 1)
    cis = CartesianIndices((dim, dim))

    # basis idx, row, col -> coefficient
    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()

    n_basis = length(unsymmetrized_basis)

    for ci in cis
        for (coeff, alpha) in terms(localizing_mtx[ci])
            # `terms(::Polynomial)` yields `(coefficient, Monomial)` pairs; internal SOS basis
            # for moment problems is still in `NormalMonomial` space.
            alpha_key = alpha
            idx = searchsortedfirst(unsymmetrized_basis, alpha_key)
            # Skip monomials not in basis (e.g., odd-parity fermionic monomials)
            # These have expectation value 0 and don't contribute
            if idx <= n_basis && unsymmetrized_basis[idx] == alpha_key
                dictionary_of_keys[(idx, ci.I[1], ci.I[2])] = coeff
            end
        end
    end

    return dictionary_of_keys
end


# =============================================================================
# SOS Dualization (Unified)
# =============================================================================

"""
    sos_dualize(mp::MomentProblem{A,T,M,P}) where {A,T,M,P} -> SOSProblem

Convert a symbolic moment problem into its dual SOS (Sum of Squares) problem.

# Arguments
- `mp::MomentProblem{A,T,M,P}`: The symbolic primal moment problem to dualize

# Keyword Arguments
- `hermitian_representation::Symbol=:real_lift`: use `:real_lift` for a real
  PSD embedding or `:native` for native Hermitian PSD cones.

# Returns
- `SOSProblem`: The dual SOS problem with matrix variables and constraints

# Description
The dualization process involves:
1. Creating dual variables for each constraint:
   - `:Zero` constraints -> equality multipliers
   - `:PSD` constraints -> PSDCone (real positive semidefinite)
   - `:HPSD` constraints -> lifted real PSDCone of size `2n × 2n`
2. Building coefficient expressions for `objective - Aᴴ(dual)`
3. Maximizing the affine identity coefficient directly
4. Constraining all remaining coefficient expressions to zero

For complex algebras (Pauli, Fermionic, Bosonic), Hermitian PSD constraints
are embedded as real 2n x 2n PSD constraints using the standard construction:
```
H in HPSD <=> [Re(H), -Im(H); Im(H), Re(H)] in PSD
```

# Examples
```julia
mp = moment_relax(pop, corr_sparsity, term_sparsities)
sos = sos_dualize(mp)
set_optimizer(sos.model, Clarabel.Optimizer)
optimize!(sos.model)
obj = objective_value(sos.model)
```

See also: [`MomentProblem`](@ref), [`moment_relax`](@ref), [`SOSProblem`](@ref)
"""
function sos_dualize(
    mp::MomentProblem{A,TI,M,P};
    hermitian_representation::Symbol=:real_lift,
    coefficient_scaling::Symbol=:none,
    coefficient_scaling_floor::Real=0,
) where {A<:AlgebraType, TI<:Integer, M<:NormalMonomial{A,TI}, P<:Polynomial{A,TI}}
    if _is_real_moment_problem(mp)
        return _sos_dualize_real(mp; coefficient_scaling, coefficient_scaling_floor)
    elseif _is_complex_problem(A)
        return _sos_dualize_hermitian(
            mp;
            hermitian_representation,
            coefficient_scaling,
            coefficient_scaling_floor,
        )
    else
        return _sos_dualize_real(mp; coefficient_scaling, coefficient_scaling_floor)
    end
end

"""
    sos_dualize(L::MomentLinearData) -> SOSProblem

Convert finalized linear moment data into its dual SOS problem without requiring
a symbolic `MomentProblem` wrapper.
"""
function sos_dualize(
    L::MomentLinearData;
    hermitian_representation::Symbol=:real_lift,
    coefficient_scaling::Symbol=:none,
    coefficient_scaling_floor::Real=0,
)
    if _sos_coeff_type(L) <: Real
        return _sos_dualize_real(L; coefficient_scaling, coefficient_scaling_floor)
    else
        return _sos_dualize_hermitian(
            L;
            hermitian_representation,
            coefficient_scaling,
            coefficient_scaling_floor,
        )
    end
end


# =============================================================================
# SOS Dualization Helpers for MomentLinearData
# =============================================================================

@inline _sos_coeff_type(::MomentLinearData{K,C,M}) where {K,C,M} = C
@inline _sos_real_type(::MomentLinearData{K,C,M}) where {K,C,M} = typeof(real(one(C)))

function _sos_moment_index(L::MomentLinearData{K}, key::K) where {K}
    return _get_key_value(L.moment_index, key, "moment index")
end

function _check_real_sos_cones!(mp::MomentProblem)
    for (cone, _) in mp.constraints
        (cone == :Zero || cone == :PSD) || error("Unexpected cone type $cone for real SOS dualization")
    end
    return _check_real_sos_cones!(mp.linear)
end

function _check_real_sos_cones!(L::MomentLinearData)
    for block in L.psd_blocks_lin
        block.meta.cone == :PSD || error("Unexpected cached cone type $(block.meta.cone) for real SOS dualization")
    end
    return nothing
end

function _check_hermitian_sos_cones!(mp::MomentProblem)
    for (cone, mat) in mp.constraints
        if cone == :Zero
            _is_hermitian_poly_matrix(mat) || throw(ArgumentError(
                "Complex SOS dualization requires Hermitian zero-cone matrices. " *
                "Pass moment problems built through `moment_relax`, which splits " *
                "non-Hermitian equalities into Hermitian real/imaginary components first."
            ))
        elseif cone != :HPSD
            error("Unexpected cone type $cone for complex SOS dualization")
        end
    end
    return _check_hermitian_sos_cones!(mp.linear)
end

function _check_hermitian_sos_cones!(L::MomentLinearData)
    for block in L.psd_blocks_lin
        block.meta.cone == :HPSD || error("Unexpected cached cone type $(block.meta.cone) for complex SOS dualization")
    end
    return nothing
end

function _accumulate_dual_contribution!(
    eqs::AbstractVector,
    idx::Integer,
    coef,
    dual_block,
    row::Integer,
    col::Integer,
    cone::Symbol,
)
    cone == :PSD || error("Real SOS dualization expected :PSD block, got $cone")
    add_to_expression!(eqs[idx], -coef, dual_block[row, col])
    return nothing
end

function _normalize_hermitian_sos_representation(representation::Symbol)
    if representation in (:real_lift, :lift, :lifted, :hermitian_lift)
        return :hermitian_lift
    elseif representation in (:native, :hermitian)
        return :hermitian
    end
    throw(ArgumentError(
        "Unsupported hermitian_representation $(repr(representation)); expected :real_lift or :native."
    ))
end

function _normalize_sos_equation_scaling(scaling::Symbol)
    scaling in (:none, :max_abs) && return scaling
    throw(ArgumentError(
        "Unsupported coefficient_scaling $(repr(scaling)); expected :none or :max_abs."
    ))
end

function _normalize_sos_equation_scaling_floor(scaling_floor::Real)
    floor = Float64(scaling_floor)
    (isfinite(floor) && floor >= 0) || throw(ArgumentError(
        "coefficient_scaling_floor must be finite and nonnegative, got $scaling_floor."
    ))
    return floor
end

function _scale_sos_coefficient_expression(
    expression,
    scaling::Symbol;
    scaling_floor::Real=0,
)
    mode = _normalize_sos_equation_scaling(scaling)
    mode == :none && return expression
    floor = _normalize_sos_equation_scaling_floor(scaling_floor)

    scale = abs(JuMP.constant(expression))
    for (coefficient, _) in JuMP.linear_terms(expression)
        scale = max(scale, abs(coefficient))
    end
    scale = max(scale, floor)
    (iszero(scale) || isone(scale)) && return expression
    return expression / scale
end

function _accumulate_native_hermitian_dual_contribution!(
    eqs_re::AbstractVector,
    eqs_im::AbstractVector,
    idx::Integer,
    coef,
    dual_block,
    row::Integer,
    col::Integer,
    cone::Symbol,
)
    cone == :HPSD || error("Hermitian SOS dualization expected :HPSD block, got $cone")

    X = dual_block[row, col]
    X_re = real(X)
    X_im = imag(X)
    c_re = real(coef)
    c_im = imag(coef)

    add_to_expression!(eqs_re[idx], -c_re, X_re)
    add_to_expression!(eqs_re[idx], +c_im, X_im)
    add_to_expression!(eqs_im[idx], -c_im, X_re)
    add_to_expression!(eqs_im[idx], -c_re, X_im)
    return nothing
end

"""
    _accumulate_dual_contribution!(eqs_re, eqs_im, idx, coef, lifted, row, col, :HPSD, dim)

Accumulate one cached linear-form term into the Hermitian SOS coefficient
matching equations. This is the single home for the Hermitian real-lift adjoint
scaling: for a lifted dual matrix `Z`, the effective complex multiplier is
`X₁ + im*X₂` with `X₁ = Z₁₁ + Z₂₂` and `X₂ = Z₂₁ - Z₁₂`. The `Z₁₁ + Z₂₂`
sum is intentional; dropping either block is the classic factor-of-2 bug.
"""
function _accumulate_dual_contribution!(
    eqs_re::AbstractVector,
    eqs_im::AbstractVector,
    idx::Integer,
    coef,
    lifted,
    row::Integer,
    col::Integer,
    cone::Symbol,
    dim::Integer,
)
    cone == :HPSD || error("Hermitian SOS dualization expected :HPSD block, got $cone")

    n = Int(dim)
    i = Int(row)
    j = Int(col)
    c_re = real(coef)
    c_im = imag(coef)

    X1 = lifted[i, j] + lifted[n + i, n + j]
    X2 = lifted[n + i, j] - lifted[i, n + j]

    # Coefficient expressions are objective - Aᴴ(dual); the real identity
    # coefficient is maximized directly and the rest are constrained to zero.
    # Real(c * (X1 + im*X2)) = c_re*X1 - c_im*X2
    # Imag(c * (X1 + im*X2)) = c_im*X1 + c_re*X2
    add_to_expression!(eqs_re[idx], -c_re, X1)
    add_to_expression!(eqs_re[idx], +c_im, X2)
    add_to_expression!(eqs_im[idx], -c_im, X1)
    add_to_expression!(eqs_im[idx], -c_re, X2)
    return nothing
end


# =============================================================================
# Real SOS Dualization
# =============================================================================

"""
    _sos_dualize_real(mp::MomentProblem{A,TI,M,P}) -> SOSProblem

Internal: Dualize a real-valued symbolic moment problem.

For real algebras (NonCommutative, Projector, Unipotent), constraints use
standard PSD cones without complex embedding.
"""
function _sos_dualize_real(
    mp::MomentProblem{A,TI,M,P};
    coefficient_scaling::Symbol=:none,
    coefficient_scaling_floor::Real=0,
) where {A<:AlgebraType, TI<:Integer, M<:NormalMonomial{A,TI}, P<:Polynomial{A,TI}}
    _check_real_sos_cones!(mp)

    return _sos_dualize_real(mp.linear; coefficient_scaling, coefficient_scaling_floor)
end

function _sos_dualize_real(
    L::MomentLinearData;
    coefficient_scaling::Symbol=:none,
    coefficient_scaling_floor::Real=0,
)
    _check_real_sos_cones!(L)
    scaling = _normalize_sos_equation_scaling(coefficient_scaling)
    scaling_floor = _normalize_sos_equation_scaling_floor(coefficient_scaling_floor)

    C = _sos_coeff_type(L)
    dual_model = GenericModel{C}()

    psd_duals = Any[]
    for block in L.psd_blocks_lin
        push!(
            psd_duals,
            @variable(dual_model, [1:block.size, 1:block.size] in PSDCone(), set_string_name=false),
        )
    end
    zero_duals = [@variable(dual_model, set_string_name=false) for _ in L.zero_constraints]

    fα_constraints = [zero(GenericAffExpr{C,VariableRef}) for _ in L.moments]

    for (key, coef) in L.objective_lin
        add_to_expression!(fα_constraints[_sos_moment_index(L, key)], coef)
    end

    for (block_idx, block) in enumerate(L.psd_blocks_lin)
        dual_block = psd_duals[block_idx]
        for i in 1:block.size, j in 1:block.size
            for (key, coef) in block.entries[i, j]
                _accumulate_dual_contribution!(
                    fα_constraints,
                    _sos_moment_index(L, key),
                    coef,
                    dual_block,
                    i,
                    j,
                    block.meta.cone,
                )
            end
        end
    end

    for (zc_idx, zc) in enumerate(L.zero_constraints)
        λ = zero_duals[zc_idx]
        for (key, coef) in zc.form
            add_to_expression!(fα_constraints[_sos_moment_index(L, key)], -coef, λ)
        end
    end

    # Eliminate the usual scalar bound variable: the real identity coefficient is b.
    identity_idx = _sos_moment_index(L, L.identity)
    @objective(dual_model, Max, fα_constraints[identity_idx])

    coefficient_indices = [
        i for i in eachindex(fα_constraints) if i != identity_idx && !iszero(fα_constraints[i])
    ]
    coefficient_expressions = [
        _scale_sos_coefficient_expression(
            fα_constraints[i],
            scaling;
            scaling_floor,
        )
        for i in coefficient_indices
    ]
    isempty(coefficient_expressions) ||
        @constraint(dual_model, coefficient_expressions .== zero(C))

    psd_dual_blocks = [
        SOSDualBlock(block.meta, dual_block, :psd)
        for (block, dual_block) in zip(L.psd_blocks_lin, psd_duals)
    ]
    return SOSProblem(dual_model, length(L.moments), psd_dual_blocks, zero_duals)
end

function _sos_dual_certificate_residual_real(mp::MomentProblem, sos::SOSProblem)
    return _sos_dual_certificate_residual_real(mp.linear, sos)
end

function _sos_dual_certificate_residual_real(L::MomentLinearData, sos::SOSProblem)
    residuals = zeros(Float64, length(L.moments))

    for (key, coef) in L.objective_lin
        residuals[_sos_moment_index(L, key)] += Float64(real(coef))
    end

    for (block_idx, block) in enumerate(L.psd_blocks_lin)
        dual_block = sos.psd_dual_blocks[block_idx]
        dual_block.representation == :psd || throw(ArgumentError(
            "Expected real PSD dual block $block_idx to use representation :psd, got $(repr(dual_block.representation))."
        ))
        dual_values = value.(dual_block.variable)
        for i in 1:block.size, j in 1:block.size
            dual_value = Float64(real(dual_values[i, j]))
            for (key, coef) in block.entries[i, j]
                residuals[_sos_moment_index(L, key)] -= Float64(real(coef)) * dual_value
            end
        end
    end

    for (zc_idx, zc) in enumerate(L.zero_constraints)
        λ = Float64(real(value(sos.zero_duals[zc_idx])))
        for (key, coef) in zc.form
            residuals[_sos_moment_index(L, key)] -= Float64(real(coef)) * λ
        end
    end

    residuals[_sos_moment_index(L, L.identity)] -= Float64(objective_value(sos.model))
    return _sos_residual_result(L, residuals)
end


# =============================================================================
# Hermitian SOS Dualization
# =============================================================================

"""
    sos_dualize(mp::StateMomentProblem{A,ST,TI,M,P}) -> SOSProblem

Convert a symbolic state moment problem into its dual SOS problem.

State polynomial optimization problems (with Unipotent, Projector algebras) use
real PSD constraints without complex embedding.

# Arguments
- `mp::StateMomentProblem`: The symbolic state moment problem

# Returns
- `SOSProblem`: The dual SOS problem
"""
function sos_dualize(mp::StateMomentProblem{A,ST,TI,M,P}) where {A<:AlgebraType, ST<:StateType, TI<:Integer, M<:NCStateWord{ST,A,TI}, P<:NCStatePolynomial}
    # State polynomial optimization uses real PSD constraints
    return _sos_dualize_state(mp)
end

"""
    _sos_dualize_state(mp::StateMomentProblem) -> SOSProblem

Internal: Dualize a state polynomial moment problem.

State polynomial optimization with algebras like UnipotentAlgebra and ProjectorAlgebra
uses real-valued SDP with standard PSD constraints.
"""
function _sos_dualize_state(mp::StateMomentProblem{A,ST,TI,M,P}) where {A<:AlgebraType, ST<:StateType, TI<:Integer, M<:NCStateWord{ST,A,TI}, P<:NCStatePolynomial}
    # Get coefficient type from polynomial
    C = eltype(coefficients(mp.objective))

    dual_model = GenericModel{C}()

    # Create matrix variables for each constraint (now includes block_basis in tuple)
    dual_variables = map(mp.constraints) do (cone, mat, _block_basis)
        G_dim = size(mat, 1)
        if cone == :Zero
            @variable(dual_model, [1:G_dim, 1:G_dim] in SymmetricMatrixSpace())
        else  # :PSD
            @variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
        end
    end

    # For state polynomial optimization, we work with expectation values (StateWords)
    # The key insight: expval(<I>*xy) = <xy> = expval(<xy>*I)
    # So we convert NCStateWords to StateWords via expval for comparison

    # Create canonical StateWord basis from total_basis
    # Apply expval to convert NCStateWord to StateWord, then symmetric_canon
    SW = StateWord{ST,A,TI}
    state_basis = _sorted_stateword_basis_from_ncsw(mp.total_basis)
    n_basis = length(state_basis)
    sw_to_idx = Dict(sw => i for (i, sw) in enumerate(state_basis))

    identity_sw = one(SW)
    identity_idx = searchsortedfirst(state_basis, identity_sw)
    if identity_idx > n_basis || state_basis[identity_idx] != identity_sw
        error("Identity StateWord not found in basis - this shouldn't happen")
    end

    # Initialize constraint expressions
    fα_constraints = [zero(GenericAffExpr{C,VariableRef}) for _ in 1:n_basis]

    # Add objective polynomial coefficients
    # Objective NCStateWords have form <xy>*I, we need to match them to basis StateWords
    missing_objective_words = SW[]
    for (coef, ncsw) in zip(coefficients(mp.objective), monomials(mp.objective))
        canon_sw = symmetric_canon(expval(ncsw))
        sym_idx = get(sw_to_idx, canon_sw, 0)
        if iszero(sym_idx)
            push!(missing_objective_words, canon_sw)
            continue
        end
        add_to_expression!(fα_constraints[sym_idx], coef)
    end
    _throw_missing_state_words(missing_objective_words, "the objective"; source="Relaxation basis")

    # Process each constraint matrix by iterating directly over block (row, col) pairs
    # This ensures all StateWord contributions are accumulated correctly across blocks
    # (matching NCTSSOS's on-the-fly loop approach)
    missing_constraint_words = SW[]
    for (i, (_, mat, _block_basis)) in enumerate(mp.constraints)
        dim = size(mat, 1)
        # Iterate over upper triangle to avoid double-counting
        for row in 1:dim
            for col in row:dim
                poly = mat[row, col]
                for (coeff, ncsw) in zip(coefficients(poly), monomials(poly))
                    # Convert NCStateWord to canonical StateWord
                    sw = symmetric_canon(expval(ncsw))
                    sw_idx = get(sw_to_idx, sw, 0)
                    if iszero(sw_idx)
                        push!(missing_constraint_words, sw)
                        continue
                    end
                    # Multiply off-diagonal entries by 2 for symmetric matrix contribution
                    effective_coeff = (row == col) ? coeff : 2 * coeff
                    add_to_expression!(fα_constraints[sw_idx], -effective_coeff, dual_variables[i][row, col])
                end
            end
        end
    end
    _throw_missing_state_words(missing_constraint_words, "the constraint matrices"; source="Relaxation basis")

    # Eliminate the usual scalar bound variable: the identity coefficient is b.
    @objective(dual_model, Max, fα_constraints[identity_idx])

    coefficient_indices = [
        i for i in eachindex(fα_constraints) if i != identity_idx && !iszero(fα_constraints[i])
    ]
    isempty(coefficient_indices) || @constraint(dual_model, fα_constraints[coefficient_indices] .== 0)

    return SOSProblem(dual_model, n_basis)
end

"""
    _sos_dualize_hermitian(mp::MomentProblem{A,TI,M,P}) -> SOSProblem

Internal: Dualize a complex/Hermitian symbolic moment problem using the cached
`mp.linear::MomentLinearData` forms.

For complex algebras (Pauli, Fermionic, Bosonic), each cached `:HPSD` block gets
a lifted real PSD dual matrix. Cached scalar zero constraints get real free
multipliers. The Hermitian factor-of-2 convention lives only in
`_accumulate_dual_contribution!`.
"""
function _sos_dualize_hermitian(
    mp::MomentProblem{A,TI,M,P};
    hermitian_representation::Symbol=:real_lift,
    coefficient_scaling::Symbol=:none,
    coefficient_scaling_floor::Real=0,
) where {A<:AlgebraType, TI<:Integer, M<:NormalMonomial{A,TI}, P<:Polynomial{A,TI}}
    _check_hermitian_sos_cones!(mp)

    return _sos_dualize_hermitian(
        mp.linear;
        hermitian_representation,
        coefficient_scaling,
        coefficient_scaling_floor,
    )
end

function _sos_dualize_hermitian(
    L::MomentLinearData;
    hermitian_representation::Symbol=:real_lift,
    coefficient_scaling::Symbol=:none,
    coefficient_scaling_floor::Real=0,
)
    representation = _normalize_hermitian_sos_representation(hermitian_representation)
    if representation == :hermitian
        return _sos_dualize_hermitian_native(
            L;
            coefficient_scaling,
            coefficient_scaling_floor,
        )
    end
    return _sos_dualize_hermitian_lift(
        L;
        coefficient_scaling,
        coefficient_scaling_floor,
    )
end

function _sos_dualize_hermitian_lift(
    L::MomentLinearData;
    coefficient_scaling::Symbol=:none,
    coefficient_scaling_floor::Real=0,
)
    _check_hermitian_sos_cones!(L)
    scaling = _normalize_sos_equation_scaling(coefficient_scaling)
    scaling_floor = _normalize_sos_equation_scaling_floor(coefficient_scaling_floor)

    RC = _sos_real_type(L)
    dual_model = GenericModel{RC}()

    psd_duals = Any[]
    for block in L.psd_blocks_lin
        dim = block.size
        push!(psd_duals, (
            dim = dim,
            lifted = @variable(
                dual_model,
                [1:2*dim, 1:2*dim] in PSDCone(),
                set_string_name=false,
            ),
        ))
    end
    zero_duals = [@variable(dual_model, set_string_name=false) for _ in L.zero_constraints]

    fα_constraints_re = [zero(GenericAffExpr{RC,VariableRef}) for _ in L.moments]
    fα_constraints_im = [zero(GenericAffExpr{RC,VariableRef}) for _ in L.moments]

    for (key, coef) in L.objective_lin
        idx = _sos_moment_index(L, key)
        add_to_expression!(fα_constraints_re[idx], real(coef))
        add_to_expression!(fα_constraints_im[idx], imag(coef))
    end

    for (block_idx, block) in enumerate(L.psd_blocks_lin)
        dual = psd_duals[block_idx]
        for i in 1:block.size, j in 1:block.size
            for (key, coef) in block.entries[i, j]
                _accumulate_dual_contribution!(
                    fα_constraints_re,
                    fα_constraints_im,
                    _sos_moment_index(L, key),
                    coef,
                    dual.lifted,
                    i,
                    j,
                    block.meta.cone,
                    dual.dim,
                )
            end
        end
    end

    for (zc_idx, zc) in enumerate(L.zero_constraints)
        λ = zero_duals[zc_idx]
        for (key, coef) in zc.form
            idx = _sos_moment_index(L, key)
            add_to_expression!(fα_constraints_re[idx], -real(coef), λ)
            add_to_expression!(fα_constraints_im[idx], -imag(coef), λ)
        end
    end

    # Eliminate only the real identity equation; the imaginary identity equation
    # still enforces a real scalar bound.
    identity_idx = _sos_moment_index(L, L.identity)
    @objective(dual_model, Max, fα_constraints_re[identity_idx])

    real_coefficient_indices = [
        i for i in eachindex(fα_constraints_re) if i != identity_idx && !iszero(fα_constraints_re[i])
    ]
    imag_coefficient_indices = [
        i for i in eachindex(fα_constraints_im) if !iszero(fα_constraints_im[i])
    ]
    real_coefficient_expressions = [
        _scale_sos_coefficient_expression(
            fα_constraints_re[i],
            scaling;
            scaling_floor,
        )
        for i in real_coefficient_indices
    ]
    imag_coefficient_expressions = [
        _scale_sos_coefficient_expression(
            fα_constraints_im[i],
            scaling;
            scaling_floor,
        )
        for i in imag_coefficient_indices
    ]
    isempty(real_coefficient_expressions) ||
        @constraint(dual_model, real_coefficient_expressions .== zero(RC))
    isempty(imag_coefficient_expressions) ||
        @constraint(dual_model, imag_coefficient_expressions .== zero(RC))

    psd_dual_blocks = [
        SOSDualBlock(block.meta, dual.lifted, :hermitian_lift)
        for (block, dual) in zip(L.psd_blocks_lin, psd_duals)
    ]
    return SOSProblem(dual_model, length(L.moments), psd_dual_blocks, zero_duals)
end

function _sos_dualize_hermitian_native(
    L::MomentLinearData;
    coefficient_scaling::Symbol=:none,
    coefficient_scaling_floor::Real=0,
)
    _check_hermitian_sos_cones!(L)
    scaling = _normalize_sos_equation_scaling(coefficient_scaling)
    scaling_floor = _normalize_sos_equation_scaling_floor(coefficient_scaling_floor)

    RC = _sos_real_type(L)
    dual_model = GenericModel{RC}()

    psd_duals = Any[]
    for block in L.psd_blocks_lin
        push!(
            psd_duals,
            @variable(
                dual_model,
                [1:block.size, 1:block.size] in HermitianPSDCone(),
                set_string_name=false,
            ),
        )
    end
    zero_duals = [@variable(dual_model, set_string_name=false) for _ in L.zero_constraints]

    fα_constraints_re = [zero(GenericAffExpr{RC,VariableRef}) for _ in L.moments]
    fα_constraints_im = [zero(GenericAffExpr{RC,VariableRef}) for _ in L.moments]

    for (key, coef) in L.objective_lin
        idx = _sos_moment_index(L, key)
        add_to_expression!(fα_constraints_re[idx], real(coef))
        add_to_expression!(fα_constraints_im[idx], imag(coef))
    end

    for (block_idx, block) in enumerate(L.psd_blocks_lin)
        dual_block = psd_duals[block_idx]
        for i in 1:block.size, j in 1:block.size
            for (key, coef) in block.entries[i, j]
                _accumulate_native_hermitian_dual_contribution!(
                    fα_constraints_re,
                    fα_constraints_im,
                    _sos_moment_index(L, key),
                    coef,
                    dual_block,
                    i,
                    j,
                    block.meta.cone,
                )
            end
        end
    end

    for (zc_idx, zc) in enumerate(L.zero_constraints)
        λ = zero_duals[zc_idx]
        for (key, coef) in zc.form
            idx = _sos_moment_index(L, key)
            add_to_expression!(fα_constraints_re[idx], -real(coef), λ)
            add_to_expression!(fα_constraints_im[idx], -imag(coef), λ)
        end
    end

    identity_idx = _sos_moment_index(L, L.identity)
    @objective(dual_model, Max, fα_constraints_re[identity_idx])

    real_coefficient_indices = [
        i for i in eachindex(fα_constraints_re) if i != identity_idx && !iszero(fα_constraints_re[i])
    ]
    imag_coefficient_indices = [
        i for i in eachindex(fα_constraints_im) if !iszero(fα_constraints_im[i])
    ]
    real_coefficient_expressions = [
        _scale_sos_coefficient_expression(
            fα_constraints_re[i],
            scaling;
            scaling_floor,
        )
        for i in real_coefficient_indices
    ]
    imag_coefficient_expressions = [
        _scale_sos_coefficient_expression(
            fα_constraints_im[i],
            scaling;
            scaling_floor,
        )
        for i in imag_coefficient_indices
    ]
    isempty(real_coefficient_expressions) ||
        @constraint(dual_model, real_coefficient_expressions .== zero(RC))
    isempty(imag_coefficient_expressions) ||
        @constraint(dual_model, imag_coefficient_expressions .== zero(RC))

    psd_dual_blocks = [
        SOSDualBlock(block.meta, dual_block, :hermitian)
        for (block, dual_block) in zip(L.psd_blocks_lin, psd_duals)
    ]
    return SOSProblem(dual_model, length(L.moments), psd_dual_blocks, zero_duals)
end

function _sos_dual_certificate_residual_hermitian(mp::MomentProblem, sos::SOSProblem)
    return _sos_dual_certificate_residual_hermitian(mp.linear, sos)
end

function _sos_dual_certificate_residual_hermitian(L::MomentLinearData, sos::SOSProblem)
    residuals_re = zeros(Float64, length(L.moments))
    residuals_im = zeros(Float64, length(L.moments))

    for (key, coef) in L.objective_lin
        idx = _sos_moment_index(L, key)
        residuals_re[idx] += Float64(real(coef))
        residuals_im[idx] += Float64(imag(coef))
    end

    for (block_idx, block) in enumerate(L.psd_blocks_lin)
        dual_block = sos.psd_dual_blocks[block_idx]
        dual_block.representation in (:hermitian_lift, :hermitian) || throw(ArgumentError(
            "Expected Hermitian PSD dual block $block_idx to use representation :hermitian_lift or :hermitian, got $(repr(dual_block.representation))."
        ))
        native_values = _sos_native_dual_block_value(dual_block, value.(dual_block.variable))
        for i in 1:block.size, j in 1:block.size
            multiplier = native_values[i, j]
            for (key, coef) in block.entries[i, j]
                idx = _sos_moment_index(L, key)
                contribution = coef * multiplier
                residuals_re[idx] -= Float64(real(contribution))
                residuals_im[idx] -= Float64(imag(contribution))
            end
        end
    end

    for (zc_idx, zc) in enumerate(L.zero_constraints)
        λ = Float64(real(value(sos.zero_duals[zc_idx])))
        for (key, coef) in zc.form
            idx = _sos_moment_index(L, key)
            residuals_re[idx] -= Float64(real(coef)) * λ
            residuals_im[idx] -= Float64(imag(coef)) * λ
        end
    end

    residuals_re[_sos_moment_index(L, L.identity)] -= Float64(objective_value(sos.model))
    return _sos_residual_result(L, complex.(residuals_re, residuals_im))
end
