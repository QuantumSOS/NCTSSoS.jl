module NCTSSoSManiDSDPExt

using LinearAlgebra
using NCTSSoS
using ManiDSDP

const MOI = NCTSSoS.MOI

const _MANIDSDP_DEP_TOL = 1e-9

struct _DiagonalGramFactor{R<:Real,V<:AbstractVector{R}}
    invdiag::V
end

Base.:\(factor::_DiagonalGramFactor, x::AbstractVector) = factor.invdiag .* x
Base.:\(factor::_DiagonalGramFactor, X::AbstractMatrix) = reshape(factor.invdiag, :, 1) .* X

struct _SparseAccumulator
    rows::Vector{Int}
    cols::Vector{Int}
    vals::Vector{ComplexF64}
end

_SparseAccumulator() = _SparseAccumulator(Int[], Int[], ComplexF64[])

mutable struct _PauliCoordinateDualSDP{T<:Number,
                                       R<:Real,
                                       GT<:AbstractMatrix{R},
                                       FT} <: ManiDSDP.AbstractDualSDP{T}
    C::Matrix{T}
    b::Vector{R}
    rank::Int
    gram::GT
    gram_factor::FT
    D::Matrix{T}
    rows::Vector{Int}
    cols::Vector{Int}
    vals::Vector{T}
    starts::Vector{Int}
    n::Int
    D_dot_C::Float64
end

function ManiDSDP.apply_A(problem::_PauliCoordinateDualSDP, X::AbstractMatrix)
    size(X) == size(problem.C) || throw(DimensionMismatch("X must have size $(size(problem.C))"))
    m = length(problem.b)
    out = Vector{Float64}(undef, m)
    @inbounds for k in 1:m
        acc = zero(ComplexF64)
        for ptr in problem.starts[k]:(problem.starts[k + 1] - 1)
            acc += conj(problem.vals[ptr]) * X[problem.rows[ptr], problem.cols[ptr]]
        end
        out[k] = real(acc)
    end
    return out
end

function ManiDSDP.apply_At(problem::_PauliCoordinateDualSDP{T}, y::AbstractVector) where {T}
    length(y) == length(problem.b) || throw(DimensionMismatch("length(y) must match length(b)"))
    out = zeros(T, problem.n, problem.n)
    @inbounds for k in 1:length(problem.b)
        yk = Float64(ManiDSDP._real_value(y[k], "y[$k]"))
        iszero(yk) && continue
        for ptr in problem.starts[k]:(problem.starts[k + 1] - 1)
            out[problem.rows[ptr], problem.cols[ptr]] += yk * problem.vals[ptr]
        end
    end
    return ManiDSDP._symmetrize(out)
end

function ManiDSDP.project_range(problem::_PauliCoordinateDualSDP, X::AbstractMatrix)
    coeffs = problem.gram_factor \ ManiDSDP.apply_A(problem, X)
    return ManiDSDP.apply_At(problem, coeffs)
end

function ManiDSDP.dual_y(problem::_PauliCoordinateDualSDP, S::AbstractMatrix, Xtilde::AbstractMatrix, sigma::Real)
    size(S) == size(problem.C) || throw(DimensionMismatch("S must have size $(size(problem.C))"))
    inv_sigma = inv(sigma)
    rhs = Vector{Float64}(undef, length(problem.b))
    @inbounds for k in eachindex(problem.b)
        acc = zero(ComplexF64)
        for ptr in problem.starts[k]:(problem.starts[k + 1] - 1)
            row = problem.rows[ptr]
            col = problem.cols[ptr]
            acc += conj(problem.vals[ptr]) * (S[row, col] + problem.C[row, col] + inv_sigma * Xtilde[row, col])
        end
        rhs[k] = real(acc)
    end
    return problem.gram_factor \ rhs
end

function ManiDSDP.residual(problem::_PauliCoordinateDualSDP{T}, S::AbstractMatrix, y::AbstractVector) where {T}
    length(y) == length(problem.b) || throw(DimensionMismatch("length(y) must match length(b)"))
    out = Matrix{T}(undef, problem.n, problem.n)
    @inbounds for idx in eachindex(out)
        out[idx] = -S[idx] - problem.C[idx]
    end
    @inbounds for k in eachindex(problem.b)
        yk = Float64(ManiDSDP._real_value(y[k], "y[$k]"))
        iszero(yk) && continue
        for ptr in problem.starts[k]:(problem.starts[k + 1] - 1)
            out[problem.rows[ptr], problem.cols[ptr]] += yk * problem.vals[ptr]
        end
    end
    return out
end

function ManiDSDP.gradient_matrix(problem::_PauliCoordinateDualSDP{T},
                                  S::AbstractMatrix,
                                  Xtilde::AbstractMatrix,
                                  sigma::Real) where {T}
    y = ManiDSDP.dual_y(problem, S, Xtilde, sigma)
    R = ManiDSDP.residual(problem, S, y)
    G = Matrix{T}(undef, problem.n, problem.n)
    @inbounds for idx in eachindex(G)
        G[idx] = Xtilde[idx] + problem.D[idx] - sigma * R[idx]
    end
    return G, y, R
end

function ManiDSDP.update_multiplier!(::_PauliCoordinateDualSDP,
                                     Xtilde::AbstractMatrix,
                                     R::AbstractMatrix,
                                     sigma::Real)
    @inbounds for idx in eachindex(Xtilde, R)
        Xtilde[idx] -= sigma * R[idx]
    end
    return Xtilde
end

function ManiDSDP.augmented_cost(problem::_PauliCoordinateDualSDP,
                                 Y::AbstractMatrix,
                                 Xtilde::AbstractMatrix,
                                 sigma::Real)
    S = ManiDSDP.slack_from_factor(Y)
    y = ManiDSDP.dual_y(problem, S, Xtilde, sigma)
    R = ManiDSDP.residual(problem, S, y)
    return ManiDSDP._real_inner(problem.D, S) +
           problem.D_dot_C -
           ManiDSDP._real_inner(Xtilde, R) +
           (sigma / 2) * ManiDSDP._real_inner(R, R)
end

# The extension is intentionally narrow: it bypasses JuMP only for Pauli
# SOS-dual relaxations whose PSD blocks are plain moment-matrix blocks.  Extra
# equality/localizing machinery changes the dual shape and must not be shoved
# through ManiDSDP by accident.
function NCTSSoS.solve_sdp(
    mp::NCTSSoS.MomentProblem{A},
    optimizer::Union{Type{ManiDSDP.Optimizer},ManiDSDP.Optimizer};
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
) where {A<:NCTSSoS.PauliAlgebra}
    opt = optimizer isa ManiDSDP.Optimizer ? optimizer : optimizer()
    return _solve_pauli_with_manidsdp(
        mp,
        opt;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
    )
end

function NCTSSoS.solve_sdp(
    mp::NCTSSoS.MomentProblem{A},
    optimizer::MOI.OptimizerWithAttributes;
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
) where {A<:NCTSSoS.PauliAlgebra}
    opt = _manidsdp_optimizer_from_attributes(optimizer)
    return _dispatch_or_fallback(
        mp,
        optimizer,
        opt;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
    )
end

function NCTSSoS.solve_sdp(
    mp::NCTSSoS.MomentProblem{A},
    optimizer::F;
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
) where {A<:NCTSSoS.PauliAlgebra,F<:Function}
    opt = _try_construct_manidsdp_optimizer(optimizer)
    return _dispatch_or_fallback(
        mp,
        optimizer,
        opt;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
    )
end

function _dispatch_or_fallback(
    mp,
    optimizer,
    opt;
    dualize::Bool,
    formulation::Symbol,
    representation::Symbol,
    orphan_policy::Symbol,
)
    if opt === nothing
        return invoke(
            NCTSSoS.solve_sdp,
            Tuple{Any,Any},
            mp,
            optimizer;
            dualize=dualize,
            formulation=formulation,
            representation=representation,
            orphan_policy=orphan_policy,
        )
    end
    return _solve_pauli_with_manidsdp(
        mp,
        opt;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
    )
end

function _solve_pauli_with_manidsdp(
    mp,
    opt::ManiDSDP.Optimizer;
    dualize::Bool,
    formulation::Symbol,
    representation::Symbol,
    orphan_policy::Symbol,
)
    dualize || throw(ArgumentError("ManiDSDP native Pauli path supports only dualize=true."))
    if formulation != :moment_variables || representation != :real || orphan_policy != :error
        throw(ArgumentError("Moment lowering options apply only with dualize=false; ManiDSDP native Pauli path solves the SOS dual directly."))
    end

    problem, constant = _pauli_moment_problem_to_manidsdp(mp, opt)
    objective, status = if problem === nothing
        constant, MOI.OPTIMAL
    else
        result = ManiDSDP.solve(problem; settings=opt.settings, opt.solve_kwargs...)
        constant + result.objective, _manidsdp_status(result.status)
    end

    return (
        objective = objective,
        model = NCTSSoS.GenericModel{typeof(objective)}(),
        n_unique_elements = mp.n_unique_moment_matrix_elements,
        status = status,
    )
end

function _manidsdp_optimizer_from_attributes(spec::MOI.OptimizerWithAttributes)
    opt = _try_construct_manidsdp_optimizer(spec.optimizer_constructor)
    opt === nothing && return nothing
    for (attr, value) in spec.params
        _apply_manidsdp_attribute!(opt, attr, value)
    end
    return opt
end

function _try_construct_manidsdp_optimizer(constructor)
    if constructor isa Type && constructor <: ManiDSDP.Optimizer
        return constructor()
    end
    opt = try
        constructor()
    catch
        return nothing
    end
    return opt isa ManiDSDP.Optimizer ? opt : nothing
end

function _apply_manidsdp_attribute!(opt::ManiDSDP.Optimizer, attr, value)
    if attr isa MOI.RawOptimizerAttribute
        MOI.set(opt, attr, value)
    elseif attr isa MOI.Silent
        Bool(value) && MOI.set(opt, MOI.RawOptimizerAttribute("verbose"), false)
    elseif attr isa AbstractString
        MOI.set(opt, MOI.RawOptimizerAttribute(attr), value)
    elseif attr isa Symbol
        MOI.set(opt, MOI.RawOptimizerAttribute(String(attr)), value)
    else
        throw(ArgumentError("unsupported ManiDSDP optimizer attribute $(attr) in NCTSSoS extension"))
    end
    return opt
end

function _pauli_moment_problem_to_manidsdp(mp, opt::ManiDSDP.Optimizer)
    L = mp.linear
    _check_supported_pauli_manidsdp_shape(L)

    blocks = L.psd_blocks_lin
    offsets = _block_offsets(blocks)
    n = isempty(blocks) ? 0 : offsets[end] + blocks[end].size - 1
    n > 0 || throw(ArgumentError("Pauli ManiDSDP path requires at least one HPSD moment block."))

    moment_count = length(L.moments)
    raw = Vector{Union{Nothing,_SparseAccumulator}}(nothing, moment_count)
    for (block_idx, block) in pairs(blocks)
        offset = offsets[block_idx]
        for j in 1:block.size, i in 1:block.size
            row = offset + i - 1
            col = offset + j - 1
            for (key, coef) in block.entries[i, j]
                idx = L.moment_index[key]
                acc = raw[idx]
                if acc === nothing
                    acc = _SparseAccumulator()
                    raw[idx] = acc
                end
                push!(acc.rows, row)
                push!(acc.cols, col)
                push!(acc.vals, ComplexF64(coef))
            end
        end
    end

    objective = zeros(ComplexF64, moment_count)
    for (key, coef) in L.objective_lin
        objective[L.moment_index[key]] += ComplexF64(coef)
    end

    identity = L.identity
    identity_idx = L.moment_index[identity]
    B0_rows, B0_cols, B0_vals, _ = _hermitian_coordinates(raw[identity_idx], n, one(ComplexF64))
    B0 = _dense_from_coordinates(B0_rows, B0_cols, B0_vals, n)
    constant = real(objective[identity_idx])

    _assert_unit_diagonal_manidsdp_objective(B0)

    rows = Int[]
    cols = Int[]
    vals = ComplexF64[]
    starts = Int[1]
    gram_diag = Float64[]
    rhs = Float64[]

    # The imaginary identity equation is normally zero for Pauli Hamiltonians,
    # but keep it if present and harmless: it is an affine equality, not the
    # ManiDSDP objective matrix.
    _append_coordinate_constraint!(rows, cols, vals, starts, gram_diag, rhs,
                                   raw[identity_idx], n, ComplexF64(im), imag(objective[identity_idx]))

    for idx in eachindex(L.moments)
        idx == identity_idx && continue
        C_key_raw = raw[idx]
        c_key = objective[idx]
        if C_key_raw === nothing
            abs(c_key) <= _MANIDSDP_DEP_TOL || throw(ArgumentError(
                "SOS objective references a moment with no moment-matrix coefficient; native ManiDSDP problem is infeasible."
            ))
            continue
        end
        _append_coordinate_constraint!(rows, cols, vals, starts, gram_diag, rhs,
                                       C_key_raw, n, one(ComplexF64), real(c_key))
        _append_coordinate_constraint!(rows, cols, vals, starts, gram_diag, rhs,
                                       C_key_raw, n, ComplexF64(im), imag(c_key))
    end

    isempty(rhs) && return nothing, constant

    rank = Int(get(opt.solve_kwargs, :rank, min(n, 2)))
    problem = _coordinate_pauli_dualsdp(-B0, rows, cols, vals, starts, gram_diag, rhs; rank=rank)
    return problem, constant
end

function _check_supported_pauli_manidsdp_shape(L)
    isempty(L.zero_constraints) || throw(ArgumentError(
        "ManiDSDP native Pauli path currently supports pure moment-matrix SOS relaxations only; " *
        "found $(length(L.zero_constraints)) scalar zero/equality constraints."
    ))
    isempty(L.psd_blocks_lin) && throw(ArgumentError("ManiDSDP native Pauli path found no PSD blocks."))
    for (idx, block) in pairs(L.psd_blocks_lin)
        block.meta.cone == :HPSD || throw(ArgumentError(
            "ManiDSDP native Pauli path expected HPSD block $idx, got $(block.meta.cone)."
        ))
        block.meta.origin isa NCTSSoS.MomentMatrixOrigin || throw(ArgumentError(
            "ManiDSDP native Pauli path supports moment-matrix blocks only; " *
            "block $idx has origin $(typeof(block.meta.origin))."
        ))
    end
    return nothing
end

function _block_offsets(blocks)
    offsets = Vector{Int}(undef, length(blocks))
    next = 1
    for (idx, block) in pairs(blocks)
        offsets[idx] = next
        next += block.size
    end
    return offsets
end

function _accumulate_coordinate!(dict::Dict{Int,ComplexF64}, key::Int, value::ComplexF64)
    if haskey(dict, key)
        dict[key] += value
    else
        dict[key] = value
    end
    return dict
end

function _hermitian_coordinates(acc::Union{Nothing,_SparseAccumulator}, n::Integer, phase::ComplexF64)
    dict = Dict{Int,ComplexF64}()
    if acc !== nothing
        half = 0.5
        phase_conj = conj(phase)
        @inbounds for idx in eachindex(acc.vals)
            row = acc.rows[idx]
            col = acc.cols[idx]
            value = acc.vals[idx]
            _accumulate_coordinate!(dict, row + (col - 1) * n, half * phase * conj(value))
            _accumulate_coordinate!(dict, col + (row - 1) * n, half * phase_conj * value)
        end
    end

    keys_sorted = sort!(collect(keys(dict)))
    rows = Int[]
    cols = Int[]
    vals = ComplexF64[]
    sizehint!(rows, length(keys_sorted))
    sizehint!(cols, length(keys_sorted))
    sizehint!(vals, length(keys_sorted))
    norm2 = 0.0
    for key in keys_sorted
        value = dict[key]
        abs(value) <= _MANIDSDP_DEP_TOL && continue
        col, row0 = divrem(key - 1, n)
        push!(rows, row0 + 1)
        push!(cols, col + 1)
        push!(vals, value)
        norm2 += abs2(value)
    end
    return rows, cols, vals, norm2
end

function _dense_from_coordinates(rows::Vector{Int}, cols::Vector{Int}, vals::Vector{ComplexF64}, n::Integer)
    out = zeros(ComplexF64, n, n)
    @inbounds for idx in eachindex(vals)
        out[rows[idx], cols[idx]] += vals[idx]
    end
    return Matrix(Hermitian(out))
end

function _append_coordinate_constraint!(rows, cols, vals, starts, gram_diag, rhs,
                                        acc::Union{Nothing,_SparseAccumulator},
                                        n::Integer,
                                        phase::ComplexF64,
                                        b::Real)
    c_rows, c_cols, c_vals, norm2 = _hermitian_coordinates(acc, n, phase)
    norm2 <= _MANIDSDP_DEP_TOL^2 && abs(b) <= _MANIDSDP_DEP_TOL && return false
    norm2 <= _MANIDSDP_DEP_TOL^2 && throw(ArgumentError(
        "SOS coefficient equation has zero matrix but nonzero rhs $(b); the ManiDSDP native problem would be infeasible."
    ))
    @inbounds for idx in eachindex(c_vals)
        if c_rows[idx] == c_cols[idx] && abs(c_vals[idx]) > 1e-8
            throw(ArgumentError(
                "ManiDSDP requires non-identity SOS coefficient matrices to have zero diagonal; " *
                "diagonal entry $(c_rows[idx]) = $(c_vals[idx])."
            ))
        end
        push!(rows, c_rows[idx])
        push!(cols, c_cols[idx])
        push!(vals, -c_vals[idx])
    end
    push!(starts, length(vals) + 1)
    push!(gram_diag, norm2)
    push!(rhs, -Float64(b))
    return true
end

function _coordinate_pauli_dualsdp(C::AbstractMatrix,
                                   rows::Vector{Int},
                                   cols::Vector{Int},
                                   vals::Vector{ComplexF64},
                                   starts::Vector{Int},
                                   gram_diag::Vector{Float64},
                                   b::Vector{Float64};
                                   rank::Integer)
    n = size(C, 1)
    1 <= rank <= n || throw(ArgumentError("rank must be in 1:size(C, 1)"))
    length(starts) == length(b) + 1 || throw(ArgumentError("starts must have length length(b) + 1"))
    length(gram_diag) == length(b) || throw(ArgumentError("gram_diag must match length(b)"))
    all(>(_MANIDSDP_DEP_TOL), gram_diag) || throw(ArgumentError(
        "ManiDSDP Pauli SOS coefficient matrices must be numerically nonzero."
    ))

    gram = Diagonal(gram_diag)
    factor = _DiagonalGramFactor(inv.(gram_diag))
    problem = _PauliCoordinateDualSDP(
        Matrix{ComplexF64}(C),
        b,
        Int(rank),
        gram,
        factor,
        zeros(ComplexF64, n, n),
        rows,
        cols,
        vals,
        starts,
        n,
        0.0,
    )
    problem.D .= ManiDSDP.apply_At(problem, factor \ b)
    problem.D_dot_C = ManiDSDP._real_inner(problem.D, problem.C)
    return problem
end

function _assert_unit_diagonal_manidsdp_objective(B0::AbstractMatrix)
    for i in axes(B0, 1)
        abs(real(B0[i, i]) - 1.0) <= 1e-8 && abs(imag(B0[i, i])) <= 1e-8 ||
            throw(ArgumentError("ManiDSDP requires the SOS identity matrix to have unit real diagonal; entry $i is $(B0[i, i])."))
    end
    return nothing
end

function _manidsdp_status(status::Symbol)
    status === :optimal && return MOI.OPTIMAL
    status === :max_iter && return MOI.ITERATION_LIMIT
    return MOI.OTHER_ERROR
end

end # module
